"""
Wrapper around MAFFT and other tools to perform an MSA comparison of phased variants
"""
import os
import re
import sys
import glob
import shutil
import logging
import argparse
import resource
import multiprocessing
from io import StringIO
from functools import partial

import pysam
from pysam import bcftools
import truvari

DEFAULT_MAFFT_PARAM="--auto --thread 1"

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="phab", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--region", type=str, required=True,
                        help="Bed filename or comma-separated list of chrom:start-end regions to process")
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline vcf to MSA")
    parser.add_argument("-c", "--comp", type=str, default=None,
                        help="Comparison vcf to MSA")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference")
    parser.add_argument("--buffer", type=int, default=100,
                        help="Number of reference bases before/after region to add to MSA (%(default)s)")
    parser.add_argument("-o", "--output", type=str, default="phab_out.vcf.gz",
                        help="Output VCF")
    parser.add_argument("--bSamples", type=str, default=None,
                        help="Subset of samples to MSA from base-VCF")
    parser.add_argument("--cSamples", type=str, default=None,
                        help="Subset of samples to MSA from comp-VCF")
    parser.add_argument("-m", "--mafft-params", type=str, default=DEFAULT_MAFFT_PARAM,
                        help="Parameters for mafft, wrap in a single quote (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    return args

def parse_regions(argument):
    """
    Parse the --region rgument
    returns list of regions
    """
    ret = []
    if not argument:
        return ret
    if os.path.exists(argument):
        for i in truvari.opt_gz_open(argument):
            try:
                chrom, start, end = i.strip().split('\t')[:3]
                start = int(start)
                end = int(end)
            except Exception: #pylint: disable=broad-except
                logging.error("Unable to parse bed line %s", i)
                sys.exit(1)
            ret.append((chrom, start, end))
    else:
        for i in argument.split(','):
            try:
                chrom, start, end = re.split(':|-', i)
                start = int(start)
                end = int(end)
            except ValueError:
                logging.error("Unable to parse region line %s", i)
                sys.exit(1)
            ret.append((chrom, start, end))
    return ret


def pull_variants(vcf_fn, chrom, start, end):
    """
    Given vcf and a region, grab the variants
    """
    ret = []
    start -= 100
    end += 100
    vcf = pysam.VariantFile(vcf_fn)
    for entry in vcf.fetch(chrom, start, end):
        st, ed = truvari.entry_boundaries(entry)
        if start <= st < ed < end:
            ret.append(entry)
    return ret


def build_consensus(vcf, ref, region, output, samples=None, prefix_name=False):
    """
    Make the consensus sequence - appends to output
    """
    chrom, start, end = region
    prefix = 'p:' if prefix_name else ''
    ref = output + '.ref.fa'
    cmd = "cat {ref} | bcftools consensus -H{hap} --sample {sample} --prefix {prefix}{sample}_{hap}_ {vcf} >> {output}"
    for samp in samples:
        for hap in [1, 2]:
            m_cmd = cmd.format(ref=ref,
                               chrom=chrom,
                               start=start,
                               end=end,
                               hap=hap,
                               sample=samp,
                               prefix=prefix,
                               vcf=vcf,
                               output=output)
            ret = truvari.cmd_exe(m_cmd, pipefail=True)
            if ret.ret_code != 0:
                logging.error("Unable to make haplotype %d for sample %s", hap, samp)
                logging.error(ret.stderr)
                sys.exit(1)

def run_mafft(seq_fn, output, params=DEFAULT_MAFFT_PARAM):
    """
    Run mafft - return True if successful
    """
    cmd = f"mafft {params} {seq_fn} > {output}"
    ret = truvari.cmd_exe(cmd)
    if ret.ret_code != 0:
        logging.error("Unable to run MAFFT on %s", seq_fn)
        logging.error(ret.stderr)
        return False
    return True

def make_haplotypes(variants, refseq, refstart, sample=0):
    """
    Make the phased haplotypes for a set of variants
    """
    def add_var(entry, seq, last_pos, ref, refpos):
        """
        Add upto variant's end to the sequence
        """
        seq.write(ref[last_pos: entry.start - refpos])
        seq.write(entry.alts[0])
        return entry.stop - refpos

    m_paths = [StringIO(), StringIO()]
    last_pos = [0, 0]
    for entry in variants:
        if entry.samples[sample]["GT"][0] == 1:
            last_pos[0] = add_var(entry, m_paths[0], last_pos[0], refseq, refstart)
        if entry.samples[sample]["GT"][1] == 1:
            last_pos[1] = add_var(entry, m_paths[1], last_pos[1], refseq, refstart)
    for pos, path in enumerate(m_paths):
        path.write(refseq[last_pos[pos]:])
        path.seek(0)
    return [_.read() for _ in m_paths]

def phab(var_region, base_vcf, reference, buffer=100,
        comp_vcf=None, bSamples=None, cSamples=None,
        mafft_params=DEFAULT_MAFFT_PARAM, prefix_comp=False):
    """
    Harmonize variants with MSA.

    :param `var_region`: Tuple of region's (chrom, start, end)
    :type `var_region`: :class:`tuple`
    :param `base_vcf`: VCF file name
    :type `base_vcf`: :class:`str`
    :param `reference`: Reference file name
    :type `reference`: :class:`str`
    :param `buffer`: Flanking reference bases to add to haplotypes
    :type `buffer`: :class:`int`
    :param `comp_vcf`: VCF file name
    :type `comp_vcf`: :class:`str`
    :param `bSamples`: Samples from `base_vcf` to create haplotypes
    :type `bSamples`: :class:`list`
    :param `cSamples`: Samples from `comp_vcf` to create haplotypes
    :type `cSamples`: :class:`list`
    :param `mafft_params`: Parameters for mafft
    :type `mafft_params`: :class:`str`
    :param `prefix_comp`: Ensure unique sample names by prefixing comp samples
    :type `prefix_comp`: :class:`bool`
    """
    if bSamples is None:
        bSamples = list(pysam.VariantFile(base_vcf).header.samples)
    if comp_vcf and cSamples is None:
        cSamples = list(pysam.VariantFile(comp_vcf).header.samples)
    
    b_variants = pull_variants(base_vcf, *var_region)
    c_variants = [] if comp_vcf is None else pull_variants(comp_vcf, *var_region)
    
    # No variants to harmonize?
    if len(b_variants) + len(c_variants) == 0:
        return ""

    start = var_region[1] - buffer
    end = var_region[2] + buffer
    anchor_base = 'N' # will need later
    sequences = truvari.make_temp_filename(suffix=".fa")

    ref = pysam.FastaFile(reference)
    with open(sequences, 'w') as fout:
        oseq = ref.fetch(var_region[0], start - 1, end)
        anchor_base = oseq[0]
        ref_seq = oseq[1:]
        fout.write(f">ref_{var_region[0]}_{start}_{end}\n{ref_seq}\n")
        
        for samp in bSamples:
            for i, seq in enumerate(make_haplotypes(b_variants, ref_seq, start, samp)):
                name = f"{samp}_{i + 1}_"
                fout.write(f">{name}\n{seq}\n")
                    
        if cSamples is not None:
            prefix = 'p:' if prefix_comp else ''
            for samp in cSamples:
                for i, seq in enumerate(make_haplotypes(c_variants, ref_seq, start, samp)):
                    name = f"{prefix}{samp}_{i + 1}_"
                    fout.write(f">{name}\n{seq}\n")
 
    msa_output = sequences + '.msa'
    if not run_mafft(sequences, msa_output, mafft_params):
        # hmmm... silent failures?
        return ""

    return truvari.msa2vcf(msa_output, anchor_base)

def check_requirements():
    """
    ensure external programs are in PATH
    """
    check_fail = False
    for prog in ["mafft"]:
        if not shutil.which(prog):
            logging.error("Unable to find `%s` in PATH", prog)
            check_fail = True
    return check_fail

def check_params(args):
    """
    Ensure files are okay to use
    """
    check_fail = False
    if not args.output.endswith(".vcf.gz"):
        logging.error("Output file must be a '.vcf.gz', got %s", args.output)
        check_fail = True
    if args.comp is not None and not os.path.exists(args.comp):
        logging.error("File %s does not exist", args.comp)
        check_fail = True
    if not os.path.exists(args.base):
        logging.error("File %s does not exist", args.base)
        check_fail = True
    if args.comp is not None and not args.comp.endswith(".gz"):
        logging.error(
            "Comparison vcf %s does not end with .gz. Must be bgzip'd", args.comp)
        check_fail = True
    if args.comp is not None and not os.path.exists(args.comp + '.tbi'):
        logging.error(
            "Comparison vcf index %s.tbi does not exist. Must be indexed", args.comp)
        check_fail = True
    if not args.base.endswith(".gz"):
        logging.error(
            "Base vcf %s does not end with .gz. Must be bgzip'd", args.base)
        check_fail = True
    if not os.path.exists(args.base + '.tbi'):
        logging.error(
            "Base vcf index %s.tbi does not exist. Must be indexed", args.base)
        check_fail = True
    if not os.path.exists(args.reference):
        logging.error("Reference %s does not exist", args.reference)
        check_fail = True

    return check_fail


def phab_multi(var_regions, base_vcf, reference, buf=100, comp_vcf=None,
               bSamples=None, cSamples=None, mafft_params=DEFAULT_MAFFT_PARAM,
               prefix_comp=False, threads=1):
    """
    Run phab on multiple regions
    """
    # Do a partial instead of phab_wrapper
    params = {"base_vcf": base_vcf,
              "reference": reference,
              "buffer": buf,
              "comp_vcf": comp_vcf,
              "bSamples": bSamples,
              "cSamples": cSamples,
              "mafft_params": mafft_params,
              "prefix_comp": prefix_comp}
    phab_partial = partial(phab, **params)

    result = StringIO()
    result.write('##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf = pysam.VariantFile(base_vcf)
    for ctg in vcf.header.contigs:
        result.write(str(vcf.header.contigs[ctg].header_record))
    phab_partial(var_regions[0])
    with multiprocessing.Pool(threads, maxtasksperchild=1) as pool:
        chrom_line = False
        for i in pool.imap_unordered(phab_partial, var_regions):
            if i == "": # no variants in region to phab
                continue
            # Need to write the CHROM header line once
            if not chrom_line:
                result.write(i)
                chrom_line = True
            else:
                result.write(i[i.index('\n') + 1:])
        pool.close()
        pool.join()
    result.seek(0)
    return result.read()

def phab_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    truvari.setup_logging(args.debug, show_version=True)
    if check_requirements() or check_params(args):
        logging.error("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    # Setting up the samples is tricky
    if args.bSamples is None:
        args.bSamples = list(pysam.VariantFile(args.base).header.samples)
    else:
        args.bSamples = args.bSamples.split(',')

    if args.comp:
        if args.cSamples is None:
            args.cSamples = list(pysam.VariantFile(args.comp).header.samples)
        else:
            args.cSamples = args.cSamples.split(',')

    prefix_comp = False
    if args.cSamples and set(args.cSamples) & set(args.bSamples):
        logging.warning("--cSamples intersect with --bSamples. Output vcf --comp SAMPLE names will have prefix 'p:'")
        prefix_comp = True

    all_regions = parse_regions(args.region)
    result = phab_multi(all_regions, args.base, args.reference, args.buffer, args.comp,
               args.bSamples, args.cSamples, args.mafft_params, prefix_comp, args.threads)

    with open(args.output[:-len('.gz')], 'w') as fout:
        fout.write(result)
    truvari.compress_index_vcf(args.output[:-len('.gz')], fout=args.output)
    logging.info("Finished phab")
