"""
Turn an MSA fasta into VCF. Assumes one entry is reference with name >ref_chrom:start-end
"""
from io import StringIO
from collections import defaultdict

import pysam

REFIDX = 3
ALTIDX = 4


def make_all_tests():
    """
    Every msa pattern I can think of
    ref_msa, alt_msa, variants (pos, ref, alt)
    """
    ret = []
    # Clean insertion
    ret.append(("ATCG----------ATCG",
                "ATCGATCGGATCGGATCG",
                [[3, "G", "GATCGGATCGG"]]))
    # Clean deletion
    ret.append(("ATCGATCGGATCGGATCG",
                "ATCG----------ATCG",
                [[3, "GATCGGATCGG", "G"]]))
    # SNP
    ret.append(("ATCGCG",
                "ATCACG",
                [[3, "G", "A"]]))
    # SNP on anchor base of insertion
    ret.append(("ATCG----------ATCG",
                "ATCCATCGGATCGGATCG",
                [[3, "G", "CATCGGATCGG"]]))
    # SNP on anchor base of deletion
    ret.append(("ATCGATCGGATCGGATCG",
                "ATCT----------ATCG",
                [[3, "GATCGGATCGG", "T"]]))
    # MNP
    ret.append(("ATCGCG",
                "ATAACG",
                [[2, "CG", "AA"]]))
    
    # Starting INS
    ret.append(("----------ATCG",
                "ATCGGATCGGATCG",
                [[-1, "N", "NATCGGATCGG"]]))

    # InDel - we're just going to make them overlapping
    ret.append(("ATCG-----ATCGGATCGGATCG",
                "ATCGATCGG----------ATCG",
                [[3, "GATCGGATCGG", "G"],
                 [3, "G", "GATCGG"]))

    return ret

def zip_seq_var(refseq, altseq, anchor_base='N'):
    """
    Zip the msa and create list of variants [[POS (zero based), REF, ALT]]
    """
    ret = []
    prev_mode = 0 # match - 1 mismatch - 2 insertion - 3 deletion
    cur_variant = []
    cur_pos = 0
    for ref_base, alt_base in zip(refseq, altseq):
        # Null
        if ref_base == '-' and alt_base == '-':
            continue
        
        # What are we looking at
        if ref_base == alt_base: # Mat
            cur_mode = 0
        elif ref_base == '-': # Ins
            cur_mode = 2
        elif alt_base == '-': # Del
            cur_mode = 3
        else: # Mis
            cur_mode = 1
        
        ref_base = '' if ref_base == '-' else ref_base
        alt_base = '' if alt_base == '-' else alt_base

        # Back to matching - finish previous variant
        if cur_mode == 0 and cur_mode != prev_mode and cur_variant:
            # need a 'alter with anchor base' method TODO
            ret.append(cur_variant)
            cur_variant = []
            continue
        
        # First change
        if prev_mode == 0 and cur_mode != 0:
            cur_variant = [cur_pos, anchor_base + ref_base, anchor_base + alt_base]
            prev_mode = cur_mode
            
            continue
        
        # Continue change
        # InsDel - keep extending
        if cur_mode in [2, 3]:
            cur_variant[REFIDX] += ref_base
            cur_variant[ALTIDX] += alt_base
            
        # And we'll add anchor base at the end

def msa_to_vars(msa, ref_seq, chrom, start_pos=0, abs_anchor_base='N'):
    """
    Turn MSA into VCF entries and their presence in samples
    returns list of sample names parsed and dictionary of variant : samples containing the variant
    """
    sample_names = set()

    final_vars = defaultdict(list)

    for alt_key in msa.references:
        if alt_key.startswith("ref_"):
            continue
        # gross
        cur_samp_hap = "_".join(alt_key.split('_')[:-1])
        sample_names.add("_".join(alt_key.split('_')[:-2]))
        alt_seq = msa[alt_key].upper()

        anchor_base = ref_seq[0] if ref_seq[0] != '-' else abs_anchor_base

        cur_variant = []
        cur_pos = start_pos + 1 # still have a problem here with anchor base.
        # This is too long. need to have a separate zip method
        for ref_base, alt_base in zip(ref_seq, alt_seq):
            is_ref = ref_base != '-'
            if ref_base == '-':
                ref_base = ""
            if alt_base == '-':
                alt_base = ""

            # nothing to compare
            if not ref_base and not alt_base:
                continue

            if ref_base == alt_base: # No variant
                if cur_variant and is_ref: # back to matching reference
                    # Anchor base correction
                    if len(cur_variant[REFIDX]) == len(cur_variant[ALTIDX]):
                        cur_variant[REFIDX] = cur_variant[REFIDX][1:]
                        cur_variant[ALTIDX] = cur_variant[ALTIDX][1:]
                        cur_variant[1] += 1
                    key = "\t".join([str(_) for _ in cur_variant])
                    # this is a weird edge check
                    # sometimes reference bases aren't aligned
                    if cur_variant[REFIDX] != cur_variant[ALTIDX]:
                        final_vars[key].append(cur_samp_hap)
                    cur_variant = []
            else:
                if not cur_variant:
                    # We don't need to correct cur_pos for the anchor base
                    # because we're now assuming that the ref coordinates are zero based
                    cur_variant = [chrom, cur_pos - 1, '.', anchor_base + ref_base, \
                                   anchor_base + alt_base, '.', '.', '.', 'GT']
                else:
                    cur_variant[REFIDX] += ref_base
                    cur_variant[ALTIDX] += alt_base
            if is_ref:
                cur_pos += 1
                anchor_base = ref_base
        # End Zipping
        if cur_variant:
            # Anchor base correction
            if len(cur_variant[REFIDX]) == len(cur_variant[ALTIDX]):
                cur_variant[REFIDX] = cur_variant[REFIDX][1:]
                cur_variant[ALTIDX] = cur_variant[ALTIDX][1:]
                cur_variant[1] += 1
            key = "\t".join([str(_) for _ in cur_variant])
            # this is a weird edge check
            # sometimes reference bases aren't aligned
            if cur_variant[REFIDX] != cur_variant[ALTIDX]:
                final_vars[key].append(cur_samp_hap)
        # End alignment
    sample_names = sorted(list(sample_names))
    return sample_names, final_vars

def make_vcf(variants, sample_names):
    """
    Write VCF - building GTs
    """
    out = StringIO()
    out.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\t")
    out.write("\t".join(sample_names) + '\n')

    for var in variants:
        out.write(var)
        for sample in sample_names:
            out.write('\t')
            gt = ["0", "0"]
            if sample + '_1' in variants[var]:
                gt[0] = "1"
            if sample + '_2' in variants[var]:
                gt[1] = "1"
            out.write("/".join(gt))
        out.write('\n')
    out.seek(0)
    return out.read()

def msa2vcf(fn, anchor_base=None):
    """
    Parse an MSA file and returns its VCF as a string

    Assumes one entry in the MSA has the name `ref_${chrom}:${start}-${end}` which gives VCF entries coordinates

    Provide anchor_base to prevent 'N' from being used as an anchor base
    Returns a string
    """
    msa = pysam.FastaFile(fn)

    # The ref_key identifies reference
    ref_key = [_ for _ in msa.references if _.startswith("ref_")][0] # pylint: disable=not-an-iterable
    ref_seq = msa[ref_key].upper()
    _, chrom, start_pos, _ = ref_key.split('_')
    start_pos = int(start_pos)

    sample_names, variants = msa_to_vars(msa, ref_seq, chrom, start_pos, anchor_base)

    return make_vcf(variants, sample_names)
