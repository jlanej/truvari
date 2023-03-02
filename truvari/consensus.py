"""
Utilities to produce haplotypes from sets of variants
"""
import logging
from io import StringIO
import pysam


def create_allele(ref, variants, refpos = 0):
    """
    Create haplotypes from phased variants

    ref is the reference sequence
    Assumes all variants are contained within the ref region

    returns the sequence produced
    """
    seq = StringIO()
    last_pos = 0
    for entry in variants:
        seq.write(ref[last_pos: entry.start - refpos])
        seq.write(entry.alts[0])
        last_pos = entry.stop - refpos
    seq.write(ref[last_pos:])
    seq.seek(0)
    return seq.read()


def make_paths(variants, refseq, refstart, max_paths=2):
    """
    Given a set of variants, find upto max_paths consensus sequences through the variants
    """
    m_paths = [[] for _ in range(max_paths)]
    cur_path = 0
    for var in variants:
        is_placed = False
        # if homozygous, needs to be applied to both...
        while not is_placed and cur_path < max_paths:
            # empty path or non-overlapping - can place
            if not m_paths[cur_path] or var.start > m_paths[cur_path][-1].stop:
                m_paths[cur_path].append(var)
                is_placed = True
            cur_path += 1
        #cur_path %= max_paths # try to separate the variants evenly?
        cur_path = 0
        if not is_placed:
            logging.error("Cannot make non-overlapping paths for variants. Try increasing max_paths")
    return [create_allele(refseq, _, refstart) for _ in m_paths]


def test():
    chrom, start, end = "chr20", 35539135, 35539719
    # Buffer..
    ref = pysam.FastaFile("/Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa")
    variants = pysam.VariantFile("chr20.HG002.strawman.vcf.gz")
    
    m_vars = variants.fetch(chrom, start, end) # just assume, whatever
    BUFFER = 100
    start -= BUFFER
    end += BUFFER
    m_seq = ref.fetch(chrom, start, end)
    seqs = make_paths(m_vars, m_seq, start)
    for pos, i in enumerate(seqs):
        pos += 1
        print(f">samp_{pos}_\n{i}")
    print(f">ref_{chrom}_{start + 1}_{end}\n{m_seq}")

    
if __name__ == '__main__':
    test()
