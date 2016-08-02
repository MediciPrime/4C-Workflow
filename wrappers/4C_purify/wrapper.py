from __future__ import print_function
from snakemake.shell import shell
from Bio import SeqIO, Seq

__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"


try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""


def primer_coords(fn, primer):
    """
    Given fasta file `fn`, extract headers of the form "chr1:1-100" into tuples
    of (chr1, 1, 100), but only for cases where the primer is found in the
    sequence.
    """
    primer_forward = primer.lower()
    primer_reverse = str(Seq.Seq(primer).reverse_complement()).lower()

    for i, rec in enumerate(SeqIO.parse(fn, 'fasta')):
        rec_seq = rec.seq.lower()
        if (primer_forward in rec_seq) or (primer_reverse in rec_seq):
            coord = rec.name
            chrom, startstop = coord.split(':')
            start, stop = startstop.split('-')
            yield i, (chrom, start, stop)


def parse_fasta(fasta, bed, primer):
    f = open('temp.bed', 'w')
    res = list(primer_coords(fasta, primer))
    if len(res) == 0:
        raise ValueError("Primer %s was not found in %s" % (primer, fasta))

    if len(res) > 1:
        raise ValueError("More than one primer found: %s" % res)
    pos, coords = res[0]
    for i, line in enumerate(open(bed)):
        if (pos - 1) <= i <= (pos + 1):
            if i == pos:
                assert tuple(line.strip().split('\t')[:3]) == coords
            f.write(line)
        continue
    f.close()


# run to purify bed
parse_fasta(
    fasta=snakemake.input.fasta,
    bed=snakemake.input.bed,
    primer=extra)


"""
parse_fasta(
    fasta = '../../reduced_genome/dm6_hindiii_flanking_sequences_37_unique_2.fa',
    bed = '../../reduced_genome/dm6_hindiii_flanking_sites_37_unique_2.bed',
    primer = 'ACCGTTTTTGACAACAGCAGCTGTA'
)
"""
