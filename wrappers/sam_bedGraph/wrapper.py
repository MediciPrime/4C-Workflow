__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

shell("""
awk '{{print $3}}' {snakemake.input.sam} \
        | grep "^ch" | sort | uniq -c | sed -e 's/_\|:\| \|-/\\t/g' \
        |  awk '{{OFS="\t";print $2,$3,$4,$1}}' > \
        {snakemake.output}"""
)
