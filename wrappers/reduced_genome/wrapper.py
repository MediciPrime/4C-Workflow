__author__ = "Behram Radmanesh"
__copyright__ = "Copyright 2016, Behram Radmanesh"
__email__ = "behram.radmanesh@nih.gov"
__license__ = "MIT"

# include any imports that will be needed, example below
from snakemake.shell import shell

# Use this block to support arbitrary arguments to be passed in to the shell
# call without raising an error if no params.extra were provided. If there
# aren't any additional params possible then you can delete this.
try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

# Use this block to redirect stdout and stderr to a log if it was provided.
# Adjust the stdout/stderr redirection as needed.
if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""

# Formatting: put logical blocks of arguments together on one line, but use
# multiple lines for clarity. Note the space at the end of each line for when
# the lines are concatenated. Also note the use of {extra} and {log} as
# specified above.
shell("""
fl={snakemake.params.fragLen}
enzyme={snakemake.params.enzyme}
genome={snakemake.params.genome}
enzymeLoc={snakemake.input.enzyme}
genomeLoc={snakemake.input.genome}
oligoMatch ${{enzymeLoc}} ${{genomeLoc}} \
	   ${{genome}}_${{enzyme}}_restriction_sites_oligomatch.bed
#get coordinates of upstream fragments
awk -v fl=$fl '{{print $1"\\t"$2-fl"\\t"$2}}' \
    ${{genome}}_${{enzyme}}_restriction_sites_oligomatch.bed >up.txt
#get coordinates of downstream fragments
awk -v fl=$fl '{{print $1"\\t"$3"\\t"$3+fl}}' \
    ${{genome}}_${{enzyme}}_restriction_sites_oligomatch.bed > down.txt
#combine up and downstream fragments
cat up.txt down.txt > ${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_2.bed
#remove any fragments with negative coordinates (incude as prior step)
awk '{{if($2 >= 0 && $3 >=0) print $0}}' \
    ${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_2.bed \
    | grep -v -E 'random|JH|GL|Un' - | sort \
    | uniq  > ${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_unique_2.bed
#get the sequence of unique flanking coordinates
fastaFromBed -fi ${{genomeLoc}} -bed \
	     ${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_unique_2.bed -fo \
	     ${{genome}}_${{enzyme}}_flanking_sequences_${{fl}}_unique_2.fa
#get only sequences from FASTA file
grep -v '^>' ${{genome}}_${{enzyme}}_flanking_sequences_${{fl}}_unique_2.fa \
    | sort | uniq -i -u \
    | grep -xF -f - -B 1 \
	   ${{genome}}_${{enzyme}}_flanking_sequences_${{fl}}_unique_2.fa \
    | grep -v '^--' > ${{genome}}_${{enzyme}}_flanking_sequences_${{fl}}_unique.fa

#remove unwanted intermediate files
rm up.txt
rm down.txt

#make a BED file of unique sequences
grep '^>' ${{genome}}_${{enzyme}}_flanking_sequences_${{fl}}_unique.fa > \
     ${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_unique.bed
sed -i 's/>//g' ${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_unique.bed
sed -i 's/:\|-/\\t/g' ${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_unique.bed

# move all dm6* files to reduced_genome directory
rm ./reduced_genome/${{genome}}.fa.fai
rm ./${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_2.bed
rm ./${{genome}}_${{enzyme}}_flanking_sites_${{fl}}_unique.bed
rm ./${{genome}}_${{enzyme}}_restriction_sites_oligomatch.bed
mv ./dm6* reduced_genome/

exit 0;
""")
