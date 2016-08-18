configfile: "config.yaml"

targets = []
for comparison in config['comparisons']:
    targets.extend(
        expand(
            'output/{comparison}/{bait}/{bait}_stats.txt',
            comparison=comparison,
            bait=config['comparisons'][comparison]
        )
    )
    
HERE = srcdir('')
WRAPPER = './wrappers'

def wrapper_for(tool):
    """
    Returns the wrapper directory for the specified tool
    """
    return os.path.join(HERE, WRAPPER, tool)  

rule all:
    input:
        targets
        
rule reduced_genome:
    input:
        enzyme = "reduced_genome/{enzyme}.fa",
        genome = "reduced_genome/{genome}.fa"
    output:
        "reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique.fa",
        "reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique_2.fa",
        "reduced_genome/{genome}_{enzyme}_flanking_sites_{fragLen}_unique_2.bed"
    params:
        fragLen = "{fragLen}",
        enzyme = "{enzyme}",
        genome = "{genome}"
    wrapper:
        wrapper_for('reduced_genome')

        
rule bowtie2_build:
    input:
        reduced_genome="reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique.fa"
    output:
        "reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique.{num}.bt2"
    params:
        prefix = "reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique"
    shell:
        "bowtie2-build {input.reduced_genome} {params.prefix}"

genome=None
enzyme=None
fragLen=None

def get_bowtie(wc):
     genome=config['baits'][wc.bait]['reduced_genome']
     enzyme=config['baits'][wc.bait]['primary_enz_name']
     fragLen=config['baits'][wc.bait]['fragment_len']
     num=[1, 2, 3, 4]
     return expand("reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique.{num}.bt2")
        
rule bowtie2:
    input:
        bowtie=get_bowtie,
        fastq=lambda wc: config["samples"][wc.bait][wc.sample]
    output:
        aligned_sam="sam_files/{wildcards.bait}/{sample}_aligned.sam",
        unaligned_sam="sam_files/{wildcards.bait}/{sample}_unaligned.sam"
    threads:
        8
    params:
        index="reduced_genome/%s_%s_flanking_sequences_%s_unique" % (genome, enzyme, fragLen)
    shell:
        "bowtie2 -p {threads} -N 0 -5 {fragLen} \
        --un {output.unaligned_sam} \
        -x {params.index} \
        -U {input.fastq} \
        -S {output.aligned_sam}"
   

rule bedGraph_Counts:
    input:
        sam="sam_files/{wildcards.bait}/{sample}_aligned.sam"
    output:
        "sam_files/{wildcards.bait}/{sample}_aligned.bedGraph"
    wrapper:
        wrapper_for("sam_bedGraph")

def get_fasta(wc):
    genome=config['baits'][wc.bait]['reduced_genome']
    enzyme=config['baits'][wc.bait]['primary_enz_name']
    fragLen=config['baits'][wc.bait]['fragment_len']
    return expand('reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique_2.fa')

def get_bed(wc):
    genome=config['baits'][wc.bait]['reduced_genome']
    enzyme=config['baits'][wc.bait]['primary_enz_name']
    fragLen=config['baits'][wc.bait]['fragment_len']
    return expand('reduced_genome/{genome}_{enzyme}_flanking_sites_{fragLen}_unique_2.bed')

rule purify_helper:
    input:
        fasta=get_fasta,
        bed=get_bed
    output:
        temp("temp.bed")
    params:
        extra=lambda wc: config["baits"][wc.bait]["primer"]
    wrapper:
        wrapper_for("4C_purify")


rule purify_aligned_bedGraph:
    input:
        aligned_bedG="sam_files/{wildcards.bait}/{sample}_aligned.bedGraph",
        temp_bed = "temp.bed"
    output:
        "bedGraphs/{wildcards.bait}/{sample}_aligned_rm_self_und.bedGraph"
    shell:
        "bedtools substract -a input.aligned_bedG -b input.temp_bed"

def get_bedgraphs(wc):
    samples = config['comparisons'][wc.comparison][wc.bait]
    return expand('bedGraphs/{wildcards.bait}/{sample}_aligned_rm_self_und.bedGraph', sample=samples)

        
rule R_script:
    input:
        bedGraph=get_bedgraphs
    output:
        "output/{comparison}/{wildcards.bait}/{wildcards.bait}_stats.txt"
    shell:
        "Rscript 4C.R --bait {wildcards.bait} --comparison {wildcards.comparison}"
