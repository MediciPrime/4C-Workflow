configfile: "config.yaml"

targets = []
for comparison in config['comparisons']:
    targets.extend(
        expand(
            'output/{comparison}/{bait}_stats.txt',
            comparison=comparison,
            bait=config['comparisons'][comparison]['baits']
        )
    )

HERE = srcdir('')
WRAPPER = './wrappers'

def wrapper_for(tool):
    """
    Returns the wrapper directory for the specified tool
    """
    return os.path.join(HERE, WRAPPER, tool)  

"""
primer= config[comparison]["primer"]
genome= config[comparison]["reduced_genome"]
enzyme= config[comparison]["primary_enz_name"]
fragLen= config[comparison]["fragment_len"]
"""

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

rule bowtie2:
    input:
        bowtie=expand("reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique.{num}.bt2",
                      genome=config[comparison]["reduced_genome"],
                      enzyme=config[comparison]["primary_enz_name"],
                      fragLen=config[comparison]["fragment_len"],
                      num=[1, 2, 3, 4]),
        fastq=lambda wc: config["samples"][wc.sample]
    output:
        aligned_sam="sam_files/{sample}_aligned.sam",
        unaligned_sam="sam_files/{sample}_unaligned.sam"
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
        sam="sam_files/{sample}_aligned.sam"
    output:
        "sam_files/{sample}_aligned.bedGraph"
    wrapper:
        wrapper_for("sam_bedGraph")


rule purify_helper:
    input:
        fasta="reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique_2.fa",
        bed="reduced_genome/{genome}_{enzyme}_flanking_sites_{fragLen}_unique_2.bed"
    output:
        temp("temp.bed")
    params:
        extra=primer
    wrapper:
        wrapper_for("4C_purify")

rule purify_aligned_bedGraph:
    input:
        aligned_bedG="sam_files/{sample}_aligned.bedGraph",
        temp_bed = "temp.bed"
    output:
        "bedGraphs/{sample}_aligned_rm_self_und.bedGraph"
    run:
        import pybedtools
        bt = pybedtools.BedTool(input.aligned_bedG)
        bt.subtract(input.temp_bed, output=output[0])

def get_bedgraphs(wc):
    samples = config['comparisons'][wc.comparison][wc.bait]
    return expand('bedGraphs/{wildcards.bait}/{sample}_aligned_rm_self_und.bedGraph', sample=samples)
        
rule R_script:
    input:
        bedGraph=get_bedgraphs
    output:
        "output/{comparison}/{bait}_stats.txt"
    shell:
        "Rscript 4C.R {input.bedGraph} {wildcards.bait}"
