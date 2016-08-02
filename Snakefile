configfile: "config.yaml"


targets=expand("reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique.{num}.bt2",
                num=[1, 2, 3, 4],
                genome=config["other"]["reduced_genome"],
                enzyme=config["other"]["primary_enz_name"],
                fragLen=config["other"]["fragment_len"])
targets+=expand("sam_files/{sample}_aligned.sam",
                sample=config["samples"].keys())
"""
targets+=expand("sam_files/{sample}_aligned.bedGraph",
                sample=config["samples"]), "temp.bed"
"""

targets+=expand("bedGraphs/{sample}_aligned_rm_self_und.bedGraph",
                sample=config["samples"])

targets+=expand("output/{comparison}/{sample}_nearbait_norm_counts.bedGraph",
                comparison=config["comparisons"],
                sample=config["samples"])



HERE = srcdir('')
WRAPPER = './wrappers'

def wrapper_for(tool):
    """
    Returns the wrapper directory for the specified tool
    """
    return os.path.join(HERE, WRAPPER, tool)  

primer= config["other"]["primer"]
genome= config["other"]["reduced_genome"]
enzyme= config["other"]["primary_enz_name"]
fragLen= config["other"]["fragment_len"]

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
                      num=[1, 2, 3, 4], genome=genome, enzyme=enzyme, fragLen=fragLen),
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
        fasta="reduced_genome/{genome}_{enzyme}_flanking_sequences_{fragLen}_unique_2.fa"\
                     .format(**locals()),
        bed="reduced_genome/{genome}_{enzyme}_flanking_sites_{fragLen}_unique_2.bed"\
                   .format(**locals())
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

rule R_script:
    input:
        sampleTable="{comparison}.tsv",
        bedGraph="bedGraphs/{samples}_aligned_rm_self_und.bedGraph"
    output:
        "output/{comparison}/{samples}_nearbait_norm_counts.bedGraph"
    shell:
        "Rscript 4C.R {input.sampleTable}"
