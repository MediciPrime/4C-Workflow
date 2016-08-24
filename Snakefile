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
for bait in config['baits']:
    targets.extend(
        expand(
            'reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique.fa',
            bait=bait,
            enzyme=config['baits'][bait]['primary_enz_name'],
            fragLen=config['baits'][bait]['fragment_len']))
    
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
        genome = config['reference_genome']['location']
    output:
        'reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique.fa',
        "reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique_2.fa",
        "reduced_genome/{enzyme}_{bait}_flanking_sites_{fragLen}_unique_2.bed"
    params:
        fragLen = "{fragLen}",
        enzyme = "{enzyme}",
        genome = config['reference_genome']['name'],
        bait = "{bait}"
    wrapper:
        wrapper_for('reduced_genome')
        
rule bowtie2_build:
    input:
        reduced_genome="reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique.fa"
    output:
        "reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique.{num}.bt2"
    params:
        prefix = "reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique"
    shell:
        "bowtie2-build -q {input.reduced_genome} {params.prefix}"

def get_bt2(wc):
    enzyme=config['baits'][wc.bait]['primary_enz_name']
    fragLen=config['baits'][wc.bait]['fragment_len']
    return expand("reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique.{num}.bt2",
                  enzyme=enzyme, fragLen=fragLen, num=[1, 2, 3, 4],
                  bait=wc.bait)

def index(wc):
    enzyme=config['baits'][wc.bait]['primary_enz_name']
    fragLen=config['baits'][wc.bait]['fragment_len']
    return (
        "reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique"
        .format(enzyme=enzyme, fragLen=fragLen, bait=wc.bait))

rule bowtie2:
    input:
        bowtie=get_bt2,
        fastq=lambda wc: config["samples"][bait][wc.sample]
    output:
        aligned_sam="sam_files/{bait}/{sample}_aligned.sam",
        unaligned_sam="sam_files/{bait}/{sample}_unaligned.sam"
    threads:
        8
    params:
        index=index,
        fragLen=config['baits'][bait]['fragment_len']
    shell:
        "bowtie2 -p {threads} -N 0 -5 {params.fragLen} "
        "--un {output.unaligned_sam} "
        "-x {params.index} "
        "-U {input.fastq} "
        "-S {output.aligned_sam}"
   

rule bedGraph_Counts:
    input:
        sam="sam_files/{bait}/{sample}_aligned.sam"
    output:
        "sam_files/{bait}/{sample}_aligned.bedGraph"
    wrapper:
        wrapper_for("sam_bedGraph")

def get_fasta(wc):
    enzyme=config['baits'][wc.bait]['primary_enz_name']
    fragLen=config['baits'][wc.bait]['fragment_len']
    return expand('reduced_genome/{enzyme}_{bait}_flanking_sequences_{fragLen}_unique_2.fa',
                  enzyme=enzyme, fragLen=fragLen, bait=wc.bait)

def get_bed(wc):
    enzyme=config['baits'][wc.bait]['primary_enz_name']
    fragLen=config['baits'][wc.bait]['fragment_len']
    return expand('reduced_genome/{enzyme}_{bait}_flanking_sites_{fragLen}_unique_2.bed',
                  enzyme=enzyme, fragLen=fragLen, bait=wc.bait)
    
rule purify_helper:
    input:
        fasta=get_fasta,
        bed=get_bed
    output:
        temp("temp_{bait}.bed")
    params:
        extra=config['baits'][bait]['primer']
    wrapper:
        wrapper_for("4C_purify")


rule purify_aligned_bedGraph:
    input:
        aligned_bedG="sam_files/{bait}/{sample}_aligned.bedGraph",
        temp_bed = "temp_{bait}.bed"
    output:
        "bedGraphs/{bait}/{sample}_aligned_rm_self_und.bedGraph"
    shell:
        "bedtools subtract -a {input.aligned_bedG} -b {input.temp_bed} "
        "> {output[0]} "


def get_bedgraphs(wc):
    samples = []
    for v in config['comparisons'][wc.comparison][wc.bait].values():
        samples.extend(v)
    return expand('bedGraphs/{bait}/{sample}_aligned_rm_self_und.bedGraph',
                  sample=samples, bait=wc.bait)


rule R_script:
    input:
        bedGraph=get_bedgraphs
    output:
        "output/{comparison}/{bait}/{bait}_stats.txt"
    shell:
        "Rscript 4C.R --bait {wildcards.bait} --comparison {wildcards.comparison}"
