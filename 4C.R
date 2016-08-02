#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()
parser$add_argument("file", nargs=1, help="SampleTable to be used")

args <- parser$parse_args()

file <- args$file

if(file.access(file) == -1) {
    stop(sprintf("Specified sampletable ( %s ) does not exist", file))
} else {
    sampletable = file
}

library(R.4Cker)
library(yaml)

#st = read.table(as.character(sampletable), header=TRUE)

st = read.table("siAGO2_vs_dsLamin.tsv", header=TRUE)

setwd("bedGraphs/")

config = yaml.load_file("../config.yaml")

files = as.character(lapply(st$sampleID, function(x) paste0(x, "_aligned_rm_self_und.bedGraph")))

conditions = as.character(unique(unlist(lapply(strsplit(unlist(lapply(st$sampleID, as.character)), '_'), function(x) x[1]))))

replicates = c(as.vector(table(st$treatment)["control"]), as.vector(table(st$treatment)["treated"]))

samples = as.character(st$sampleID)

output_dir = sprintf("../output/%s_vs_%s/", conditions[1], conditions[2])

my_obj = createR4CkerObjectFromFiles(files = c(files[1], files[2], files[3], files[4]), bait_chr = shQuote(config$other$bait_chr), bait_coord = config$other$bait_coord, bait_name = config$other$bait_name, primary_enz = config$other$primary_enz, samples = c(samples[1], samples[2], samples[3], samples[4]), conditions = c(conditions[1], conditions[2]), replicates = c(replicates[1], replicates[2]), species = config$other$species, output_dir = output_dir)

nb_results = nearBaitAnalysis(my_obj, k=10)

cis_results = cisAnalysis(my_obj, k=10)

trans_results = transAnalysis(my_obj, k=20)
