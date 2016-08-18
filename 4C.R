#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()
parser$add_argument("files", nargs=1, help="bedGraph files to be used")
parser$add_argument("bait", nargs=1, help="bait being used")

args <- parser$parse_args()

files <- args$files
sample <- args$samples
bait <- args$baitInfo

if(file.access(files) == -1) {
    stop(sprintf("There weren't any files provided.", files))
} else {
    files = files
}

library(R.4Cker)
library(yaml)

#wd=paste0("bedGraphs/", bait)

config = yaml.load_file("config.yaml")

conditions = unlist(strsplit(comparison, '_vs_'))

replicates = c(length(config$comparisons[[comparison]][[bait]]$control), length(config$comparisons[[comparison]][[bait]]$treatment)

samples = c(unlist(config$comparisons[[comparison]][[bait]]$control), unlist(config$comparisons[[comparison]][[bait]]$treatment))

output_dir = sprintf("../output/%s/%s/", comparison, bait)

my_obj = createR4CkerObjectFromFiles(files = c(files[1], files[2], files[3], files[4]), bait_chr = config$baits[[bait]]$bait_chr, bait_coord = config$baits[[bait]]$bait_coord, bait_name = config$baits[[bait]]$bait_name, primary_enz = config$baits[[bait]]$primary_enz, samples = c(samples[1], samples[2], samples[3], samples[4]), conditions = c(conditions[1], conditions[2]), replicates = c(replicates[1], replicates[2]), species = config$other$species, output_dir = output_dir)

nb_results = nearBaitAnalysis(my_obj, k=10)

cis_results = cisAnalysis(my_obj, k=10)

trans_results = transAnalysis(my_obj, k=20)
