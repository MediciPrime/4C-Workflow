#!/usr/bin/env Rscript

library(argparse)

# determine the files based on the comparison and bait to ensure correct file and sample order

parser <- ArgumentParser()
parser$add_argument("--bait", nargs=1, help="bait being used")
parser$add_argument("--comparison", nargs=1, help="comparison being made")
parser$add_argument("--pval", nargs=1, help="p-value for DESeq2 differential analysis")

args <- parser$parse_args()

bait <- args$bait
comparison <- args$comparison
pval <- args$pval

library(R.4Cker)
library(yaml)

config = yaml.load_file("config.yaml")

conditions = unlist(strsplit(comparison, '_vs_'))

replicates = c(length(config$comparisons[[comparison]][[bait]]$control), length(config$comparisons[[comparison]][[bait]]$treatment))

samples = c(unlist(config$comparisons[[comparison]][[bait]]$control), unlist(config$comparisons[[comparison]][[bait]]$treatment))

files = unlist(unique(lapply(samples, function(x) sprintf("bedGraphs/%s/%s_aligned_rm_self_und.bedGraph", bait, samples))))

dir.create(sprintf("output/%s/", comparison), showWarnings = FALSE)

output_dir = sprintf("output/%s/%s/", comparison, bait)

my_obj = createR4CkerObjectFromFiles(files = files, bait_chr = config$baits[[bait]]$bait_chr, bait_coord = config$baits[[bait]]$bait_coord, bait_name = config$baits[[bait]]$bait_name, primary_enz = config$baits[[bait]]$primary_enz, samples = samples, conditions = c(conditions[1], conditions[2]), replicates = c(replicates[1], replicates[2]), species = config$reference_genome$species, output_dir = output_dir)

nb_results = nearBaitAnalysis(my_obj, k=10)

cis_results = cisAnalysis(my_obj, k=10)

#trans_results = transAnalysis(my_obj, k=20)

# differential analysis via DESeq2

differentialAnalysis(obj=my_obj,
                     norm_counts_avg=nb_results$norm_counts_avg,
                     windows=nb_results$window_counts,
                     conditions=c(conditions[1], conditions[2]),
                     region="nearbait",
                     coordinates=NULL,
                     pval=pval[1])

differentialAnalysis(obj=my_obj,
                     norm_counts_avg=cis_results$norm_counts_avg,
                     windows=cis_results$window_counts,
                     conditions=c(conditions[1], conditions[2]),
                     region="cis",
                     coordinates=NULL,
                     pval=pval[1])

######################## Trans Analysis is NOT recommended #################
## differentialAnalysis(obj=my_obj,
##                      norm_counts_avg=trans_results$norm_counts_avg,
##                      windows=trans_results$window_counts,
##                      conditions=c(conditions[1], conditions[2]),
##                      region="trans",
##                      coordinates=NULL,
##                      pval=pval[1])
