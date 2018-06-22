# Qiwen Hu 2018
# This script is used to processing HCA melanoma augmentation datasets
# Usage: Rscript salmon.data.aug.processing.R

library(tximport)
library(dplyr)

input.dir <- file.path("data", "real_data", "HCA_data")
output.dir <- file.path("data")

dir <- input.dir
# read sample list
samples <- read.table(file.path("codes", "sample.list.txt"), sep = "\t", header = FALSE, comment.char = "")
samples <- as.factor(samples[, 1])

files <- file.path(input.dir, "mouse_melanonma_bootstrap", samples, "quant.sf")

names(files) <- samples

txi.inf.rep <- tximport(files, type = "salmon", txOut = TRUE)
count <- dplyr::bind_cols(as.data.frame(txi.inf.rep$infReps))
rownames(count) <- rownames(txi.inf.rep$counts)
write.table(count, quote = FALSE, col.names = TRUE, 
            row.names = TRUE, sep = "\t", 
            file = file.path(input.dir, "melanoma.data.aug.cout.matrix.txt"))

