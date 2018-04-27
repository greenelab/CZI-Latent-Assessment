# Qiwen Hu 2018
# This script is used to calculate the coranking metrics for real single cell datasets
# Usage:
#
#  Rscript coranking.eval.R
#
# Output:
#  table contains q.local/q.global for all real datasets

library(NMI)
library(mclust)
library(dplyr)
library(gplots)
library(coRanking)
library(readr)
library(tibble)
source(file.path("util", "util.R"))


input.dir <- file.path("features", "real_data")
file.dir <- file.path("data", "real_data")
output.dir <- file.path("tables")

# list all files contain projected features
feature.files <- list.files(input.dir)

# list all original count files
count.files <- list.files(file.dir)
count.files <- count.files[grep(paste(".*.zeroone_rnaseq.tsv", sep = ""), count.files)] 

tybalt_d1_files <- feature.files[grep(".*depth1_tybalt_tsne.*", feature.files)]
tybalt_d2_files <- feature.files[grep("depth2_tybalt_tsne", feature.files)]
tybalt_d3_files <- feature.files[grep("depth3_tybalt_tsne", feature.files)]
tsne_files <- feature.files[grep("depth1_rnaseq_tsne", feature.files)]
umap_files <- feature.files[grep("umap", feature.files)]
zifa_files <- feature.files[grep("ZIFA", feature.files)]
pca_files <- feature.files[grep("pca", feature.files)]

q.local <- list()
q.global <- list()
for(i in 1:length(tybalt_d1_files)){
  tybalt_d1_feature <- read.table(file.path(input.dir, tybalt_d1_files[i]), sep = "\t", header = TRUE)
  tybalt_d2_feature <- read.table(file.path(input.dir, tybalt_d2_files[i]), sep = "\t", header = TRUE)
  tybalt_d3_feature <- read.table(file.path(input.dir, tybalt_d3_files[i]), sep = "\t", header = TRUE)
  tsne_feature <- read.table(file.path(input.dir, tsne_files[i]), sep = "\t", header = TRUE)
  umap_feature <- read.table(file.path(input.dir, umap_files[i]), sep = "\t", header = TRUE)
  zifa_feature <- read.table(file.path(input.dir, zifa_files[i]), sep = "\t", header = TRUE)
  pca_feature <- read.table(file.path(input.dir, pca_files[i]), sep = "\t", header = TRUE)
  pca_feature <- tibble::add_column(pca_feature, sample = rownames(pca_feature), .before = "pca.1")
  
  file.list <- list(tybalt_d1 = tybalt_d1_feature, tybalt_d2 = tybalt_d2_feature, tybalt_d3 = tybalt_d3_feature, 
                    tsne = tsne_feature, umap = umap_feature, zifa = zifa_feature, pca = pca_feature)
  
  count.matrix <- readr::read_tsv(file.path(file.dir, count.files[i]))
  
  #remove gene names
  count.matrix <- subset(count.matrix, select = -c(X1))
  
  coranking.eval <- lapply(names(file.list), function(x) coranking_metric(count.matrix, file.list[[x]]))
  
  param <- unlist(strsplit(tybalt_d1_files[i], split = "[.]"))[1]
  
  q.local[[i]] <- data.frame(param, t(unlist(lapply(coranking.eval, `[[`, 1))))
  colnames(q.local[[i]]) <- c("files", "tybalt depth1", "tybalt depth2", "tybalt depth3",
                              "tnse", "umap", "ZIFA", "pca")
  q.global[[i]] <- data.frame(param, t(unlist(lapply(coranking.eval, `[[`, 2))))
  colnames(q.global[[i]]) <- c("files", "tybalt depth1", "tybalt depth2", "tybalt depth3",
                               "tnse", "umap", "ZIFA", "pca")
  

}

q.local <- dplyr::bind_rows(q.local)
q.global <- dplyr::bind_rows(q.global)

write.table(q.local, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "real.data.q.local.txt"))

write.table(q.global, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "real.data.q.global.txt"))
