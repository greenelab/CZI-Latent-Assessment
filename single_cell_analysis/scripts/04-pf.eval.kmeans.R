# Qiwen Hu 2018
# This script is used to evaluate the performance of dimension reduction algorithms
# based on different types of metrics
# Usage: Rscript pf.eval.R
#
# Output:
#  performance evaluation tables and heatmap figs 

library(NMI)
library(mclust)
library(dplyr)
library(coRanking)
source(file.path("util","util.R"))
source(file.path("util", "vis.R"))

# k-means based
kmean_pf_metric <- function(input.dir, option = "real"){
  # This function computes the k-means based performance metrics for different approaches
  # (tybalt, tsne, umap, ZIFA and PCA) for all simulated/real single cell datasets
  # Args:
  #   input.dir - input directory contains feature files and celltype informations for each dataset
  #   option - simulated/real
  # Returns:
  #   list contains NMI, ARI and average silhouette score
  
  feature.files <- list.files(input.dir)
  
  # list files contain celltype information
  cellinfo.files <- feature.files[grep(paste(".*.celltype.txt", sep = ""), feature.files)]
  
  # tybalt depth 1
  tybalt_d1_files <- feature.files[grep("depth1_.*tybalt_tsne", feature.files)]
  # tybalt depth 2
  tybalt_d2_files <- feature.files[grep("depth2_.*tybalt_tsne", feature.files)]
  # tabalt depth 3
  tybalt_d3_files <- feature.files[grep("depth3_.*tybalt_tsne", feature.files)]
  tsne_files <- feature.files[grep("depth1_rnaseq_tsne", feature.files)]
  umap_files <- feature.files[grep("umap", feature.files)]
  zifa_files <- feature.files[grep("ZIFA", feature.files)]
  pca_files <- feature.files[grep("pca", feature.files)]
  
  nmi_score <- list()
  ari_score <- list()
  silhouette_score <- list()
  
  for(i in 1:length(tybalt_d1_files)){
    filelist <- c(tybalt_d1_files[i], tybalt_d2_files[i], tybalt_d3_files[i], tsne_files[i],
                  umap_files[i], zifa_files[i], pca_files[i])
    # read files
    datalist <- lapply(file.path(input.dir, filelist), 
                       function(x) read.table(x, header=T,  comment.char = "")) 
    cellinfo <- read.table(file.path(input.dir, cellinfo.files[i]), 
                           sep = "\t", header = TRUE, comment.char = "")
    
    # extract dataset information
    dataset.info <- unlist(strsplit(tybalt_d1_files[i], split = "[.]"))
    if(option == "real"){
      dataset <- dataset.info[1]
    } else{
      dataset <- paste(dataset.info[2], dataset.info[3], sep = ".")
    }
    
    # compute NMI and ARI
    if(option == "real"){
      cl.score <- lapply(datalist, kmeans_eval, celltype = cellinfo$celltype)
    } else{
      cl.score <- lapply(datalist, kmeans_eval, celltype = cellinfo$Group)
    }
    
    nmi_score[[i]] <- data.frame(dataset, cl.score[[1]][1], cl.score[[2]][1], cl.score[[3]][1], 
                                 cl.score[[4]][1], cl.score[[5]][1], cl.score[[6]][1], cl.score[[7]][1])
    colnames(nmi_score[[i]]) <- c("files", "tybalt depth1", "tybalt depth2", "tybalt depth3", 
                                  "tnse", "umap", "ZIFA", "pca")
    ari_score[[i]] <- data.frame(dataset, cl.score[[1]][2], cl.score[[2]][2], cl.score[[3]][2], 
                                 cl.score[[4]][2], cl.score[[5]][2], cl.score[[6]][2], cl.score[[7]][2])
    colnames(ari_score[[i]]) <- c("files", "tybalt depth1", "tybalt depth2", "tybalt depth3", 
                                  "tnse", "umap", "ZIFA", "pca")
    
    # compute average sihouette score
    if(option == "real"){
      ave.sil <- lapply(datalist, ave_sil, celltype = cellinfo$celltype)
    } else{
      ave.sil <- lapply(datalist, ave_sil, celltype = cellinfo$Group)
    }
    silhouette_score[[i]] <- data.frame(dataset, ave.sil[[1]][1], ave.sil[[2]][1], ave.sil[[3]][1], 
                                        ave.sil[[4]][1], ave.sil[[5]][1], ave.sil[[6]][1], ave.sil[[7]][1])
    colnames(silhouette_score[[i]]) <- c("files", "tybalt depth1", "tybalt depth2", "tybalt depth3", 
                                         "tnse", "umap", "ZIFA", "pca")
  }
  
  nmi_score <- dplyr::bind_rows(nmi_score)
  ari_score <- dplyr::bind_rows(ari_score)
  silhouette_score <- dplyr::bind_rows(silhouette_score)
  
  return(list(nmi = nmi_score, ari = ari_score, ave.sil = silhouette_score))
}

# k-means based performance for real single cell datasets
input.dir <- file.path("features", "real_data")
output.dir <- "tables"
km.pf.metrics <- kmean_pf_metric(input.dir, option = "real")

# write performance metrics into table
write.table(km.pf.metrics$nmi, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "real.data.all_nmi_score.txt"))
write.table(km.pf.metrics$ari, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "real.data.all_ari_score.txt"))
write.table(km.pf.metrics$ave.sil, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "real.data.silhouette_score.txt"))

# Visulize performance
clustering_perf_heatmap(km.pf.metrics$nmi, type = "NMI", option = "real")
clustering_perf_heatmap(km.pf.metrics$ari, type = "ARI", option = "real")
clustering_perf_heatmap(km.pf.metrics$ave.sil, type = "Average silhouette score", option = "real")

# k-means based performance for simulated datasets

input.dir <- file.path("features", "simulated_data")
output.dir <- "tables"
km.pf.metrics <- kmean_pf_metric(input.dir, option = "simulated")

# write performance metrics into table
write.table(km.pf.metrics$nmi, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "simulated.data.all_nmi_score.txt"))
write.table(km.pf.metrics$ari, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "simulated.data.all_ari_score.txt"))
write.table(km.pf.metrics$ave.sil, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "simulated.data.silhouette_score.txt"))
# Visulization
clustering_perf_heatmap(km.pf.metrics$nmi, type = "NMI", option = "simulated")
clustering_perf_heatmap(km.pf.metrics$ari, type = "ARI", option = "simulated")
clustering_perf_heatmap(km.pf.metrics$ave.sil, type = "Average silhouette score", option = "simulated")
