# Qiwen Hu 2018
# This script is used to evaluate the performance of dimension reduction algorithms
# based on different types of metrics
# Usage: Rscript pf.eval.knn.R
#
# Output:
#  performance evaluation tables 

library(caret)
source(file.path("util","util.R"))
source(file.path("util", "vis.R"))

knn_perf_metric <- function(input.dir, option = "real"){
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
  
  perf <- list()
  
  for(i in 1:length(tybalt_d1_files)){
    filelist <- c(tybalt_d1_files[i], tybalt_d2_files[i], tybalt_d3_files[i], tsne_files[i],
                umap_files[i], zifa_files[i], pca_files[i])
    # read files
    datalist = lapply(file.path(input.dir, filelist), 
                      function(x) read.table(x, header=T, comment.char = ""))
    cellinfo <- read.table(file.path(input.dir, cellinfo.files[i]), sep = "\t", 
                           header = TRUE, comment.char = "")
    
    # extract dataset information
    dataset.info <- unlist(strsplit(tybalt_d1_files[i], split = "[.]"))
    if (option == "real"){
      dataset <- dataset.info[1]
    } else{
      dataset <- paste(dataset.info[2], dataset.info[3], sep = ".")
    }
    
    # get knn-based eval and confusion matrix
    if(option == "real"){
      knn.eval <- lapply(datalist, knn_eval, celltype = cellinfo$celltype)
    } else{
      knn.eval <- lapply(datalist, knn_eval, celltype = cellinfo$Group)
    }
    
    # get knn-based performance
    perf[[i]] <- as.data.frame(do.call("rbind", lapply(knn.eval, "[[", 1)))
    dataset <- rep(dataset, length(filelist))
    approach <- c("tybalt_depth3", "tybalt_depth2", "tybalt_depth1",
                  "rnaseq_tsne", "umap", "ZIFA", "pca")
    perf[[i]]$dataset <- dataset
    perf[[i]]$approach <- approach
    }
  perf <- dplyr::bind_rows(perf)
  return(perf)
}

# knn performance for simulated data
input.dir <- file.path("features", "simulated_data")
output.dir <- "tables"

knn.pf.metrics <- knn_perf_metric(input.dir, option = "simulated")
write.table(knn.pf.metrics, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "simulated.knn.pf.metrics.txt"))

# visulization for model performance 
knn_perf_vis(knn.pf.metrics, type = "simulated")

# knn performance for real data
input.dir <- file.path("features", "real_data")
output.dir <- "tables"

knn.pf.metrics <- knn_perf_metric(input.dir, option = "real")
write.table(knn.pf.metrics, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "real.knn.pf.metrics.txt"))

# visulization for model performance
knn_perf_vis(knn.pf.metrics, type = "real")
