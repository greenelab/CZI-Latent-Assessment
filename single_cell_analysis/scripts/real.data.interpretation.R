# Qiwen Hu 2018
# Visulization and interpretation of cluster for real single cell datasets

library(mlearning)
source(file.path("codes", "util", "vis.R"))

input.dir <- file.path("features", "real_data")
output.dir <- file.path("figures", "real_data_interpretation")

files <- list.files(input.dir)
data.list <- c("yan", "goolam", "pollen", "patel", "nestorowa", "Petropoulos", "baron-human", "melanoma")
perf <- list()

for(j in 1:length(data.list)){
  # get projected feature files
  feature.files <- files[grep(paste(data.list[j], ".*exp.*", sep = ""), files)]
  
  # get celltype information
  cellinfofiles <- files[grep(paste(data.list[j], ".*celltype.txt", sep = ""), files)]
  
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
  
  filelist <- c(tybalt_d1_files, tybalt_d2_files, tybalt_d3_files, tsne_files,
                umap_files, zifa_files, pca_files)
  # read files
  datalist <- lapply(file.path(input.dir, filelist), function(x) read.table(x, header=T))
  cellinfo <- read.table(file.path(input.dir, cellinfofiles), sep = "\t", header = TRUE)
  
  # get knn-based eval and confusion matrix
  knn.eval <- lapply(datalist, knn_eval, celltype = cellinfo$celltype)
  
  # get knn-based performance
  perf[[j]] <- as.data.frame(do.call("rbind", lapply(knn.eval, "[[", 1)))
  dataset <- rep(data.list[j], length(filelist))
  approach <- c("tybalt_depth3", "tybalt_depth2", "tybalt_depth1",
                "rnaseq_tsne", "umap", "ZIFA", "pca")
  perf[[j]]$dataset <- dataset
  perf[[j]]$approach <- approach
  
  confusion.list <- lapply(knn.eval, "[[", 2)
  option.list <- c("tybalt_depth1", "tybalt_depth2", "tybalt_depth3", "rnaseq_tsne", 
                   "umap", "ZIFA", "pca")
  
  confusion.plot <- list()
  
  for(i in 1:length(option.list)){
   confusion.plot[[i]] <- confusion_matrix_vis(knn.eval[[i]]$confusion, option = option.list[i], data = data.list[j])
  }
  
  cowplot::plot_grid(plotlist = confusion.plot, ncol = 4)
  ggsave(file.path(output.dir, paste(data.list[j], "confusion.plot.pdf", sep = ".")), 
         width = 20, height = 7)
}

perf <- dplyr::bind_rows(perf)
write.table(perf, file = file.path(output.dir, "real.data.knn.perf.txt"),
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


