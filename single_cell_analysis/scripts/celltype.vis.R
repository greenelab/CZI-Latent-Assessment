# Qiwen Hu 2018
# visulize projected features into 2-D space based on simulated and real single cell datasets
#
# Usage:
#
#  Rscript celltype.vis.R
#
# Output:
#  Visulized figures for all simulated and real single cell datasets

library(ggplot2)
library(reshape2)
library(ggbiplot)
library(cowplot)
source(file.path("util", "vis.R"))

celltye_vis <- function(input.dir, output.dir, file.list, option = "simulated"){
  #Visulize projected features from different dimension reduction approaches
  # Args:
  #  input.dir: directory contains all projected feature files
  #  output.dir: directory contains output figures
  #  file.list: list contains the names of files to be analyzed
  #  option: simulated or real datasets
  # Output:
  #  visulization figures in output directory
  #
  files <- list.files(input.dir)
  for(j in 1:nrow(file.list)){
    # get projected feature files
    feature_file <- files[grep(paste(file.list[j, ], ".exp.*", sep = ""), files)]
    # get celltype information
    cellinfo <- files[grep(paste(file.list[j, ], ".celltype.txt", sep = ""), files)]
    # list of different approaches
    option.list <- list(tybalt_plot_d3 = "tybalt_tsne_depth3", 
                        tybalt_plot_d2 = "tybalt_tsne_depth2",
                        tybalt_plot_d1 = "tybalt_tsne_depth1", 
                        tnse_plot = "rnaseq_tsne", umap_plot = "umap",
                        zifa_plot = "ZIFA", pca_plot = "pca")
    
    plot <- lapply(names(option.list), 
                   function(x) feature_vis(feature_file, cellinfo, type = option.list[[x]], data = "simulated"))
    vis.plot <- list(plot[[1]][[1]], plot[[2]][[1]], plot[[3]][[1]], plot[[4]][[1]], 
                     plot[[5]][[1]], plot[[6]][[1]])
    
    pdf_size <- length(feature_file) 
    cowplot::plot_grid(plotlist = vis.plot, ncol = 6)
    ggsave(file.path (output.dir, paste(file.list[j, ], "feature.vis.pdf", sep = ".")), 
           width = 32, height = 5)
  }
}

# visulization for simulated data
input.dir <- file.path("features", "simulated_data")
output.dir <- file.path("figures", "features_vis", "simulated_data")

# list of simulated files
sim.file.list <- read.table(file.path(input.dir, "sim.file.list.txt"))
celltye_vis(input.dir, output.dir, sim.file.list, option = "simulated")

# visulization for real data
input.dir <- file.path("features", "real_data")
output.dir <- file.path("figures", "features_vis", "real_datasets")

# list of real files
real.file.list <- read.table(file.path(input.dir, "file.list.txt"))
celltye_vis(input.dir, output.dir, real.file.list, option = "real")
