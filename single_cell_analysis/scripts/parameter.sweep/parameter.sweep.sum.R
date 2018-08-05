# Qiwen Hu 2018
# This script is used to calculate k-means performance based on parameter sweep result
# Usage Rscript parameter.sweep.sum.R

source("pf.eval.kmeans.R")

# read parameter sweep files
input.dir <- file.path("features", "param.sweep")
output.dir <- file.path("tables")

feature.files <- list.files(input.dir)

cellinfo.files <- feature.files[grep(paste("sim.ncell.30000.*.cell.info.txt", sep = ""), feature.files)]
count.files <- feature.files[grep(paste("sim.ncell.30000.*.features.tsv", sep = ""), feature.files)] 

nmi_score <- list()
ari_score <- list()
for(i in 1:length(count.files)){
  tybalt_feature <- read.table(file.path(input.dir, count.files[i]), sep = "\t",
                               header = TRUE, comment.char = "")
  
  cellinfo <- read.table(file.path(input.dir, cellinfo.files), sep = "\t", 
                         header = TRUE, comment.char = "")
  
  tybalt_feature <- merge(tybalt_feature, cellinfo, by = c(1))
  
  tybalt_score <- round(clustering_eval(tybalt_feature[, 2:3], tybalt_feature$Group), 2)
  
  param <- count.files[i]
  
  nmi_score[[i]] <- data.frame(param, tybalt_score[1])
  colnames(nmi_score[[i]]) <- c("param", "nmi")
  
  ari_score[[i]] <- data.frame(param, tybalt_score[2])
  colnames(ari_score[[i]]) <- c("param", "ari")
  
}

nmi_score <- dplyr::bind_rows(nmi_score)
ari_score <- dplyr::bind_rows(ari_score)

write.table(nmi_score, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "sim.ncell.30000.nmi.txt"))

write.table(ari_score, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep = "\t", 
            file = file.path(output.dir, "sim.ncell.30000.ari.txt"))

# visulization
input.dir <- file.path("features", "best.param.com")
output.dir <- file.path("tables")

feature.files <- list.files(input.dir)

cellinfo.files <- feature.files[grep(paste("sim.ncell.30000.*.cell.info.txt", sep = ""), feature.files)]
count.files <- feature.files[grep(paste("sim.ncell.30000.*.features.tsv", sep = ""), feature.files)] 

feature_plot <- list()
for(i in 1:length(count.files)){
  file <- read.table(file.path(input.dir, count.files[i]), header = TRUE,
                       sep = "\t", comment.char = "")
  celltype <- read.table(file.path(input.dir, cellinfo.files), header = TRUE,
                           sep = "\t", comment.char = "")
  file <- cbind(file, celltype[, 3])
  colnames(file) <- c("id", "X1", "X2", "Group")
  # remove outlier for visulization
  file <- file[file$X1 > -50 & file$X2 > -50, ]
  file <- file[file$X1 < 100 & file$X2 < 100, ]
  
  gtitle <- ""
  if(i == 1){
    gtitle <- "depth2, before tuning, ncell = 30,000"
  } else if(i == 2){
    gtitle <- "depth3, before tuning, ncell = 30,000"
  } else if(i == 3){
    gtitle <- "depth2, after tuning, ncell = 30,000"
  } else{
    gtitle <- "depth3, after tuning, ncell = 30,000"
  }
  
  feature_plot[[i]] <- ggplot2::ggplot(file, aes(X1, X2)) + 
      geom_point(aes(colour = factor(Group)), size = 1) + 
      ggtitle(gtitle) + xlab("") + ylab("") + 
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = element_text(hjust = 0.5))
}

pdf_size <- length(feature_plot)
cowplot::plot_grid(plotlist = feature_plot, ncol = 4)
ggsave(file.path (output.dir, "ncell.30000.feature.vis.pdf"), 
       width = 16, height = 5)

