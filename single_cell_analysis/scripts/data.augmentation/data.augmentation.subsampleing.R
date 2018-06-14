# Qiwen Hu 2018
# This script is unsed ot subsampling melanoma datasets (augmentated)
# Usage: Rscript data.augmentation.subsampling.R

library(readr)
library(dplyr)
set.seed("12345")

cell_sampling <- function(ncells, celltype, seed = 12345){
  # This function is used to subsample number of cells from single-cell dataset
  # Args:
  #  ncells: number of cells to sampling for each celltype
  #  celltype: a file contains sample id and celltype information
  # Returns
  #  a data frame of subsampling of samples for each celltype
  set.seed(seed)
  celltype_id <- unique(celltype$celltype)
  sampling <- list()
  for(i in 1:length(celltype_id)){
    sub.cell.type <- celltype[celltype$celltype == celltype_id[i], ]
    sampling[[i]] <- data.frame(id = sample(sub.cell.type[, 1], ncells), 
                            celltype = rep(celltype_id[i], ncells))
  }
  sampling <- dplyr::bind_rows(sampling)
  return(sampling)
} 

subsampling <- function(c.matrix, celltype, subsample.num, aug.num, seed = seed){
  # This function is used to do subsampling from data augmentation
  # Args:
  #  c.matrix: originan count matrix for single-cell data
  #  celltype: a dataframe contains samples and celltype information
  #  subsample.num: number of cells for subsampling
  #  aug.num: subsampling augmentation
  # Returns:
  #  a list contains celltype information and its correspondent count matrix
  celltype_id <- unique(celltype$celltype)
  ncells <- round(subsample.num / length(celltype_id), 0)
  cell_sample <- cell_sampling(ncells, celltype)
  
  total.aug <- as.vector(1:10)
  subsample.result <- list()
  for(i in 1:nrow(cell_sample)){
    subsample.result[[i]] <- data.frame(id = paste(cell_sample[i, ]$id, sample(total.aug, aug.num), sep = "."),
                                        celltype = rep(cell_sample[i, ]$celltype, aug.num))
  }
  subsample.result <- dplyr::bind_rows(subsample.result)
  count.matrix <- c.matrix[, which(colnames(c.matrix) %in% subsample.result$id)]
  rownames(count.matrix) <- c.matrix$id
  return(list(celltype.info = subsample.result, count = count.matrix))
}

input.dir <- "data/real_data/HCA_data"
output.dir <- "data/data_aug"

subsample.num <- c(100, 500, 1000, 2000)
aug.num <- c(1, 2, 5, 10)

# read melanoma file with augmentation
file <- read_tsv(file.path(input.dir, "melanoma.data.aug.txtscaled_zeroone_rnaseq.tsv"))
#change colnames of file
colnames(file)[1] <- "id"
colnames(file)[-1] <- gsub("X", "", colnames(file)[-1])
colnames(file)[-1] <- gsub("_2", "", colnames(file)[-1])
colnames(file)[-1] <- sub("[.]", "#", colnames(file)[-1])
sample_id <- data.frame(aug.id = colnames(file)[-1])

celltype <- read.table(file.path(input.dir, "melanoma.celltype.txt"), sep = "\t", 
                       header = TRUE, comment.char = "")

for(i in subsample.num){
  for(j in aug.num){
    sample_data <- subsampling(file, celltype, i, j)
    
    # write celltype information to file
    write.table(sample_data$celltype.info, 
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t", 
                file = file.path(output.dir, paste("melanoma", "ncell", i, "aug", j, "celltype.txt", sep = ".")))
    
    # write count matrix to file
    write.table(sample_data$count, 
                quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t", 
                file = file.path(output.dir, paste("melanoma", "ncell", i, "aug", j, "count.txt", sep = ".")))
  }
}



