# Qiwen Hu 2018
# This script is unsed ot subsampling melanoma datasets (augmentated)
# Usage: Rscript data.augmentation.subsampling.R
# Output: The output of this script is the subsampling count matrix and
#          its correspondent celltype information under different sampling
#          parameters

library(readr)
library(dplyr)
set.seed("12345")

cell_sampling <- function(ncells, celltype, seed = 12345){
  # This function is used to sample number of cells from single-cell dataset
  # Args:
  #  ncells: number of cells to sample for each celltype
  #  celltype: a file contains sample id and celltype information
  # Returns
  #  a data frame of subsampling of samples for each celltype
  
  set.seed(seed)
  celltype_id <- unique(celltype$celltype)
  
  sampling <- celltype %>% dplyr::group_by(cell_type) %>% dplyr::sample_n(ncells)
  
  return(sampling)
} 

subsampling <- function(c.matrix, celltype, subsample.num, aug.num, seed = 12345){
  # This function is used to do subsampling from data augmentation
  # Args:
  #  c.matrix: original count matrix for single-cell data
  #  celltype: a dataframe contains samples and celltype information
  #  subsample.num: number of cells for subsampling
  #  aug.num: subsampling augmentation
  #  seed: random seed 
  # Returns:
  #  a list contains celltype information and its correspondent count matrix
  
  set.seed(seed)
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
  count.matrix <- c.matrix %>% dplyr::select(subsample.result$id)
  rownames(count.matrix) <- c.matrix$id
  return(list(celltype.info = subsample.result, count = count.matrix))
}

input.dir <- file.path("data", "real_data", "HCA_data")
output.dir <- file.path("data", "data_aug")

subsample.num <- c(100, 500, 1000, 2000)
aug.num <- c(1, 2, 5, 10)

# read melanoma file with augmentation
file <- read_tsv(file.path(input.dir, "melanoma.data.aug.txtscaled_zeroone_rnaseq.tsv"))

# change colnames of file
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
                file = file.path(output.dir, 
                                 paste("melanoma", "ncell", i, "aug", j, "celltype.txt", sep = ".")))
    
    # write count matrix to file
    write.table(sample_data$count, 
                quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t", 
                file = file.path(output.dir, 
                                 paste("melanoma", "ncell", i, "aug", j, "count.txt", sep = ".")))
  }
}
