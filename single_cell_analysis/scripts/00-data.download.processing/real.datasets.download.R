# Qiwen Hu 2018
# This script is used to download and process real single cell datasets
# datasets were downloaded from https://hemberg-lab.github.io/scRNA.seq.datasets/
# Usage:
#
#   Rscript single.cell.data.processing.R
#
# Output: count matrix and celltype informaiton files for each dataset


library(SingleCellExperiment)

input.dir <- file.path("scripts", "00-data.download.processing")
output.dir <- file.path("data", "real_data")
url.link <- "https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects"

# read real data list
data.list <- read.table(file.path(input.dir, "real.data.list.txt"),
                        header = FALSE)

for(i in 1:nrow(data.list)){
  # Read rds from website
  rds <- readRDS(url(paste(url.link, data.list[i, ], sep = "/")))
  
  # Extract count matrix
  rds.count <- SingleCellExperiment::logcounts(rds)
  
  # Extract celltype information
  celltype <- data.frame(sample = rownames(SingleCellExperiment::colData(rds)), 
                         celltype = SingleCellExperiment::colData(rds)$cell_type1)
  
  write.table(rds.count, file = file.path(output.dir, paste(data.list[i], "exp.matrix.txt", sep = ".")), 
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
  write.table(celltype, file = file.path(output.dir, paste(data.list[i], "celltype.txt", sep = ".")),
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}


