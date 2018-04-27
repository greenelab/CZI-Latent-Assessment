# Qiwen Hu 2018
# This script is used to generate simulated single cell datasets
#   based on splatter simulation tool (https://github.com/Oshlack/splatter)
# USAGE: Rscript splatter.simulation.R
#
# Output:
# Simulated count matrix for different simulation parameters
# Normalized gene expression matrix for different simulation parameters
# Cell type information 

library(splatter)
library(scater)
library(SingleCellExperiment)

seed = 628540
set.seed(seed)

simulation <- function(ncells = 600, nGenes = 10000, celltypes = 5,
                       batchsize = 1, out.prob = 0.05, seed = 628540){
  # Generate simulated single cell data from splatter
  # Args:
  #  ncells: # of cells in each sample, default = 600
  #  nGenes: # of genes/features in each sample, default = 10000
  #  celltypes: # of celltypes, default = 5
  #  bathsize: # of bathsize, default = 1
  #  out.prob:  probability of a gene that become an outlier, default = 0.05
  # Returns:
  #  matrix of gene counts/expression and celltype information
  
  # calculate the probabilities that cells come from particular celltypes
  # currently even distributed
  celltype_prob = 1/celltypes
  celltype_prob_distr = rep(celltype_prob, celltypes)
  
  # set batch parameter 
  batch <- ncells / batchsize
  batchcells <- rep(batch, batchsize)
  
  # get initial parameters
  params <- splatter::newSplatParams()
  params <- splatter::setParams(params, update = list(batchCells = ncells, 
                                                      nGenes = nGenes, 
                                                      out.prob = out.prob,
                                                      seed = seed))
  
  # generate simulation data
  sim.groups <- splatter::splatSimulate(params, 
                                        batchCells = batchcells, 
                                        group.prob = celltype_prob_distr, 
                                        method = "groups", 
                                        verbose = FALSE)
  # gene length
  sim.groups <- splatter::addGeneLengths(sim.groups)
  
  # TMP normalization
  tpm(sim.groups) <- scater::calculateTPM(sim.groups, rowData(sim.groups)$Length)
  
  # celltype information
  cellinfo <- SingleCellExperiment::colData(sim.groups)
  
  # get normalized expression matrix
  exp.matrix <- SingleCellExperiment::tpm(sim.groups)
  
  # gene count matrix - log normalized
  count.matrix <- SingleCellExperiment::logcounts(scater::normalize(sim.groups))
  
  return(list(count.matrix, exp.matrix, cellinfo))
}

singlecell_simulation <- function(option = "ncells", paramvector, output.dir){
  # Generate batch simulated data for different parameters
  # and output the results into data file
  # Args:
  #
  # option: 
  #  -ncell - number of cells in a sample
  #  -nGenes - number of genes in a sample
  #  -celltypes - number of celltypes in a sample
  #  -outlier - prob of a gene that become an outlier
  # paramvetor:
  #  vector of parameters
  # output.dir:
  #  directory of output file
  #
  # Returns: null
  
  for(i in 1:length(paramvector)){
    if(option == "ncells"){
      sim_result <- simulation(ncells = paramvector[i])
      
    }
    else if(option == "nGenes"){
      sim_result <- simulation(nGenes = paramvector[i])
    }
    else if(option == "celltypes"){
      sim_result <- simulation(celltypes = paramvector[i])
    }
    else if(option == "batchsize"){
      sim_result <- simulation(batchsize = paramvector[i])
    }
    else if(option == "outlier"){
      sim_result <- simulation(out.prob = paramvector[i])
    }
  }
  write.table(sim_result[[1]], 
              file = file.path(output.dir, paste("sim", 
                                                 option, paramvector[i], 
                                                 "count.matrix.txt", 
                                                 sep = ".")), 
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
  write.table(sim_result[[2]], 
              file.path(output.dir, paste("sim", 
                                          option, paramvector[i],
                                          "exp.matrix.txt", 
                                          sep = ".")),
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  write.table(sim_result[[3]], 
              file.path(output.dir, paste("sim", 
                                          option, paramvector[i], 
                                          "celltype.txt", sep = ".")),
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}

# Generate simulated data
output.dir <- file.path("data", "simulated")

# simulate different ncells
ncells_sim <- c(seq(500, 5000, 300), seq(10000, 50000, 10000))
singlecell_simulation(option = "ncells", ncells_sim, output.dir)

# nGenes
ngenes_sim <- seq(20000, 60000, 5000)
singlecell_simulation(option = "nGenes", ngenes_sim, output.dir)

# celltypes
celltypes_sim <- seq(5, 15, 1)
singlecell_simulation(option = "celltypes", celltypes_sim, output.dir)

# outliers
outlier_sim <- seq(0.1, 0.5, 0.1)
singlecell_simulation(option = "outlier", outlier_sim, output.dir)
