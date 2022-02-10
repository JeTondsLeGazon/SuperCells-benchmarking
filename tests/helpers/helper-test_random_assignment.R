# 10 February 2022

# Data input for testthat functions (helper file)

library(Seurat)
library(SuperCellBM)
library(SuperCell)


single_data <- readRDS('../data/hagai_mouse_lps_data/singleCellClusteredNormalized.rds')
set.seed(0)
gamma = 10

super <- SuperCell::SCimplify(X = single_data@assays$RNA@data, 
                              cell.annotation = single_data$label, 
                              genes.use = NULL, 
                              genes.exclude = NULL, 
                              n.var.genes = 1000, 
                              k.knn = 5, 
                              gamma = gamma, 
                              n.pc = 10, 
                              fast.pca = TRUE)

random <- SCimple2Random(SC = super,
                         gamma = gamma,
                         seed = 0)