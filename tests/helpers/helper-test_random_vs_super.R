# 10 February 2022

# Helper file for test_random_vs_super at gamma = 1

library(Seurat)
library(SuperCellBM)
library(SuperCell)


single_data <- readRDS('../data/hagai_mouse_lps_data/singleCellData.rds')
set.seed(0)

super <- SuperCell::SCimplify(X = single_data@assays$RNA@data, 
                              cell.annotation = single_data$label, 
                              genes.use = NULL,
                              genes.exclude = NULL, 
                              n.var.genes = 1000, 
                              k.knn = 5, 
                              gamma = 1, 
                              n.pc = 10, 
                              fast.pca = TRUE)


random <- SCimple2Random(SC = super,
                         gamma = 1, 
                         seed = 0)


ge.exact <- supercell_GE(single_data@assays$RNA@data, super$membership)
colnames(ge.exact) <- super$SC.cell.annotation.
ge.random <- supercell_GE(single_data@assays$RNA@data, random$membership)
colnames(ge.random) <- random$SC.cell.annotation.
