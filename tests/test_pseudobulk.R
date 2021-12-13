source('../processing.R')

library(testthat)
library(Seurat)

test_mat <- data.frame(ctrl1 = c(1, 2, 3), ctrl2 = c(2, 3, 4), ctrl3 = c(1, 8, 3), ctrl4 = c(7, 12, 0),
                       treat1 = c(5, 1, 8), treat2 = c(10, 0, 2), treat3 = c(1, 1, 9), treat4 = c(0, 0, 4),
                       row.names = c('G1', 'G2', 'G3'))
meta <- data.frame(row.names = colnames(test_mat), 
                   label = c(rep('ctrl', 4), rep('treat', 4)),
                   sample = c(rep('mouse1_ctrl', 2), rep('mouse2_ctrl', 2), rep('mouse1_treat', 2), rep('mouse2_treat', 2)))
sc <- CreateSeuratObject(counts = test_mat, meta.data = meta)

gt <- data.frame(ctrl1 = test_mat$ctrl1 + test_mat$ctrl2,
                 treat1 = test_mat$treat1 + test_mat$treat2,
                 ctrl2 = test_mat$ctrl3 + test_mat$ctrl4,
                 treat2 = test_mat$treat3 + test_mat$treat4,
                 row.names = c('G1', 'G2', 'G3'))
got <- create_pseudobulk(sc)
    
test_that('Correct Aggregation',
          expect_equal(matrix(unlist(gt)), matrix(GetAssayData(got))))