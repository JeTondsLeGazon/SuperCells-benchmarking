# 8 February 2022

# Tests for different functions of SuperCells-benchmarking and implemented 
# modifications of SupercellsBM
library(Seurat)
library(SuperCellBM)
library(SuperCell)
library(testthat)

source('src/utility.R')
source('src/supercells.R')
source('src/analysis.R')
source('src/compute_de.R')

single_data <- readRDS('data/hagai_mouse_lps_data/singleCellClusteredNormalized.rds')
set.seed(0)


SC.list <- compute_supercells(
    sc = single_data,
    ToComputeSC = T,
    data.folder = 'data',
    filename = 'tmp',
    gamma.seq = 1,
    n.var.genes = 1000,
    k.knn = 5,
    n.pc = 10,
    approx.N = 1000,
    fast.pca = TRUE,
    genes.use = NULL, 
    genes.exclude = NULL,
    seed.seq = c(0),
    split.by = 'label'
)
super.random <- SC.list$Random$`1`$`0`
super.exact <- SC.list$Exact$`1`$`0`

test_that("Unique membership", {
    expect_equal(length(super.random$membership), length(unique(super.random$membership)))
    expect_equal(length(super.exact$membership), length(unique(super.exact$membership)))
})

test_that("Same Number of memberships", {
    expect_equal(length(super.random$membership), length(super.exact$membership))
})

# test_that("Correct membership group", {
#     id.treat <- which(super.exact$SC.cell.annotation. == 'treat')
#     id.ctrl <- which(super.exact$SC.cell.annotation. == 'ctrl')
#     
#     id.treat <- which(super.random$SC.cell.annotation. == 'treat')
#     id.ctrl <- which(super.random$SC.cell.annotation. == 'ctrl')
#     
#     expect_equal(sum(super.random$), length(super.exact$membership))
# })

ge.exact <- supercell_GE(single_data@assays$RNA@data, super.exact$membership)
colnames(ge.exact) <- super.exact$SC.cell.annotation.
ge.random <- supercell_GE(single_data@assays$RNA@data, super.random$membership)
colnames(ge.random) <- super.random$SC.cell.annotation.

test_that("Same rowsums", {
    expect_equal(rowSums(ge.exact), rowSums(ge.random))
})

test_that("Same Mean per group", {
    expect_equal(rowMeans(ge.exact[, colnames(ge.exact) == 'ctrl']), 
                 rowMeans(ge.random[, colnames(ge.random) == 'ctrl']))
    expect_equal(rowMeans(ge.exact[, colnames(ge.exact) == 'treat']), 
                 rowMeans(ge.random[, colnames(ge.random) == 'treat']))
})
