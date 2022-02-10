# 8 February 2022

# Tests to check similarity at gamma = 1 of random grouping and supercells

library(testthat)

source("helpers/helper-test_random_vs_super.R")


test_that("Unique membership", {
    expect_equal(length(random$membership), length(unique(random$membership)))
    expect_equal(length(super$membership), length(unique(super$membership)))
})


test_that("Same Number of memberships", {
    expect_equal(length(random$membership), length(super$membership))
})


test_that("Same rowsums", {
    expect_equal(rowSums(ge.exact), rowSums(ge.random))
})


test_that("Same Mean per group", {
    expect_equal(rowMeans(ge.exact[, colnames(ge.exact) == 'ctrl']), 
                 rowMeans(ge.random[, colnames(ge.random) == 'ctrl']))
    expect_equal(rowMeans(ge.exact[, colnames(ge.exact) == 'treat']), 
                 rowMeans(ge.random[, colnames(ge.random) == 'treat']))
})
