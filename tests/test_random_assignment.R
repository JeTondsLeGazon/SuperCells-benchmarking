# 10 February 2022

# Test for the random grouping function SCimplify2random which was modified to 
# incorporate the random assignement to independent groups (samples, conditions, ...)

library(testthat)
source("helpers/helper-test_random_assignment.R")


test_that("All groups are filled", {
    expect_equal(length(unique(super$membership)), length(unique(random$membership)))
})


test_that("Group sizes may vary", {
    super.table <- table(super$membership)
    random.table <- table(random$membership)
    sorted.names <- sort(names(super.table))
    super.table <- super.table[sorted.names]
    random.table <- random.table[sorted.names]
    expect_false(all(random.table == super.table))
})


test_that("Group sizes close to gamma", {
    random.mean <- mean(table(random$membership))
    # expect mean supercell size to be within 20% of gamma (arbitrary)
    expect_true(random.mean >= gamma * 4 / 5)
    expect_true(random.mean <= gamma * 6 / 5)
})
