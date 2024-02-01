#test functions
library(testthat)

source("~/Projects/typewriter_analysis/src/cell_culture/understand_growth/functions.R")

# unit tests for prob_zero_descendants ----------------------------------
test_that("prob_zero_descendants returns correct values", {
  expect_equal(prob_zero_descendants(2, 1, 1), 0.3873002, tolerance = 1e-7)
})

test_that("prob_zero_descendants handles division by zero", {
  expect_error(prob_zero_descendants(1, 1, 1))
})

# unit tests for prob_one_descendant ----------------------------------
test_that("prob_one_descendant returns correct values", {
  expect_equal(prob_one_descendant(2, 1, 1), 0.1381023, tolerance = 1e-7)
})

test_that("prob_one_descendant handles division by zero", {
  expect_error(prob_one_descendant(1, 1, 1))
})

# unit tests for prob_n_descendants ----------------------------------
test_that("prob_n_descendants returns correct values", {
  lambda <- 2
  mu <- 1
  n <- 2
  t <- 1
  p_1 <- prob_one_descendant(lambda, mu, t)
  p_0 <- prob_zero_descendants(lambda, mu, t)

  expected_value <- (lambda / mu * p_0)^(n - 1) * p_1
  expect_equal(prob_n_descendants(lambda, mu, t, n), expected_value)
})


test_that("p_n_t handles division by zero", {
  expect_error(prob_zero_descendants(1, 1, 2, 1))
  expect_error(prob_zero_descendants(1, 0, 2, 1))
})

