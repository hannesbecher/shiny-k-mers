# tests/testthat/test-evalModel.R

library(Tetmer)

# Helper: create a minimal input-like list for a given model
makeInput <- function(mod) {
  list(mod = mod)
}

test_that("evalModel returns correct length vector for diploid model", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result <- evalModel(probs, factors,
                      xmin = 1, xmax = 100,
                      kcov = 30, bias = -1.8,
                      theta = 0.04, gs = 200)

  expect_length(result, 100)
})

test_that("evalModel output is numeric and finite", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result <- evalModel(probs, factors,
                      xmin = 1, xmax = 100,
                      kcov = 30, bias = -1.8,
                      theta = 0.04, gs = 200)

  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("evalModel output is non-negative", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result <- evalModel(probs, factors,
                      xmin = 1, xmax = 100,
                      kcov = 30, bias = -1.8,
                      theta = 0.04, gs = 200)

  expect_true(all(result >= 0))
})

test_that("evalModel output scales linearly with genome size", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result1 <- evalModel(probs, factors,
                       xmin = 1, xmax = 100,
                       kcov = 30, bias = -1.8,
                       theta = 0.04, gs = 100)

  result2 <- evalModel(probs, factors,
                       xmin = 1, xmax = 100,
                       kcov = 30, bias = -1.8,
                       theta = 0.04, gs = 200)

  expect_equal(result2, result1 * 2)
})

test_that("evalModel diploid and autotetraploid models give different results", {
  probs_dip <- getProbs(makeInput("d"))
  factors_dip <- getFactors(makeInput("d"))
  probs_tet <- getProbs(makeInput("tau"))
  factors_tet <- getFactors(makeInput("tau"))

  result_dip <- evalModel(probs_dip, factors_dip,
                          xmin = 1, xmax = 100,
                          kcov = 30, bias = -1.8,
                          theta = 0.04, gs = 200)

  result_tet <- evalModel(probs_tet, factors_tet,
                          xmin = 1, xmax = 100,
                          kcov = 30, bias = -1.8,
                          theta = 0.04, gs = 200)

  expect_false(identical(result_dip, result_tet))
})

test_that("evalModel allotetraploid model requires diverg and gives different result to autotetraploid", {
  probs_aut <- getProbs(makeInput("tau"))
  factors_aut <- getFactors(makeInput("tau"))
  probs_all <- getProbs(makeInput("tal"))
  factors_all <- getFactors(makeInput("tal"))

  result_aut <- evalModel(probs_aut, factors_aut,
                          xmin = 1, xmax = 100,
                          kcov = 30, bias = -1.8,
                          theta = 0.04, gs = 200)

  result_all <- evalModel(probs_all, factors_all,
                          xmin = 1, xmax = 100,
                          kcov = 30, bias = -1.8,
                          theta = 0.04, gs = 200,
                          diverg = 30)

  expect_false(identical(result_aut, result_all))
})

test_that("evalModel sum increases with higher theta", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result_low <- evalModel(probs, factors,
                          xmin = 1, xmax = 200,
                          kcov = 30, bias = -1.8,
                          theta = 0.01, gs = 200)

  result_high <- evalModel(probs, factors,
                           xmin = 1, xmax = 200,
                           kcov = 30, bias = -1.8,
                           theta = 0.5, gs = 200)

  expect_gt(sum(result_high), sum(result_low))
})

test_that("evalModel works with named numeric inputs (strips names correctly)", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  named_kcov <- c(cov = 30)
  named_bias <- c(bias = -1.8)

  result <- evalModel(probs, factors,
                      xmin = 1, xmax = 50,
                      kcov = named_kcov, bias = named_bias,
                      theta = 0.04, gs = 200)

  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})
