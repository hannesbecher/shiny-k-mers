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
                      kcov = 30, vf = 1.5,
                      theta = 0.04, gs = 200)

  expect_length(result, 100)
})

test_that("evalModel output is numeric and finite", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result <- evalModel(probs, factors,
                      xmin = 1, xmax = 100,
                      kcov = 30, vf = 1.5,
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
                      kcov = 30, vf = 1.5,
                      theta = 0.04, gs = 200)

  expect_true(all(result >= 0))
})

test_that("evalModel output scales linearly with genome size", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result1 <- evalModel(probs, factors,
                       xmin = 1, xmax = 100,
                       kcov = 30, vf = 1.5,
                       theta = 0.04, gs = 100)

  result2 <- evalModel(probs, factors,
                       xmin = 1, xmax = 100,
                       kcov = 30, vf = 1.5,
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
                          kcov = 30, vf = 1.5,
                          theta = 0.04, gs = 200)

  result_tet <- evalModel(probs_tet, factors_tet,
                          xmin = 1, xmax = 100,
                          kcov = 30, vf = 1.5,
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
                          kcov = 30, vf = 1.5,
                          theta = 0.04, gs = 200)

  result_all <- evalModel(probs_all, factors_all,
                          xmin = 1, xmax = 100,
                          kcov = 30, vf = 1.5,
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
                          kcov = 30, vf = 1.5,
                          theta = 0.01, gs = 200)

  result_high <- evalModel(probs, factors,
                           xmin = 1, xmax = 200,
                           kcov = 30, vf = 1.5,
                           theta = 0.5, gs = 200)

  expect_gt(sum(result_high), sum(result_low))
})

test_that("evalModel works with named numeric inputs (strips names correctly)", {
  input <- makeInput("d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  named_kcov <- c(cov = 30)
  named_vf <- c(vf = 1.5)

  result <- evalModel(probs, factors,
                      xmin = 1, xmax = 50,
                      kcov = named_kcov, vf = named_vf,
                      theta = 0.04, gs = 200)

  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

# -------------------------------------------------------------------
# Tests using built-in example datasets
# E030 -- diploid Euphrasia anglica, k=21
# E028 -- allotetraploid Euphrasia arctica, k=27
# -------------------------------------------------------------------

test_that("evalModel diploid output has same length a E030 spectrum range", {
  input   <- list(mod = "d")
  probs   <- getProbs(input)
  factors <- getFactors(input)
  xmin    <- 45
  xmax    <- 200

  result  <- evalModel(probs, factors,
                       xmin = xmin, xmax = xmax,
                       kcov = 15, vf = 1.5,
                       theta = 0.04, gs = 200)

  expect_length(result, xmax - xmin + 1)
})

test_that("evalModel diploid peak is near the sequencing coverage of E030", {
  input   <- list(mod = "d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  # E030 has known haploid coverage around 25-30x
  # The homozygous (2x) peak should be the large peak
  result <- evalModel(probs, factors,
                      xmin = 1, xmax = 200,
                      kcov = 27, vf = 1.5,
                      theta = 0.01, gs = 200)

  # Peak of the homozygous component should be near 2*kcov = 54
  peak_position <- which.max(result)
  expect_gt(peak_position, 20)
  expect_lt(peak_position, 100)
})

test_that("evalModel allotetraploid output has same length as E028 spectrum range", {
  input   <- list(mod = "tal")
  probs   <- getProbs(input)
  factors <- getFactors(input)
  xmin    <- 45
  xmax    <- 200

  result <- evalModel(probs, factors,
                      xmin = xmin, xmax = xmax,
                      kcov = 15, vf = 1.5,
                      theta = 0.04, gs = 200,
                      diverg = 30)

  expect_length(result, xmax -xmin + 1)
})

test_that("evalModel allotetraploid with high divergence has larger 2x peak than 1x peak", {
  # This is the key biological property of allotetraploids vs autotetraploids
  # A prominent 2x peak relative to 1x peak is the signature of allotetraploidy
  input   <- list(mod = "tal")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result  <- evalModel(probs, factors,
                       xmin = 1, xmax = 200, 
                       kcov = 15, vf = 1.5,
                       theta = 0.04, gs = 200,
                       diverg = 50)

  peak_1x <- result[15]    # value at 1x coverage
  peak_2x <- result[30]    # value at 2x coverage

  expect_gt(peak_2x, peak_1x)
})

test_that("evalModel diploid output is of the same order of magnitude as E030 data", {
  input   <- list(mod = "d")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  result  <- evalModel(probs, factors,
                       xmin = 45, xmax = 200,
                       kcov = 57, vf = 1.5,
                       theta = 0.004, gs = 200)

  observed <- E030@data$count[45:200]

  # Model and data should be in the same order of magnitude
  # i.e. their ratio should be between 0.01 and 100
  ratio <- mean(result) / mean(observed)
  expect_gt(ratio, 0.01)
  expect_lt(ratio, 100)
})

test_that("evalModel allotetraploid fitted parameters reproduce E028 spectrum shape", {
  input   <- list(mod = "tal")
  probs   <- getProbs(input)
  factors <- getFactors(input)

  # Use published-approximate parameters for E028
  # sub-genome divergence ~5% per nucleotide, k=27
  # theta ~ 0.04 per k-mer, T ~38, coverage ~15x, gs ~200 Mbp
  result <- evalModel(probs, factors,
                      xmin = 45, xmax = 200,
                      kcov = 15, vf = 1.5,
                      theta = 0.04, gs = 200,
                      diverg = 38)

  # Compare shape to actual E028 data in the same range
  observed <- E028@data$count[45:200]

  # Pearson correlation between model and data should be positive
  expect_gt(cor(result, observed), 0)
})

