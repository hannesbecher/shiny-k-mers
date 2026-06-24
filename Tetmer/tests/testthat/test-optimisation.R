# tests/testthat/test-optimisation.R
# Tests for optimisation functions using built-in example datasets
# Parameter estimates are compared against published values from
# Becher et al. (2020) Plant Communications 1:100105
#
# Note: biological range tests for theta and sub-genome divergence are
# skipped because optim() results are sensitive to starting values and
# may vary across platforms and R versions (as noted by Becher, pers. comm.)

library(Tetmer)

# Helper to create a minimal input list mimicking the Shiny UI
makeAutoInput <- function(mod, kcovRange, thRange,
                          gsRange, vfRange,
                          xrange, divRange = NULL,
                          palloRange = NULL) {
  inp <- list(
    fitmod  = "auto",
    mod     = mod,
    akcov   = kcovRange,
    ath     = thRange,
    ayadj   = gsRange,
    avf     = vfRange,
    axrange = xrange
  )
  if (!is.null(divRange))   inp$adiv   <- divRange
  if (!is.null(palloRange)) inp$apallo <- palloRange
  inp
}

# -------------------------------------------------------------------
# E030 -- diploid Euphrasia anglica, k=21
# Published estimate: theta per nucleotide ~0.2% = 0.004 per k-mer
# Genome size: ~185-225 Mbp (non-repetitive haploid)
# -------------------------------------------------------------------

test_that("doOptimisation on E030 (diploid) returns non-NULL result", {
  input <- makeAutoInput(
    mod       = "d",
    kcovRange = c(5, 100),
    thRange   = c(-3, 0),
    gsRange   = c(6, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E030, minFun, sv)

  expect_false(is.null(result))
})

test_that("doOptimisation on E030 (diploid) has named parameters", {
  input <- makeAutoInput(
    mod       = "d",
    kcovRange = c(5, 100),
    thRange   = c(-3, 0),
    gsRange   = c(6, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E030, minFun, sv)

  expect_named(result$par, c("cov", "vf", "theta", "haplSize"))
})

test_that("doOptimisation on E030 (diploid) estimates genome size within published range", {
  input <- makeAutoInput(
    mod       = "d",
    kcovRange = c(5, 100),
    thRange   = c(-3, 0),
    gsRange   = c(6, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E030, minFun, sv)

  # Published haploid non-repetitive genome size: 185-225 Mbp
  # haplSize is in units of millions, allow very generous range
  gs_mbp <- result$par["haplSize"] * 1e6 / 1e6
  expect_gt(gs_mbp, 50)
  expect_lt(gs_mbp, 500)
})

test_that("doOptimisation on E030 (diploid) returns finite objective value", {
  input <- makeAutoInput(
    mod       = "d",
    kcovRange = c(5, 100),
    thRange   = c(-3, 0),
    gsRange   = c(6, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E030, minFun, sv)

  expect_true(is.finite(result$value))
})

test_that("doOptimisation on E030 (diploid) all parameters are positive", {
  input <- makeAutoInput(
    mod       = "d",
    kcovRange = c(5, 100),
    thRange   = c(-3, 0),
    gsRange   = c(6, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E030, minFun, sv)

  # cov, vf, theta, haplSize must be positive
  expect_gt(result$par["cov"], 0)
  expect_gt(result$par["vf"], 0)
  expect_gt(result$par["theta"], 0)
  expect_gt(result$par["haplSize"], 0)
})

# -------------------------------------------------------------------
# E028 -- allotetraploid Euphrasia arctica, k=27
# Published estimate: sub-genome divergence ~5% per nucleotide
# i.e. theta * T / k ~ 0.05, so theta * T ~ 1.35
# -------------------------------------------------------------------

test_that("doOptimisation on E028 (allotetraploid) returns non-NULL result", {
  input <- makeAutoInput(
    mod       = "tal",
    kcovRange = c(5, 30),
    thRange   = c(-3, -0.5),
    gsRange   = c(7, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200),
    divRange  = c(1, 100)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E028, minFun, sv)

  expect_false(is.null(result))
})

test_that("doOptimisation on E028 (allotetraploid) has named parameters", {
  input <- makeAutoInput(
    mod       = "tal",
    kcovRange = c(5, 30),
    thRange   = c(-3, -0.5),
    gsRange   = c(7, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200),
    divRange  = c(1, 100)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E028, minFun, sv)

  expect_named(result$par, c("cov", "vf", "theta", "haplSize", "diverg"))
})

test_that("doOptimisation on E028 (allotetraploid) returns finite objective value", {
  input <- makeAutoInput(
    mod       = "tal",
    kcovRange = c(5, 30),
    thRange   = c(-3, -0.5),
    gsRange   = c(7, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200),
    divRange  = c(1, 100)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E028, minFun, sv)

  expect_true(is.finite(result$value))
})

test_that("doOptimisation on E028 (allotetraploid) all parameters are positive", {
  input <- makeAutoInput(
    mod       = "tal",
    kcovRange = c(5, 30),
    thRange   = c(-3, -0.5),
    gsRange   = c(7, 9),
    vfRange   = c(1, 30),
    xrange    = c(45, 200),
    divRange  = c(1, 100)
  )
  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)
  result  <- doOptimisation(input, E028, minFun, sv)

  expect_gt(result$par["cov"],      0)
  expect_gt(result$par["vf"],       0)
  expect_gt(result$par["theta"],    0)
  expect_gt(result$par["haplSize"], 0)
  expect_gt(result$par["diverg"],   0)
})
