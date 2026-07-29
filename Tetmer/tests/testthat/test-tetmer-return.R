library(Tetmer)

test_that("resolveTetmerAppResult returns original spectrum when app result is NULL", {
  original <- E028
  resolved <- Tetmer:::resolveTetmerAppResult(NULL, original)

  expect_s4_class(resolved, "spectrum")
  expect_identical(resolved, original)
})

test_that("resolveTetmerAppResult returns app result when non-NULL", {
  original <- E028
  fit <- methods::new("tetmerFit",
    fitType = "auto",
    model = "d",
    par = c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200),
    ranges = list(),
    xrange = c(10, 80),
    convergence = 0,
    value = 1,
    k = original@k,
    version = as.character(utils::packageVersion("Tetmer"))
  )
  updated <- addFit(original, fit)

  resolved <- Tetmer:::resolveTetmerAppResult(updated, original)

  expect_identical(resolved, updated)
  expect_equal(length(resolved@fits), length(original@fits) + 1)
})

test_that("makeFitRecord creates a manual fit from current input values", {
  input <- list(
    fitmod = "man",
    mod = "d",
    tkcov = 30,
    tvf = 1.5,
    tth = 0.04,
    tyadj = 200,
    txmax = 80
  )

  fit <- Tetmer:::makeFitRecord(input, optimised = NULL, k = 21)

  expect_s4_class(fit, "tetmerFit")
  expect_identical(fit@fitType, "manual")
  expect_identical(fit@model, "d")
  expect_equal(fit@par, c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200))
  expect_equal(fit@ranges, list())
  expect_equal(fit@xrange, c(1, 80))
  expect_true(is.na(fit@convergence))
  expect_true(is.na(fit@value))
})
