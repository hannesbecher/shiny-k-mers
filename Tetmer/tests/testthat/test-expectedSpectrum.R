library(Tetmer)

makeExpectedSpectrumTestFit <- function(model = "d") {
  par <- c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200)
  if (model %in% c("traab", "tal", "tse")) {
    par <- c(par, diverg = 30)
  }
  if (model == "tse") {
    par <- c(par, pallo = 0.6)
  }

  methods::new("tetmerFit",
    model = model,
    par = par,
    ranges = list(),
    xrange = c(10, 80),
    convergence = 0,
    value = 1,
    k = 21,
    version = as.character(utils::packageVersion("Tetmer"))
  )
}

test_that("expectedSpectrum builds a populated spectrum from model and parameters", {
  spec <- expectedSpectrum(
    "d",
    par = c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200),
    xrange = c(10, 15),
    name = "synthetic",
    k = 21
  )

  expect_s4_class(spec, "spectrum")
  expect_identical(spec@name, "synthetic")
  expect_identical(spec@k, 21)
  expect_identical(spec@data$mult, 10:15)
  expect_length(spec@fits, 0)
  expect_true(all(spec@data$count > 0))
})

test_that("expectedSpectrum uses fit xrange and k by default", {
  fit <- makeExpectedSpectrumTestFit("tal")

  spec <- expectedSpectrum(fit)

  expect_s4_class(spec, "spectrum")
  expect_identical(spec@name, fit@model)
  expect_identical(spec@k, fit@k)
  expect_equal(range(spec@data$mult), fit@xrange)
  expect_true(all(spec@data$count > 0))
})

test_that("expectedSpectrum allows overriding xrange for tetmerFit input", {
  fit <- makeExpectedSpectrumTestFit("d")

  spec <- expectedSpectrum(fit, xrange = c(5, 8), name = "override")

  expect_identical(spec@name, "override")
  expect_identical(spec@data$mult, 5:8)
})

test_that("expectedSpectrum character method requires xrange", {
  expect_error(
    expectedSpectrum("d", par = c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200)),
    "`xrange` must be supplied"
  )
})
