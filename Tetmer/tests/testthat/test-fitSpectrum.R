library(Tetmer)

test_that("fitSpectrum returns a tetmerFit with ranges from its public arguments", {
  fit <- fitSpectrum(
    E030,
    model = "d",
    kcov = c(5, 100),
    vf = c(1, 30),
    log10theta = c(-3, 0),
    log10Mbp = c(0, 3),
    xrange = c(45, 200),
    maxAttempts = 1
  )

  expect_s4_class(fit, "tetmerFit")
  expect_identical(fit@fitType, "auto")
  expect_identical(fit@model, "d")
  expect_equal(fit@ranges$kcov, c(5, 100))
  expect_equal(fit@ranges$vf, c(1, 30))
  expect_equal(fit@ranges$log10theta, c(-3, 0))
  expect_equal(fit@ranges$log10Mbp, c(0, 3))
  expect_equal(fit@xrange, c(45, 200))
  expect_named(fit@par, c("cov", "vf", "theta", "haplSize"))
  expect_match(paste(capture.output(show(fit)), collapse = "\n"), "Fit type: auto")
  expect_match(paste(capture.output(show(fit)), collapse = "\n"), "log10 theta")
  expect_match(paste(capture.output(show(fit)), collapse = "\n"), "log10 GS \\(Mbp\\)")
})
