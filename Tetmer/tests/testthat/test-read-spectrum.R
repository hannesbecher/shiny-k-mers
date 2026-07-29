library(Tetmer)

makeRoundTripFit <- function() {
  methods::new("tetmerFit",
    fitType = "auto",
    model = "d",
    par = c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200),
    ranges = list(
      kcov = c(10, 50),
      vf = c(1, 5),
      log10theta = c(-3, 0),
      log10Mbp = c(0, 3)
    ),
    xrange = c(5, 80),
    convergence = 0,
    value = 123.4,
    k = 21,
    version = "2.3.2"
  )
}

test_that("read.spectrum restores fits written by write.spectrum", {
  spec <- methods::new("spectrum",
    name = "Round trip",
    data = data.frame(mult = 1:3, count = c(10, 20, 30)),
    k = 21,
    fits = list(makeRoundTripFit())
  )
  tf <- tempfile(fileext = ".tsp")
  on.exit(unlink(tf))

  write.spectrum(spec, tf)
  header <- readLines(tf, n = 5)
  lines <- readLines(tf)
  restored <- read.spectrum(tf)

  expect_true(all(startsWith(header, "#TETMER")))
  expect_true(any(grepl("#TETMER fit.1.fitType: auto", lines, fixed = TRUE)))
  expect_true(any(grepl("#TETMER fit.1.range.log10theta:", lines, fixed = TRUE)))
  expect_identical(restored@name, spec@name)
  expect_identical(restored@k, spec@k)
  expect_length(restored@fits, 1)
  expect_s4_class(restored@fits[[1]], "tetmerFit")
  expect_identical(restored@fits[[1]]@fitType, spec@fits[[1]]@fitType)
  expect_identical(restored@fits[[1]]@model, spec@fits[[1]]@model)
  expect_equal(restored@fits[[1]]@par, spec@fits[[1]]@par)
  expect_equal(restored@fits[[1]]@ranges, spec@fits[[1]]@ranges)
  expect_equal(restored@fits[[1]]@xrange, spec@fits[[1]]@xrange)
  expect_equal(restored@fits[[1]]@convergence, spec@fits[[1]]@convergence)
  expect_equal(restored@fits[[1]]@value, spec@fits[[1]]@value)
  expect_equal(restored@fits[[1]]@k, spec@fits[[1]]@k)
  expect_identical(restored@fits[[1]]@version, spec@fits[[1]]@version)
})

test_that("read.spectrum round-trips manual fits", {
  manualFit <- methods::new("tetmerFit",
    fitType = "manual",
    model = "d",
    par = c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200),
    ranges = list(),
    xrange = c(1, 80),
    convergence = NA_real_,
    value = NA_real_,
    k = 21,
    version = "2.3.2"
  )
  spec <- methods::new("spectrum",
    name = "Manual round trip",
    data = data.frame(mult = 1:3, count = c(10, 20, 30)),
    k = 21,
    fits = list(manualFit)
  )
  tf <- tempfile(fileext = ".tsp")
  on.exit(unlink(tf))

  write.spectrum(spec, tf)
  lines <- readLines(tf)
  restored <- read.spectrum(tf)

  expect_true(any(grepl("#TETMER fit.1.fitType: manual", lines, fixed = TRUE)))
  expect_identical(restored@fits[[1]]@fitType, "manual")
  expect_equal(restored@fits[[1]]@ranges, list())
  expect_true(is.na(restored@fits[[1]]@convergence))
  expect_true(is.na(restored@fits[[1]]@value))
})

test_that("read.spectrum arguments override Tetmer header name and k", {
  spec <- methods::new("spectrum",
    name = "Original",
    data = data.frame(mult = 1:3, count = c(10, 20, 30)),
    k = 21,
    fits = list(makeRoundTripFit())
  )
  tf <- tempfile(fileext = ".tsp")
  on.exit(unlink(tf))

  write.spectrum(spec, tf)
  restored <- read.spectrum(tf, nam = "Override", k = 31)

  expect_identical(restored@name, "Override")
  expect_identical(restored@k, 31)
  expect_length(restored@fits, 1)
})

test_that("read.spectrum ignores unprefixed comment metadata", {
  tf <- tempfile(fileext = ".hist")
  on.exit(unlink(tf))
  writeLines(c(
    "# Tetmer spectrum file -- written by Tetmer v2.3.2",
    "# name: External header",
    "# k: 99",
    "# fits: 1",
    "1 10",
    "2 20"
  ), tf)

  restored <- read.spectrum(tf, no0 = TRUE)

  expect_identical(restored@name, "MySpectrum")
  expect_identical(restored@k, 0)
  expect_length(restored@fits, 0)
})
