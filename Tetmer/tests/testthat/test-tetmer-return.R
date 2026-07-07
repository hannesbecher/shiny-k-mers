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