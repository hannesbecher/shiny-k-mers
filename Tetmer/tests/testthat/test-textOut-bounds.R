# tests/testthat/test-textOut-bounds.R
# Tests for HTML output that highlights parameters at optimisation bounds

library(Tetmer)

makeAutoInput_simple <- function(mod = "d") {
  list(
    fitmod  = "auto",
    mod     = mod,
    akcov   = c(10, 100),
    avf     = c(1, 10),
    ath     = c(-2, 0.6),
    ayadj   = c(0, 3),
    axrange = c(45, 200)
  )
}

makeTinySpectrum <- function(k = 21) {
  methods::new(
    "spectrum",
    name = "x",
    data = data.frame(mult = 1, count = 1),
    k = k,
    fits = list()
  )
}

test_that("textOutHtml highlights parameters at bounds", {
  input  <- makeAutoInput_simple("d")
  spec   <- makeTinySpectrum()
  bounds <- Tetmer:::buildOptimisationSpec(input)

  optimised <- list(par = c(
    cov      = unname(bounds$lower["cov"]),
    vf       = 5,
    theta    = 0.01,
    haplSize = 100
  ))

  html <- as.character(Tetmer:::textOutHtml(input, optimised, spec))

  expect_match(html, "monoploid k-mer cov: <span")
  expect_match(html, "title=\"at lower bound\"")
})

test_that("textOutHtml does not add highlight markup when no parameters are at bounds", {
  input <- makeAutoInput_simple("d")
  spec  <- makeTinySpectrum()

  optimised <- list(par = c(
    cov      = 55,
    vf       = 5,
    theta    = 0.02,
    haplSize = 100
  ))

  html <- as.character(Tetmer:::textOutHtml(input, optimised, spec))
  expect_false(grepl("<span", html, fixed = TRUE))
})
