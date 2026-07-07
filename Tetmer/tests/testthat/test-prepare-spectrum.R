library(Tetmer)

test_that("prepareSpectrum fills missing multiplicities with zero counts", {
  raw <- methods::new("spectrum",
    name = "gappy",
    data = data.frame(mult = c(1, 3, 5), count = c(10, 30, 50)),
    k = 21,
    fits = list()
  )

  prepared <- Tetmer:::prepareSpectrum(raw)

  expect_s4_class(prepared, "spectrum")
  expect_identical(names(prepared@data), c("mult", "count"))
  expect_equal(prepared@data$mult[1:5], 1:5)
  expect_equal(prepared@data$count[1:5], c(10, 0, 30, 0, 50))
  expect_equal(nrow(prepared@data), 500)
})
