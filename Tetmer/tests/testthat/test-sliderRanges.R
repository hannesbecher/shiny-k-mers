library(Tetmer)

test_that("sliderRanges returns the expected default structure", {
  old_opt <- getOption("Tetmer.sliderRanges")
  on.exit(options(Tetmer.sliderRanges = old_opt), add = TRUE)
  options(Tetmer.sliderRanges = NULL)

  ranges <- sliderRanges()

  expected_names <- c(
    "gsMin", "gsMax",
    "kcovMin", "kcovMax",
    "vfMin", "vfMax",
    "thMin", "thMax",
    "divMin", "divMax",
    "xrangeMin", "xrangeMax",
    "palloMin", "palloMax",
    "ymax"
  )

  expect_true(is.list(ranges))
  expect_identical(names(ranges), expected_names)
  expect_true(all(vapply(ranges, is.numeric, logical(1))))
  expect_true(all(vapply(ranges, length, integer(1)) == 1L))
})

test_that("setSliderRanges applies partial overrides and preserves defaults", {
  old_opt <- getOption("Tetmer.sliderRanges")
  on.exit(options(Tetmer.sliderRanges = old_opt), add = TRUE)
  options(Tetmer.sliderRanges = NULL)

  setSliderRanges(list(kcovMax = 500, ymax = 3))
  ranges <- sliderRanges()

  expect_equal(ranges$kcovMax, 500)
  expect_equal(ranges$ymax, 3)
  expect_equal(ranges$gsMin, 3)
  expect_equal(ranges$palloMax, 1)
})

test_that("setSliderRanges rejects unknown names", {
  old_opt <- getOption("Tetmer.sliderRanges")
  on.exit(options(Tetmer.sliderRanges = old_opt), add = TRUE)
  options(Tetmer.sliderRanges = NULL)

  expect_error(
    setSliderRanges(list(notARealRange = 123)),
    "Unknown slider range name\\(s\\): notARealRange"
  )
})

test_that("setSliderRanges requires a named list", {
  old_opt <- getOption("Tetmer.sliderRanges")
  on.exit(options(Tetmer.sliderRanges = old_opt), add = TRUE)
  options(Tetmer.sliderRanges = NULL)

  expect_error(setSliderRanges(3), "named list")
  expect_error(setSliderRanges(list(3)), "named list")
})

test_that("setSliderRanges enforces scalar numeric values", {
  old_opt <- getOption("Tetmer.sliderRanges")
  on.exit(options(Tetmer.sliderRanges = old_opt), add = TRUE)
  options(Tetmer.sliderRanges = NULL)

  expect_error(
    setSliderRanges(list(ymax = c(1, 2))),
    "must have length 1"
  )
  expect_error(
    setSliderRanges(list(ymax = "banana")),
    "must be numeric"
  )
})

test_that("sliderRanges validates externally set invalid option state", {
  old_opt <- getOption("Tetmer.sliderRanges")
  on.exit(options(Tetmer.sliderRanges = old_opt), add = TRUE)

  options(Tetmer.sliderRanges = list(ymax = "banana"))
  expect_error(sliderRanges(), "Slider range `ymax` must be numeric")

  options(Tetmer.sliderRanges = list(unknown = 1))
  expect_error(sliderRanges(), "Unknown slider range name\\(s\\): unknown")
})

test_that("setSliderRanges accepts numeric-like values and stores numerics", {
  old_opt <- getOption("Tetmer.sliderRanges")
  on.exit(options(Tetmer.sliderRanges = old_opt), add = TRUE)
  options(Tetmer.sliderRanges = NULL)

  setSliderRanges(list(ymax = "3"))
  ranges <- sliderRanges()

  expect_type(ranges$ymax, "double")
  expect_equal(ranges$ymax, 3)
})
