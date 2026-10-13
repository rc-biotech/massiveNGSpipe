test_that("pshift extra-length hook is a no-op unless opted in, and propagates failures", {
  expect_identical(pshift_extra_lengths_hook(NULL, 1, NULL, NULL), 1)            # default config: unchanged
  err <- simpleError("shift failed")
  expect_identical(pshift_extra_lengths_hook(NULL, err, 22:24, NULL), err)       # shift error passes through
  testthat::local_mocked_bindings(pshift_add_extra_lengths = function(...) stop("boom"), name = function(x) "exp")
  r <- suppressMessages(pshift_extra_lengths_hook(list(), 1, 22:24, NULL))
  expect_true(inherits(r, "error"))                                              # flag not set -> retried
  testthat::local_mocked_bindings(pshift_add_extra_lengths = function(...) invisible(NULL))
  expect_identical(pshift_extra_lengths_hook(list(), 1, 22:24, NULL), 1)
})

test_that("extend_shift_table keeps existing rows and only adds lengths not already present", {
  tab <- list(a.ofst = data.frame(fraction = c(28, 29), offsets_start = c(-12, -12)))
  testthat::local_mocked_bindings(
    detectRibosomeShifts = function(reads, df, accepted.lengths, ...) data.frame(fraction = c(23, 29, 34), offsets_start = c(-9, -13, -14)),
    fimport = function(x) NULL, .package = "ORFik")
  r <- extend_shift_table(data.frame(x = 1), tab, extra_lengths = c(23, 29, 34))
  expect_identical(r$a.ofst$fraction, c(23, 28, 29, 34))
  expect_identical(r$a.ofst$offsets_start[r$a.ofst$fraction == 29], -12)          # existing offset kept
  expect_identical(nrow(attr(r, "added")), 2L)
})
