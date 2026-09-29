# Unit tests for R/pipeline_config.R -- config validation helpers not
# already covered by test-pipeline_flags.R (flag_grouping/preset_grouping
# are exercised there since they're really a pipeline_flags.R concept).

test_that("is_config accepts a config-shaped list and rejects anything else", {
  expect_true(is_config(list(preset = "Ribo-seq")))
  expect_false(is_config(list(no_preset_here = TRUE)))
  expect_false(is_config("not a list"))
  expect_false(is_config(NULL))
})

test_that("get_fun_name finds the exported name of a value in a namespace", {
  expect_identical(get_fun_name(sum, "base"), "sum")
})

test_that("get_fun_name returns character(0) when the value isn't in that namespace", {
  not_in_base <- function(x) x
  expect_identical(get_fun_name(not_in_base, "base"), character(0))
})
