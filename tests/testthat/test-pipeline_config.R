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

test_that("bpparam_from_config() uses config$threads[[step]] by default", {
  config <- fake_config()
  bp <- bpparam_from_config(config, "default")
  expect_s4_class(bp, "SerialParam") # fake_config()'s thread_type default
})

test_that("bpparam_from_config() uses an explicit workers override instead of config$threads[[step]], when given", {
  # Real use case: memory_safe_worker_count() computes a dynamic worker
  # count at runtime (e.g. for pipeline_collapse()) that can be lower
  # than the static config value -- bpparam_from_config() must actually
  # build the BPPARAM from THAT number, not silently fall back to the
  # config default.
  config <- fake_config(extra = list(thread_type = BiocParallel::MulticoreParam))
  bp <- bpparam_from_config(config, "collapse", workers = 3)
  expect_equal(BiocParallel::bpnworkers(bp), 3)
})

test_that("bpparam_from_config() errors if the requested step has no entry in config$threads", {
  config <- fake_config()
  expect_error(bpparam_from_config(config, "not_a_real_step"))
})
