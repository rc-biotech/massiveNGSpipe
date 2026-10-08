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

test_that("pipeline_config() caps pshifted/valid_pshift/pcounts at threads_blas_cap, but leaves trim/collapse's own smaller caps untouched", {
  # Confirmed live, 2026-10-08: running threads_default (46) forked
  # workers for one of these steps, itself inside an already-forked
  # main-level stage-group worker, multiplies out to hundreds of OS
  # processes+threads -- each ALSO starting its own uncapped OpenBLAS
  # thread pool -- and can exhaust the container's cgroup pids.max
  # outright (pthread_create() -> EAGAIN). A real-data benchmark
  # (shift_qc_cached(), 30 samples) also found 32 workers fastest in
  # practice, not just safest -- see threads_blas_cap's own roxygen.
  testthat::local_mocked_bindings(
    blas_set_num_threads = function(...) NULL,
    omp_set_num_threads = function(...) NULL,
    .package = "RhpcBLASctl"
  )
  config <- pipeline_config(project_dir = tempfile("mNGSp_test_"), mode = "local",
                            preset = "empty", verbose = FALSE, google_url = NULL,
                            threads_default = 46)
  expect_identical(config$threads$pshifted, 32)
  expect_identical(config$threads$valid_pshift, 32)
  expect_identical(config$threads$pcounts, 32)
  expect_identical(config$threads$trim, 8)       # unaffected
  expect_identical(config$threads$collapse, 16)  # unaffected
  expect_identical(config$threads$default, 46)   # unaffected
})

test_that("pipeline_config() never caps below the real threads_default, when that's already smaller than threads_blas_cap", {
  testthat::local_mocked_bindings(
    blas_set_num_threads = function(...) NULL,
    omp_set_num_threads = function(...) NULL,
    .package = "RhpcBLASctl"
  )
  config <- pipeline_config(project_dir = tempfile("mNGSp_test_"), mode = "local",
                            preset = "empty", verbose = FALSE, google_url = NULL,
                            threads_default = 4)
  expect_identical(config$threads$valid_pshift, 4)
})

test_that("pipeline_config() sets BLAS/OpenMP thread count to threads_blas_cap once, at config-creation time", {
  blas_calls <- list(); omp_calls <- list()
  testthat::local_mocked_bindings(
    blas_set_num_threads = function(threads) blas_calls[[length(blas_calls) + 1]] <<- threads,
    omp_set_num_threads = function(threads) omp_calls[[length(omp_calls) + 1]] <<- threads,
    .package = "RhpcBLASctl"
  )
  pipeline_config(project_dir = tempfile("mNGSp_test_"), mode = "local",
                  preset = "empty", verbose = FALSE, google_url = NULL, threads_blas_cap = 32)
  expect_equal(blas_calls, list(32))
  expect_equal(omp_calls, list(32))
})
