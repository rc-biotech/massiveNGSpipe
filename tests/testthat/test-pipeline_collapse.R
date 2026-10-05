# Unit tests for R/pipeline_collapse.R -- memory_safe_worker_count()'s
# pure logic, and pipeline_collapse()'s per-sample resume behavior.
# ORFik::collapse.fastq() itself (real fastq I/O) is mocked throughout;
# these tests cover massiveNGSpipe's own orchestration.

test_that("memory_safe_worker_count defers to max_workers when everything fits comfortably in memory", {
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Memory_Total_GB = 64, Memory_Usage_GB = 10)
  )
  files <- tempfile(); writeLines("x", files) # tiny file, far under any memory budget
  expect_equal(memory_safe_worker_count(files, max_workers = 8), 8)
})

test_that("memory_safe_worker_count caps below max_workers when memory is the limiting factor", {
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Memory_Total_GB = 10, Memory_Usage_GB = 0)
  )
  testthat::local_mocked_bindings(
    file.size = function(f) rep(3e9, length(f)), # 3GB "files" -- cumsum: 3,6,9,12,15 GB
    .package = "base"
  )
  # free memory = 10GB; cumsum first exceeds it at the 4th file (12GB > 10GB);
  # memory_derived = 4 - 2 (default safety margin) = 2.
  expect_equal(memory_safe_worker_count(rep("f", 5), max_workers = 16), 2)
})

test_that("memory_safe_worker_count never returns less than 1, even under severe memory pressure", {
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Memory_Total_GB = 1, Memory_Usage_GB = 0.9)
  )
  testthat::local_mocked_bindings(
    file.size = function(f) rep(5e9, length(f)), .package = "base"
  )
  expect_equal(memory_safe_worker_count(rep("f", 3), max_workers = 16), 1)
})

test_that("memory_safe_worker_count treats a zero-length files argument as unconstrained", {
  expect_equal(memory_safe_worker_count(character(), max_workers = 8), 8)
})

test_that("memory_safe_worker_count's safety_margin_workers is configurable", {
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Memory_Total_GB = 10, Memory_Usage_GB = 0)
  )
  testthat::local_mocked_bindings(
    file.size = function(f) rep(3e9, length(f)), .package = "base"
  )
  expect_equal(memory_safe_worker_count(rep("f", 5), max_workers = 16, safety_margin_workers = 0), 4)
})

test_that("pipeline_collapse() skips a sample already marked done, collapsing only the rest", {
  # run_experiment_subprocess() genuinely spawns a real callr subprocess
  # (see test-pipeline_async_exec.R) -- a mocked collapse.fastq() in
  # THIS process would not be visible there. Mocking
  # run_experiment_subprocess() itself to a synchronous in-process
  # pass-through tests pipeline_collapse()'s own orchestration (resume
  # filtering, args construction, flag setting) without re-testing
  # run_experiment_subprocess()'s own subprocess mechanics, which
  # test-pipeline_async_exec.R already covers.
  config <- fake_config(preset = "Ribo-seq")
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(
    bam_dir = bam_dir,
    runs = data.table::data.table(Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp_name <- "PRJNA000001-homo_sapiens"
  fake_mark_all_done(config, "trim", exp_name) # prerequisite step, so step_is_next_not_done("collapsed") is TRUE
  trimmed_dir <- file.path(bam_dir, "trim")
  dir.create(trimmed_dir, recursive = TRUE)
  file.create(file.path(trimmed_dir, "trimmed_SRR001.fastq"))
  file.create(file.path(trimmed_dir, "trimmed_SRR002.fastq"))

  # SRR001 already marked done for "collapsed" -- must be skipped.
  set_sample_flag(config, "collapsed", exp_name, "SRR001")

  collapsed_runs <- character()
  testthat::local_mocked_bindings(
    run_experiment_subprocess = function(func, args, ...) do.call(func, args)
  )
  testthat::local_mocked_bindings(
    collapse.fastq = function(filename, outdir, ...) {
      collapsed_runs <<- c(collapsed_runs, basename(filename))
      invisible(NULL)
    }, .package = "ORFik"
  )

  pipeline_collapse(pipelines[["PRJNA000001"]], config)

  expect_identical(collapsed_runs, "trimmed_SRR002.fastq")
  expect_setequal(samples_done(config, "collapsed", exp_name), c("SRR001", "SRR002"))
  expect_true(step_is_done(config, "collapsed", exp_name))
})

test_that("pipeline_collapse() does nothing (no error) when every sample is already done", {
  config <- fake_config(preset = "Ribo-seq")
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(
    bam_dir = bam_dir,
    runs = data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp_name <- "PRJNA000001-homo_sapiens"
  fake_mark_all_done(config, "trim", exp_name)
  trimmed_dir <- file.path(bam_dir, "trim")
  dir.create(trimmed_dir, recursive = TRUE)
  file.create(file.path(trimmed_dir, "trimmed_SRR001.fastq"))
  set_sample_flag(config, "collapsed", exp_name, "SRR001")

  called <- FALSE
  testthat::local_mocked_bindings(
    run_experiment_subprocess = function(func, args, ...) do.call(func, args)
  )
  testthat::local_mocked_bindings(
    collapse.fastq = function(...) { called <<- TRUE }, .package = "ORFik"
  )
  expect_no_error(pipeline_collapse(pipelines[["PRJNA000001"]], config))
  expect_false(called)
  expect_true(step_is_done(config, "collapsed", exp_name))
})
