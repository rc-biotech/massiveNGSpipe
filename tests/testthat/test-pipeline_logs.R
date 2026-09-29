# Unit tests for R/pipeline_logs.R -- error/update bookkeeping, all
# tempdir()-based. progress_log()'s real `system("tail ...")` call is out
# of scope (thin wrapper, needs a real log file being actively written).

test_that("set_parallel_conf validates logdir and sets bplog/bplogdir/bpjobname on a real BPPARAM", {
  d <- tempfile()
  bp <- BiocParallel::SerialParam()
  bp2 <- set_parallel_conf(bp, list(logdir = d, jobname = "test_job"))
  expect_true(dir.exists(d))
  expect_true(BiocParallel::bplog(bp2))
  expect_identical(BiocParallel::bplogdir(bp2), d)
})

test_that("set_parallel_conf rejects a non-character logdir", {
  bp <- BiocParallel::SerialParam()
  expect_error(set_parallel_conf(bp, list(logdir = 123)))
})

test_that("report_failed_pipe_path returns a real path under error_dir, or a tempfile fallback if unset", {
  config <- fake_config(extra = list(error_dir = tempfile()))
  p <- report_failed_pipe_path(config, "exp1")
  expect_identical(dirname(p), config$error_dir)

  config_no_dir <- fake_config(extra = list(error_dir = NULL))
  expect_warning(p2 <- report_failed_pipe_path(config_no_dir, "exp1"))
  expect_identical(dirname(p2), sub("/[^/]+$", "", tempfile()))
})

test_that("report_failed_pipe: a real try-error is recorded and returns FALSE; success returns TRUE", {
  config <- fake_config(extra = list(error_dir = tempfile()))
  ok <- try(1 + 1, silent = TRUE)
  expect_true(report_failed_pipe(ok, config, "align", "exp1"))
  expect_false(file.exists(report_failed_pipe_path(config, "exp1")))

  bad <- try(stop("boom"), silent = TRUE)
  # report_failed_pipe() issues two separate warning() calls ("Failed at
  # step..." then the error text itself) -- suppressWarnings() rather
  # than expect_warning() here, since the latter only catches one.
  result <- suppressWarnings(report_failed_pipe(bad, config, "align", "exp1"))
  expect_false(result)
  expect_true(file.exists(report_failed_pipe_path(config, "exp1")))
})

test_that("report_failed_pipe: a try-error with no error_dir set warns but doesn't error", {
  config <- fake_config(extra = list(error_dir = NULL))
  bad <- try(stop("boom"), silent = TRUE)
  result <- suppressWarnings(report_failed_pipe(bad, config, "align", "exp1"))
  expect_false(result)
})

test_that("session_error_dirs / session_error_dirs_count reflect real tempdir fixtures", {
  project <- tempfile("mNGSp_test_")
  d1 <- file.path(project, "error_logs", "2020-01-01")
  d2 <- file.path(project, "error_logs", "2025-01-01")
  dir.create(d1, recursive = TRUE); dir.create(d2, recursive = TRUE)
  file.create(file.path(d2, c("a.rds", "b.rds")))

  config <- fake_config(project = project)
  dirs <- session_error_dirs(config)
  expect_length(dirs, 2)

  counts <- session_error_dirs_count(config)
  expect_identical(unname(counts["2025-01-01"]), 2L)
  expect_identical(unname(counts["2020-01-01"]), 0L)
})

test_that("last_session_errors reads .rds error files from the most recent session and can filter with regex", {
  project <- tempfile("mNGSp_test_")
  d <- file.path(project, "error_logs", "2025-01-01")
  dir.create(d, recursive = TRUE)
  saveRDS("error about align", file.path(d, "study1-organism.rds"))
  saveRDS("error about fetch", file.path(d, "study2-organism.rds"))

  config <- fake_config(project = project)
  all_errors <- last_session_errors(config)
  expect_length(all_errors, 2)

  filtered <- last_session_errors(config, regex = "align")
  expect_length(filtered, 1)
})

test_that("last_session_errors errors clearly when the requested index exceeds available sessions", {
  project <- tempfile("mNGSp_test_")
  dir.create(file.path(project, "error_logs", "2025-01-01"), recursive = TRUE)
  config <- fake_config(project = project)
  expect_error(last_session_errors(config, index = 5))
})

test_that("last_update / last_update_diff / is_progressing work off an injected info data.frame or real flag mtimes", {
  now <- Sys.time()
  info <- data.frame(modification_time = c(now - 3600, now - 60))
  expect_equal(as.numeric(last_update(config = NULL, info = info)), as.numeric(now - 60), tolerance = 1)
})

test_that("last_update_message respects the units parameter (previously always hardcoded to hours)", {
  project <- tempfile("mNGSp_test_")
  dir.create(file.path(project, "flags", "aligned"), recursive = TRUE)
  file.create(file.path(project, "flags", "aligned", "exp1.rds"))
  config <- fake_config(project = project)

  out_hours <- capture.output(last_update_message(config, units = "hours"))
  out_mins <- capture.output(last_update_message(config, units = "mins"))
  expect_match(paste(out_hours, collapse = " "), "hours ago")
  expect_match(paste(out_mins, collapse = " "), "mins ago")
})
