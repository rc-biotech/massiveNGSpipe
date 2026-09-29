# Unit tests for R/pipeline_async_exec.R -- the log-sanitizer is pure and
# tempfile()-based. run_experiment_subprocess() genuinely spawns a real
# callr subprocess, but the functions it runs here are trivial/instant
# (no STAR/fastp/genome), so this stays fast while still exercising the
# real callr plumbing (success, error-reraise, and log-file creation).

test_that("sanitize_progress_log converts bare \\r into real newlines", {
  f <- tempfile()
  writeLines_raw <- function(txt) writeChar(txt, f, eos = NULL)
  writeLines_raw("a\rb\rc")
  sanitize_progress_log(f)
  expect_identical(readChar(f, file.info(f)$size), "a\nb\nc")
})

test_that("sanitize_progress_log collapses \\r\\n first, so no double newlines appear", {
  # A file with ONLY proper \r\n sequences and no bare \r is intentionally
  # left untouched entirely (see the function's own grepl("\r(?!\n)", ...)
  # guard) -- so this needs at least one bare \r alongside the \r\n
  # sequences to actually exercise both gsub passes in one call, matching
  # a realistic mixed log (normal lines plus one progress-bar update).
  f <- tempfile()
  writeChar("a\r\nb\rc", f, eos = NULL)
  sanitize_progress_log(f)
  expect_identical(readChar(f, file.info(f)$size), "a\nb\nc")
})

test_that("sanitize_progress_log leaves a file with only proper \\r\\n untouched (no bare \\r present)", {
  f <- tempfile()
  writeChar("a\r\nb\r\nc", f, eos = NULL)
  before <- readChar(f, file.info(f)$size)
  sanitize_progress_log(f)
  expect_identical(readChar(f, file.info(f)$size), before)
})

test_that("sanitize_progress_log is a no-op for a missing, empty, or already-clean file", {
  missing <- tempfile()
  expect_false(file.exists(missing))
  expect_no_error(sanitize_progress_log(missing))
  expect_false(file.exists(missing)) # doesn't create it either

  empty <- tempfile(); file.create(empty)
  sanitize_progress_log(empty)
  expect_identical(file.info(empty)$size, 0)

  clean <- tempfile(); writeLines(c("a", "b"), clean)
  before <- readChar(clean, file.info(clean)$size)
  sanitize_progress_log(clean)
  expect_identical(readChar(clean, file.info(clean)$size), before)
})

test_that("run_experiment_subprocess requires distinct stdout/stderr paths", {
  same <- tempfile()
  expect_error(
    run_experiment_subprocess(function() 1, logfile_out = same, logfile_err = same),
    "logfile_out != logfile_err"
  )
})

test_that("run_experiment_subprocess runs a trivial function for real and returns its value", {
  out <- tempfile(); err <- tempfile()
  result <- run_experiment_subprocess(
    function(x) x + 1, args = list(x = 41),
    logfile_out = out, logfile_err = err, poll_interval = 0.2
  )
  expect_identical(result, 42)
  expect_true(file.exists(out))
  expect_true(file.exists(err))
})

test_that("run_experiment_subprocess re-raises an error from inside the subprocess", {
  out <- tempfile(); err <- tempfile()
  expect_error(
    run_experiment_subprocess(
      function() stop("boom from child"),
      logfile_out = out, logfile_err = err, poll_interval = 0.2
    ),
    "boom from child"
  )
})

test_that("run_experiment_subprocess calls on_poll at least once for a function that outlives one poll interval", {
  out <- tempfile(); err <- tempfile()
  poll_env <- new.env()
  poll_env$n <- 0
  run_experiment_subprocess(
    function() Sys.sleep(0.5),
    logfile_out = out, logfile_err = err,
    poll_interval = 0.1,
    on_poll = function() poll_env$n <- poll_env$n + 1
  )
  expect_gte(poll_env$n, 1)
})
