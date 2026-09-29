# Unit tests for R/pipeline_abort.R -- the cooperative stop-signal and
# graceful drain-and-abort mechanism. The heavy end (real subprocess
# kill via a live process group) is intentionally out of scope here --
# these test the file-based signal/marker plumbing and the pure
# early-return branches of abort_session_when_next_study_done_per_step(),
# which is what a fixture-based end-to-end run already validated live
# this session.

test_that("stop_requested/request_stop round-trip via a real session_dir", {
  config <- fake_config(session_dir = tempfile("session_"))
  expect_false(stop_requested(config))
  request_stop(config)
  expect_true(stop_requested(config))
})

test_that("mark_stage_exited writes a marker recording stopped vs not, and no-ops for NULL stage_name", {
  config <- fake_config(session_dir = tempfile("session_"))
  mark_stage_exited(config, "pipe_trim_collapse", stopped = TRUE)
  mark_stage_exited(config, "pipe_align_clean", stopped = FALSE)

  d <- file.path(config$session_dir, "stage_stopped")
  expect_true(readRDS(file.path(d, "pipe_trim_collapse.rds")))
  expect_false(readRDS(file.path(d, "pipe_align_clean.rds")))

  mark_stage_exited(config, NULL, stopped = TRUE)
  expect_identical(list.files(d), c("pipe_align_clean.rds", "pipe_trim_collapse.rds"))
})

test_that("abort_session_when_next_study_done_per_step errors for a nonexistent session dir", {
  config <- fake_config()
  expect_error(
    abort_session_when_next_study_done_per_step(config, "does-not-exist"),
    "No such session directory"
  )
})

test_that("abort_session_when_next_study_done_per_step returns FALSE immediately for an already-finished session", {
  project <- tempfile("mNGSp_test_")
  session <- "2026-01-01 00:00:00"
  session_dir <- file.path(project, "session_logs", session)
  dir.create(session_dir, recursive = TRUE)
  saveRDS(list(pipeline_names = "x", pgid = 1L, init_time = Sys.time(), status = "Completed"),
         file.path(session_dir, "session_info.rds"))

  config <- fake_config(project = project)
  expect_message(
    result <- abort_session_when_next_study_done_per_step(config, session),
    "already finished"
  )
  expect_false(result)
  # A finished session should never have a stop request written against it.
  expect_false(file.exists(file.path(session_dir, "stop_requested.rds")))
})

test_that("abort_session_when_next_study_done_per_step drains and returns TRUE once every stage-group has exited", {
  project <- tempfile("mNGSp_test_")
  session <- "2026-01-01 00:00:01"
  session_dir <- file.path(project, "session_logs", session)
  dir.create(session_dir, recursive = TRUE)
  saveRDS(list(pipeline_names = "x", pgid = 999999999L, init_time = Sys.time(), status = "started"),
         file.path(session_dir, "session_info.rds"))

  config <- fake_config(project = project)
  # Full replacement, not modifyList() -- flag_steps is itself a named
  # list, and modifyList() merges nested lists recursively rather than
  # replacing them outright.
  config$flag_steps <- list(pipe_a = "x", pipe_b = "y")

  # Pre-populate both stage markers as already-exited, as if a fast
  # SerialParam run had already finished draining by the time we poll.
  d <- file.path(session_dir, "stage_stopped")
  dir.create(d, recursive = TRUE)
  saveRDS(TRUE, file.path(d, "pipe_a.rds"))
  saveRDS(FALSE, file.path(d, "pipe_b.rds"))

  result <- suppressMessages(
    abort_session_when_next_study_done_per_step(config, session, poll_interval = 1, timeout = 5)
  )
  expect_true(result)
  expect_true(file.exists(file.path(session_dir, "stop_requested.rds")))
  info <- readRDS(file.path(session_dir, "session_info.rds"))
  expect_identical(info$status, "Aborted")
})

test_that("abort_session_when_next_study_done_per_step times out cleanly if a stage never drains", {
  project <- tempfile("mNGSp_test_")
  session <- "2026-01-01 00:00:02"
  session_dir <- file.path(project, "session_logs", session)
  dir.create(session_dir, recursive = TRUE)
  saveRDS(list(pipeline_names = "x", pgid = 999999998L, init_time = Sys.time(), status = "started"),
         file.path(session_dir, "session_info.rds"))

  config <- fake_config(project = project)
  config$flag_steps <- list(pipe_a = "x", pipe_never_exits = "y")

  expect_warning(
    result <- abort_session_when_next_study_done_per_step(config, session, poll_interval = 1, timeout = 1),
    "Timed out waiting"
  )
  expect_false(result)
  info <- readRDS(file.path(session_dir, "session_info.rds"))
  expect_identical(info$status, "started") # left untouched, not falsely marked Aborted
})
