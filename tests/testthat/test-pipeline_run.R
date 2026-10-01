# Unit tests for R/pipeline_run.R -- pipelines_names() (pure),
# parallel_wrap() (the per-stage-group poll loop), and run_pipeline()'s
# final-status checklist write, all exercised with trivial stub
# pipeline_steps functions, avoiding any real pipeline work.

#' Minimal config able to actually run run_pipeline() end to end with
#' SerialParam + stub steps -- fake_config() alone doesn't set
#' thread_type/parallel_conf/threads$main/pipeline_steps, since most
#' tests only need parallel_wrap() or pipeline_checklist() in isolation.
#' @param stubs a list, same length as config$flag_steps, of
#' function(pipelines, config) stage-group bodies
fake_runnable_config <- function(stubs, project = tempfile("proj_")) {
  config <- fake_config(project = project)
  config$thread_type <- BiocParallel::SerialParam
  config$parallel_conf <- list(log = FALSE, logdir = NA_character_,
                               jobname = "test", stop.on.error = TRUE)
  config$threads$main <- 1
  config$pipeline_steps <- stats::setNames(stubs, names(config$flag_steps))
  config
}

#' Newest session_logs/<init_time>/checklist.txt's first line, after a
#' run_pipeline() call (its own session_dir is computed fresh internally
#' from config$project + a new init_time, so the caller's own
#' config$session_dir going in is irrelevant/overwritten).
newest_checklist_title <- function(project) {
  newest <- sort(list.dirs(file.path(project, "session_logs"), recursive = FALSE),
                 decreasing = TRUE)[1]
  readLines(file.path(newest, "checklist.txt"))[1]
}

test_that("pipelines_names flattens experiment names; recursive=FALSE keeps nested structure", {
  pipelines <- fake_pipelines()
  expect_identical(pipelines_names(pipelines), "PRJNA000001-homo_sapiens")

  nested <- pipelines_names(pipelines, recursive = FALSE)
  expect_true(is.list(nested))
})

test_that("parallel_wrap exits immediately (no stub call) if steps are already all done", {
  config <- fake_config(session_dir = tempfile("session_"))
  pipelines <- fake_pipelines()
  exp <- pipelines_names(pipelines)
  steps <- config$flag_steps$pipe_align_clean
  for (s in steps) set_flag(config, s, exp)

  call_count <- 0
  stub <- function(pipelines, config) call_count <<- call_count + 1

  expect_message(
    parallel_wrap(stub, pipelines, config, steps, wait = 5, stage_name = "pipe_align_clean"),
    "Done for step pipeline"
  )
  expect_equal(call_count, 0)
  marker <- readRDS(file.path(config$session_dir, "stage_stopped", "pipe_align_clean.rds"))
  expect_false(marker) # exited because done, not because stopped
})

test_that("parallel_wrap calls function_call until steps become done, then marks exited (not stopped)", {
  config <- fake_config(session_dir = tempfile("session_"))
  pipelines <- fake_pipelines()
  exp <- pipelines_names(pipelines)
  steps <- config$flag_steps$pipe_align_clean

  stub <- function(pipelines, config) for (s in steps) set_flag(config, s, exp)

  expect_message(
    parallel_wrap(stub, pipelines, config, steps, wait = 5, stage_name = "pipe_align_clean"),
    "Done for step pipeline"
  )
  marker <- readRDS(file.path(config$session_dir, "stage_stopped", "pipe_align_clean.rds"))
  expect_false(marker)
})

test_that("parallel_wrap honors a pre-existing stop request without calling function_call even once", {
  config <- fake_config(session_dir = tempfile("session_"))
  pipelines <- fake_pipelines()
  request_stop(config)

  call_count <- 0
  stub <- function(pipelines, config) call_count <<- call_count + 1

  expect_message(
    parallel_wrap(stub, pipelines, config, config$flag_steps$pipe_align_clean,
                  wait = 5, stage_name = "pipe_align_clean"),
    "Stopped \\(graceful shutdown requested\\)"
  )
  expect_equal(call_count, 0)
  marker <- readRDS(file.path(config$session_dir, "stage_stopped", "pipe_align_clean.rds"))
  expect_true(marker) # exited because stopped
})

test_that("run_pipeline() writes '(done after N hours)' when every step completes normally", {
  project <- tempfile("proj_")
  pipelines <- fake_pipelines()
  n <- length(fake_config(project = project)$flag_steps)
  config <- fake_runnable_config(rep(list(function(pipelines, config) NULL), n), project)
  exp <- pipelines_names(pipelines)
  for (s in names(config$flag)) set_flag(config, s, exp) # already done -> each stage exits immediately

  invisible(capture.output(suppressMessages(run_pipeline(pipelines, config, wait = 1))))
  expect_match(newest_checklist_title(project), "\\(done after [0-9.]+ hours\\)$")
})

test_that("run_pipeline() writes '(stopped gracefully after N hours)' on a deliberate request_stop()", {
  project <- tempfile("proj_")
  pipelines <- fake_pipelines()
  n <- length(fake_config(project = project)$flag_steps)
  stubs <- rep(list(function(pipelines, config) request_stop(config)), n)
  config <- fake_runnable_config(stubs, project)

  invisible(capture.output(suppressMessages(run_pipeline(pipelines, config, wait = 1))))
  expect_match(newest_checklist_title(project), "\\(stopped gracefully after [0-9.]+ hours\\)$")
})

test_that("run_pipeline() writes '(aborted after N hours)' and still re-throws on an uncaught stage error", {
  # The core "abort catch" property: a real error must both (a) leave a
  # labeled trace in checklist.txt, via on.exit(), and (b) still
  # propagate to the caller -- never silently swallowed.
  project <- tempfile("proj_")
  pipelines <- fake_pipelines()
  n <- length(fake_config(project = project)$flag_steps)
  stubs <- rep(list(function(pipelines, config) NULL), n)
  stubs[[1]] <- function(pipelines, config) stop("simulated crash")
  config <- fake_runnable_config(stubs, project)

  expect_error(
    invisible(capture.output(suppressWarnings(suppressMessages(run_pipeline(pipelines, config, wait = 1))))),
    "simulated crash"
  )
  expect_match(newest_checklist_title(project), "\\(aborted after [0-9.]+ hours\\)$")
})
