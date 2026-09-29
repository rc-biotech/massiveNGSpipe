# Unit tests for R/pipeline_run.R -- pipelines_names() (pure) and
# parallel_wrap() (the per-stage-group poll loop) exercised with a
# trivial stub function_call, avoiding any real pipeline work.

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
