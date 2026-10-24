# pipe_fetch() must re-check the download caps before every study, not only once per call.
test_that("pipe_fetch stops downloading when the unprocessed-download cap is reached mid-loop", {
  pipelines <- lapply(paste0("PRJ", 1:5), function(a) list(accession = a))
  config <- list(stop_downloading_new_data_at_drive_usage = 92, max_unprocessed_downloads = 2,
                 project = tempfile())
  downloaded <- character()
  backlog <- 0L
  local_mocked_bindings(
    get_system_usage = function(...) list(Drive_Usage_Percent = "50%"),
    unprocessed_downloads_count = function(pipelines, config) backlog,
    report_failed_pipe_path = function(config, exp) tempfile(),
    stop_requested = function(config) FALSE,
    report_failed_pipe = function(try, config, step, exp) TRUE,
    pipeline_download = function(pipeline, config, pipelines) {
      downloaded <<- c(downloaded, pipeline$accession); backlog <<- backlog + 1L })
  suppressMessages(pipe_fetch(pipelines, config))
  expect_equal(downloaded, c("PRJ1", "PRJ2", "PRJ3"))   # backlog 3 > cap 2 after the third study
})

test_that("pipe_fetch stops when the drive cap is reached mid-loop", {
  pipelines <- lapply(paste0("PRJ", 1:4), function(a) list(accession = a))
  config <- list(stop_downloading_new_data_at_drive_usage = 92, max_unprocessed_downloads = 30, project = tempfile())
  usage <- 80; downloaded <- character()
  local_mocked_bindings(
    get_system_usage = function(...) list(Drive_Usage_Percent = paste0(usage, "%")),
    unprocessed_downloads_count = function(pipelines, config) 0L,
    report_failed_pipe_path = function(config, exp) tempfile(),
    stop_requested = function(config) FALSE,
    report_failed_pipe = function(try, config, step, exp) TRUE,
    pipeline_download = function(pipeline, config, pipelines) {
      downloaded <<- c(downloaded, pipeline$accession); usage <<- usage + 6 })
  suppressMessages(pipe_fetch(pipelines, config))
  expect_equal(downloaded, c("PRJ1", "PRJ2"))           # 80 -> 86 -> 92 (>= cap)
})
