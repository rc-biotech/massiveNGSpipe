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

test_that("per-organism stages honour a study-level error record", {
  d <- tempfile(); dir.create(d); config <- list(error_dir = d, flag = c(fetch = 1, trim = 2))
  pipelines <- list(list(accession = "PRJ1", organisms = list(`Homo sapiens` = list(conf = c(exp = "PRJ1-homo_sapiens")),
                                                                `Mus musculus` = list(conf = c(exp = "PRJ1-mus_musculus")))),
                    list(accession = "PRJ2", organisms = list(`Homo sapiens` = list(conf = c(exp = "PRJ2-homo_sapiens")))))
  acc <- exp_accessions(pipelines)
  expect_equal(unname(acc[c("PRJ1-mus_musculus", "PRJ2-homo_sapiens")]), c("PRJ1", "PRJ2"))
  saveRDS("fetch error", file.path(d, "PRJ1.rds"))
  expect_true(exp_or_study_failed(config, "PRJ1-homo_sapiens", acc[["PRJ1-homo_sapiens"]]))
  expect_false(exp_or_study_failed(config, "PRJ2-homo_sapiens", acc[["PRJ2-homo_sapiens"]]))
  saveRDS("align error", file.path(d, "PRJ2-homo_sapiens.rds"))
  expect_true(exp_or_study_failed(config, "PRJ2-homo_sapiens", "PRJ2"))
  # failed studies sitting at "fetch" do not count toward the download backlog
  local_mocked_bindings(progress_report = function(...) c(`PRJ1-homo_sapiens` = 1, `PRJ1-mus_musculus` = 1,
                                                          `PRJ2-homo_sapiens` = 1))
  expect_equal(unprocessed_downloads_count(pipelines, config), 0L)
  unlink(file.path(d, c("PRJ1.rds", "PRJ2-homo_sapiens.rds")))
  expect_equal(unprocessed_downloads_count(pipelines, config), 3L)
})
