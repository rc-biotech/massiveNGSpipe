# Unit tests for R/pipeline_checklist.R -- the Nextflow-style progress
# checklist, its formatting helpers, and the system-usage summary line.
# get_system_usage()/detect_drive() (from ORFik, real system calls) are
# mocked via local_mocked_bindings() so these stay fast and deterministic
# regardless of the machine running the tests.

test_that("stage_marker_step maps known stages and falls through to NA otherwise", {
  expect_identical(stage_marker_step("pipe_fetch"), "fetch")
  expect_identical(stage_marker_step("pipe_trim_collapse"), "trim")
  expect_identical(stage_marker_step("pipe_align_clean"), "aligned")
  expect_identical(stage_marker_step("pipe_exp_ofst"), "ofst")
  expect_identical(stage_marker_step("pipe_convert"), "bigwig")
  expect_identical(stage_marker_step("totally_unknown_stage"), NA_character_)
})

test_that("experiment_sample_counts counts rows per organism from a fake pipelines list", {
  pipelines <- fake_pipelines(runs = data.table::data.table(
    Run = c("SRR001", "SRR002", "SRR003"),
    ScientificName = "Homo sapiens"
  ))
  counts <- experiment_sample_counts(pipelines)
  expect_identical(unname(counts["PRJNA000001-homo_sapiens"]), 3L)
})

test_that("pipeline_log_base prefers session_dir when set, else falls back to project/log_pipeline", {
  config_with_session <- fake_config(session_dir = "/tmp/some/session")
  expect_identical(pipeline_log_base(config_with_session), "/tmp/some/session")

  config_no_session <- fake_config(session_dir = NULL)
  expect_identical(pipeline_log_base(config_no_session),
                   file.path(config_no_session$project, "log_pipeline"))
})

test_that("checklist_path builds off pipeline_log_base", {
  config <- fake_config(session_dir = "/tmp/some/session")
  expect_identical(checklist_path(config), "/tmp/some/session/checklist.txt")
})

test_that("format_checklist renders the expected marker/detail per state", {
  tab <- data.table::data.table(
    stage = c("pipe_a", "pipe_b", "pipe_c"),
    done = c(2L, 1L, 0L), total = c(2L, 2L, 2L),
    state = c("done", "running", "queued"),
    active_experiment = c(NA_character_, "expX", NA_character_),
    active_done = c(NA_integer_, 3L, NA_integer_),
    active_total = c(NA_integer_, 5L, NA_integer_)
  )
  txt <- format_checklist(tab)
  lines <- strsplit(txt, "\n")[[1]]
  expect_match(lines[1], "^\\[✔\\] pipe_a", perl = TRUE)
  expect_match(lines[2], "active: expX \\(3/5 samples\\)")
  expect_match(lines[3], "^\\[-\\] pipe_c")
})

test_that("format_checklist: a running stage with no active_experiment shows the no-detail fallback", {
  tab <- data.table::data.table(
    stage = "pipe_a", done = 0L, total = 2L, state = "running",
    active_experiment = NA_character_, active_done = NA_integer_, active_total = NA_integer_
  )
  expect_match(format_checklist(tab), "no per-sample detail available")
})

test_that("format_system_usage_line reports the cap note only when at/above the cap", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1.2, Memory_Usage_Percent = 3.4,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "50%")
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92))
  u <- format_system_usage_line(config)
  expect_identical(u[["line"]], "CPU (1.2%), Memory (3.4%), Drive /dev/fake (50%%)")
  expect_identical(u[["cap_note"]], "")
})

test_that("format_system_usage_line: cap note fires at exactly the boundary (>=, not >)", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "92%")
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92))
  u <- format_system_usage_line(config)
  expect_match(u[["cap_note"]], "DRIVE AT/ABOVE 92% CAP")
  expect_match(u[["cap_note"]], "pipe_fetch\\(\\) is pausing new downloads")
})

test_that("format_system_usage_line: unparseable drive percent never errors, cap note stays empty", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) NA_character_,
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = NA_character_, Drive_Usage_Percent = NA_character_)
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92))
  expect_no_error(u <- format_system_usage_line(config))
  expect_identical(u[["cap_note"]], "")
})

test_that("pipeline_checklist writes a fixed-position header (usage/cap-note/2 blanks) before the stage table", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "99%")
  )
  session_dir <- tempfile("session_")
  config <- fake_config(session_dir = session_dir,
                        extra = list(stop_downloading_new_data_at_drive_usage = 92))
  pipelines <- fake_pipelines()

  tab <- suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  expect_s3_class(tab, "data.table")

  lines <- readLines(checklist_path(config))
  expect_match(lines[1], "^Pipeline status as of")
  expect_match(lines[2], "^CPU \\(1%\\)")
  expect_match(lines[3], "DRIVE AT/ABOVE 92% CAP")
  expect_identical(lines[4], "")
  expect_identical(lines[5], "")
  expect_match(lines[6], "^\\[.\\] ") # stage table always starts on line 6
})

test_that("pipeline_checklist: an experiment with zero started samples stays 'queued', not 'running'", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(session_dir = tempfile("session_"))
  pipelines <- fake_pipelines()
  tab <- suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  first_stage_with_marker <- tab[stage == "pipe_align_clean"]
  expect_identical(first_stage_with_marker$state, "queued")
})
