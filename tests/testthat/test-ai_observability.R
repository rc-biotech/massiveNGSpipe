# Unit tests for R/ai_observability.R -- the "what's going on right
# now" helpers built from this session's own repeated need to
# manually re-derive this picture by hand (ps aux, tailing
# session_logs/error_logs, /proc/<pid>/wchan) while debugging a live
# hang (2026-10-07/08).

test_that("studies_currently_processing() returns the empty shape, not an error, when ps finds nothing relevant", {
  testthat::local_mocked_bindings(
    system2 = function(...) character(), .package = "base"
  )
  result <- studies_currently_processing()
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 0)
  expect_identical(names(result), c("pid", "accession", "cmd_snippet"))
})

test_that("studies_currently_processing() extracts an accession from a matching process line, and NA when none is present", {
  testthat::local_mocked_bindings(
    system2 = function(...) c(
      "PID COMMAND",
      "  1234 Rscript resume_PRJNA637713.R",
      "  5678 /usr/bin/fastp --in1 foo.fastq",
      "  9999 some_totally_unrelated_process"
    ),
    .package = "base"
  )
  result <- studies_currently_processing()
  expect_equal(nrow(result), 2) # the unrelated process must NOT match
  expect_identical(result[pid == 1234]$accession, "PRJNA637713")
  expect_true(is.na(result[pid == 5678]$accession)) # fastp matches by keyword, no accession visible
})

test_that("checklist_health() reports no_checklist when the file doesn't exist yet", {
  config <- fake_config(session_dir = tempfile("session_"))
  health <- checklist_health(config)
  expect_identical(health$status, "no_checklist")
  expect_false(health$exists)
})

test_that("checklist_health() recognizes a cleanly labeled title regardless of staleness", {
  session_dir <- tempfile("session_")
  dir.create(session_dir, recursive = TRUE)
  path <- file.path(session_dir, "checklist.txt")
  writeLines(c("Pipeline (done after 1.2 hours)", "some body text"), path)
  Sys.setFileTime(path, Sys.time() - 60 * 60) # 1 hour old -- would be "stale" if unlabeled

  config <- fake_config(session_dir = session_dir)
  health <- checklist_health(config, stale_after_mins = 30)

  expect_identical(health$status, "done")
  expect_true(health$labeled)
})

test_that("checklist_health() flags an unlabeled, old checklist as possibly_hung", {
  session_dir <- tempfile("session_")
  dir.create(session_dir, recursive = TRUE)
  path <- file.path(session_dir, "checklist.txt")
  writeLines(c("Pipeline (running for 13.0 hours)", "some body text"), path)
  Sys.setFileTime(path, Sys.time() - 60 * 60) # 1 hour since last update, no final label

  config <- fake_config(session_dir = session_dir)
  health <- checklist_health(config, stale_after_mins = 30)

  expect_identical(health$status, "possibly_hung")
  expect_false(health$labeled)
})

test_that("checklist_health() treats an unlabeled but recently-updated checklist as still running", {
  session_dir <- tempfile("session_")
  dir.create(session_dir, recursive = TRUE)
  path <- file.path(session_dir, "checklist.txt")
  writeLines("Pipeline (running for 0.1 hours)", path)

  config <- fake_config(session_dir = session_dir)
  health <- checklist_health(config, stale_after_mins = 30)

  expect_identical(health$status, "running")
})

test_that("summarize_last_errors() flags a known install-race message as likely_transient, and a real error as not", {
  project <- tempfile("mNGSp_test_")
  d <- file.path(project, "error_logs", "2025-01-01")
  dir.create(d, recursive = TRUE)
  saveRDS("exp1 ( align ) cannot open file '/path/massiveNGSpipe.rdb': No such file or directory",
         file.path(d, "exp1.rds"))
  saveRDS("exp2 ( trim ) some genuinely real error message", file.path(d, "exp2.rds"))

  config <- fake_config(project = project)
  summary <- summarize_last_errors(config)

  expect_equal(nrow(summary), 2)
  expect_true(summary[grepl("exp1", error_key, fixed = TRUE)]$likely_transient)
  expect_false(summary[grepl("exp2", error_key, fixed = TRUE)]$likely_transient)
})

test_that("summarize_last_errors() returns the empty shape when there are no recorded errors", {
  project <- tempfile("mNGSp_test_")
  dir.create(file.path(project, "error_logs", "2025-01-01"), recursive = TRUE)
  config <- fake_config(project = project)

  summary <- summarize_last_errors(config)
  expect_equal(nrow(summary), 0)
  expect_identical(names(summary), c("error_key", "likely_transient", "message"))
})

test_that("ai_context_snapshot() returns every expected piece without erroring, even with nothing going on", {
  config <- fake_config(session_dir = tempfile("session_"))
  testthat::local_mocked_bindings(system2 = function(...) character(), .package = "base")

  snap <- ai_context_snapshot(config)

  expect_true(is.data.frame(snap$active_processes))
  expect_identical(snap$checklist$status, "no_checklist")
  expect_true(is.data.frame(snap$recent_errors))
  expect_identical(snap$threads, config$threads)
  expect_identical(snap$mode, config$mode)
  expect_identical(snap$project, config$project)
})

test_that("stage_triage() reflects a step's own done/sample markers and filters errors to the right experiment", {
  config <- fake_config(preset = "Ribo-seq")
  exp <- "study1-organism"
  set_sample_flag(config, "aligned", exp, "SRR001")
  set_flag(config, "trim", exp)

  d <- file.path(config$project, "error_logs", "2025-01-01")
  dir.create(d, recursive = TRUE)
  saveRDS("study1-organism ( aligned ) boom", file.path(d, "study1-organism.rds"))
  saveRDS("study2-organism ( aligned ) unrelated boom", file.path(d, "study2-organism.rds"))

  testthat::local_mocked_bindings(system2 = function(...) character(), .package = "base")

  triage <- stage_triage(exp, "trim", config)
  expect_true(triage$step_done)
  expect_false(triage$currently_processing)
  expect_length(triage$recent_errors, 1)
  expect_match(unname(triage$recent_errors), "study1-organism", fixed = TRUE)

  aligned_triage <- stage_triage(exp, "aligned", config)
  expect_identical(aligned_triage$samples_done, "SRR001")
  expect_identical(aligned_triage$n_samples_done, 1L)
})

test_that("pipeline_status_all() returns one row per project dir, including ones that never ran", {
  project1 <- tempfile("mNGSp_test_")
  session1 <- file.path(project1, "session_logs", "2025-01-01")
  dir.create(session1, recursive = TRUE)
  writeLines("Pipeline (done after 2.0 hours)", file.path(session1, "checklist.txt"))

  project2 <- tempfile("mNGSp_test_") # never run -- no session_logs at all

  testthat::local_mocked_bindings(system2 = function(...) character(), .package = "base")

  status <- pipeline_status_all(c(project1, project2))
  expect_equal(nrow(status), 2)
  expect_identical(status[project_dir == project1]$status, "done")
  expect_identical(status[project_dir == project2]$status, "no_checklist")
})
