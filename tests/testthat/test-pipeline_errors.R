local_err_config <- function(env = parent.frame()) {
  proj <- withr::local_tempdir(.local_envir = env)
  sess <- "2026-10-17 10:00:00"
  list(project = proj, error_dir = file.path(proj, "error_logs", sess), session_dir = file.path(proj, "session_logs", sess))
}

test_that("error_hint maps known messages and falls back to unknown", {
  expect_match(error_hint("EXITING because of FATAL ERROR: could not open genome file x/genomeParameters.txt"), "STAR index")
  expect_match(error_hint("Library contained no periodic or accepted read-lengths"), "periodicity")
  expect_match(error_hint("fasterq-dump err: storage exhausted while writing"), "disk full")
  expect_match(error_hint("  12 unevaluated and other errors"), "worker died")
  expect_identical(error_hint("something new"), "unknown")
})

test_that("error_one_line keeps the first line and pulls out the first BiocParallel remote error", {
  bp <- "Error : X: P-shifting failed: BiocParallel errors\n  1 remote errors, element index: 1\n  first remote error:\nError in detectRibosomeShifts(rfp): no periodic read-lengths\n"
  out <- error_one_line(bp)
  expect_match(out, "^X: P-shifting failed")
  expect_match(out, "first remote error: no periodic read-lengths")
  expect_identical(error_one_line(""), "(empty error message)")
})

test_that("report_failed_pipe writes the flag file (old format + hint), an event line, and runs the hook", {
  config <- local_err_config()
  seen <- NULL
  config$error_hook <- function(ev) seen <<- ev
  res <- suppressWarnings(report_failed_pipe(try(stop("STAR alignment step failed"), silent = TRUE), config, "align", "PRJ1-x"))
  expect_false(res)
  rec <- readRDS(file.path(config$error_dir, "PRJ1-x.rds"))
  expect_match(rec, "^PRJ1-x \\( align \\) ")                       # old format kept for existing readers
  expect_match(rec, "Likely cause: STAR failed")
  ev <- pipeline_error_events(config)
  expect_identical(nrow(ev), 1L)
  expect_identical(ev$event, "error"); expect_identical(ev$exp, "PRJ1-x"); expect_identical(ev$step, "align")
  expect_identical(seen$exp, "PRJ1-x")
  # the existing checklist reader still parses the record
  em <- experiment_error_message(config, "PRJ1-x")
  expect_identical(em$step, "align")
  expect_true(report_failed_pipe(TRUE, config, "align", "PRJ1-x"))
})

test_that("a failing hook never stops error recording", {
  config <- local_err_config()
  config$error_hook <- function(ev) stop("hook broken")
  expect_warning(error_record_write(config, "PRJ2", "fetch", "boom"), "error_hook failed")
  expect_true(file.exists(file.path(config$error_dir, "PRJ2.rds")))
})

test_that("shell hook from MNGSP_ERROR_HOOK gets the event in environment variables", {
  config <- local_err_config()
  out <- withr::local_tempfile()
  withr::local_envvar(MNGSP_ERROR_HOOK = paste0("echo \"$MNGSP_EVENT|$MNGSP_EXP|$MNGSP_STEP\" > ", shQuote(out)))
  error_record_write(config, "PRJ3", "pshift", "x")
  expect_identical(readLines(out), "error|PRJ3|pshift")
})

test_that("checklist error line counts current flag files and drops cleared ones", {
  config <- local_err_config()
  expect_identical(pipeline_error_summary_line(config), "")
  suppressWarnings({
    report_failed_pipe(try(stop("a"), silent = TRUE), config, "align", "PRJ1")
    report_failed_pipe(try(stop("b"), silent = TRUE), config, "pshift", "PRJ2")
  })
  line <- pipeline_error_summary_line(config)
  expect_match(line, "^ERRORS this session: 2: ")
  expect_match(line, "PRJ1 \\(align\\)"); expect_match(line, "PRJ2 \\(pshift\\)")
  expect_identical(nrow(pipeline_errors(config)), 2L)
  clear_pipeline_error(config, "PRJ1", note = "index rebuilt")
  expect_false(file.exists(file.path(config$error_dir, "PRJ1.rds")))
  expect_length(list.files(file.path(config$error_dir, "resolved")), 1)
  expect_match(pipeline_error_summary_line(config), "^ERRORS this session: 1: PRJ2")
  expect_false(all_pending_failed(config, "PRJ1-homo_sapiens"))          # stage loops retry the cleared study
  ev <- pipeline_error_events(config)
  expect_identical(ev$event, c("error", "error", "resolved"))
  expect_error(clear_pipeline_error(config, "nope"), "No error record")
})

test_that("the resolved/ folder is not counted as an error by session helpers", {
  config <- local_err_config()
  suppressWarnings(report_failed_pipe(try(stop("a"), silent = TRUE), config, "align", "PRJ1"))
  clear_pipeline_error(config, "PRJ1")
  expect_identical(unname(session_error_dirs_count(config)), 0L)
  expect_length(last_session_errors(config), 0)
})

test_that("pipeline_crash_record writes a crash record with call stack and a crash event", {
  config <- local_err_config()
  path <- NULL
  expect_message(tryCatch(withCallingHandlers(stop("cannot allocate vector of size 5 Gb"),
    error = function(e) path <<- pipeline_crash_record(config, "pipe_align_clean", e, sys.calls())),
    error = function(e) NULL), "unexpected")
  rec <- readRDS(path)
  expect_identical(rec$where, "pipe_align_clean")
  expect_match(rec$hint, "out of memory")
  expect_true(length(rec$calls) > 0)
  expect_match(pipeline_error_summary_line(config), "1 stage crash record")
  expect_identical(pipeline_error_events(config)$event, "crash")
})

test_that("parallel_wrap records an escaped stage error before it propagates", {
  config <- local_err_config()
  local_mocked_bindings(all_substeps_done_all = function(...) FALSE, stop_requested = function(...) FALSE,
                        all_pending_failed = function(...) FALSE, pending_exps = function(...) "E1",
                        mark_stage_exited = function(...) NULL)
  expect_error(suppressMessages(parallel_wrap(function(p, c) stop("subscript out of bounds"), list(), config,
                                              steps = "x", wait = 0, stage_name = "pipe_x")), "subscript")
  expect_true(file.exists(file.path(config$session_dir, "crash_pipe_x.rds")))
})

test_that("error records include console log error lines and use them for the hint", {
  config <- local_err_config()
  dir.create(file.path(config$session_dir, "console", "align"), recursive = TRUE)
  writeLines(c("starting", "EXITING because of FATAL ERROR: could not open genome file /x/genomeParameters.txt", "done"),
             file.path(config$session_dir, "console", "align", "PRJ9-saccharomyces_x.err.log"))
  suppressWarnings(report_failed_pipe(try(stop("STAR alignment step failed, see error above"), silent = TRUE),
                                      config, "align", "PRJ9"))
  rec <- readRDS(file.path(config$error_dir, "PRJ9.rds"))
  expect_match(rec, "Console log excerpt:\nEXITING because of FATAL ERROR: could not open genome file")
  expect_match(rec, "Likely cause: STAR index missing")
  expect_identical(console_error_excerpt(config, "NOPE"), character())
})

test_that("metadata_legacy_columns renames UNIQUE to COMMENT once and keeps values", {
  dt <- data.table::data.table(Run = c("A", "B"), UNIQUE = c("", "duplicate of A"))
  metadata_legacy_columns(dt)
  expect_identical(names(dt), c("Run", "COMMENT"))
  expect_identical(dt$COMMENT, c("", "duplicate of A"))
  both <- data.table::data.table(Run = "A", UNIQUE = "x", COMMENT = "y")
  metadata_legacy_columns(both)                     # already has COMMENT: untouched
  expect_identical(names(both), c("Run", "UNIQUE", "COMMENT"))
})
