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

test_that("format_checklist names the active study with its real sample TOTAL (not a fabricated fraction) when active_done is unknown", {
  tab <- data.table::data.table(
    stage = "pipe_c", done = 1L, total = 3L, state = "running",
    active_experiment = "PRJNA000002-homo_sapiens", active_done = NA_integer_, active_total = 5L
  )
  line <- format_checklist(tab)
  expect_match(line, "active: PRJNA000002-homo_sapiens \\(5 samples, no per-sample progress tracked for this stage\\)")
})

test_that("pipeline_checklist names the active study for a MARKERLESS stage (cigar_collapse), with a real sample total but no per-sample done-count", {
  # Built from a real request (Håkon, 2026-10-09): markerless stages
  # (cigar_collapse, merge_study, counts, convert's covRLE phase) used
  # to show just "active (no per-sample detail available for this
  # stage)" with no study name at all -- but the SAME "first
  # not-yet-done candidate" inference already used everywhere else in
  # this checklist (see active_run_id()) is enough to name it.
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(preset = "RNA-seq", mode = "online", session_dir = tempfile("session_"))
  pipelines <- c(
    fake_pipelines(accession = "PRJNA000001", runs = data.table::data.table(
      Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE", LIBRARYTYPE = "RNA",
      ScientificName = "Homo sapiens")),
    fake_pipelines(accession = "PRJNA000002", runs = data.table::data.table(
      Run = c("SRR003", "SRR004", "SRR005"), LibraryLayout = "SINGLE", LIBRARYTYPE = "RNA",
      ScientificName = "Homo sapiens"))
  )
  # PRJNA000001 already finished cigar_collapse (a single, whole-study
  # flag -- no per-sample marker exists for this stage at all);
  # PRJNA000002 hasn't.
  fake_mark_all_done(config, "cigar_collapse", "PRJNA000001-homo_sapiens")

  tab <- suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  row <- tab[stage == "pipe_cigar_collapse"]
  expect_identical(row$state, "running")
  expect_identical(row$active_experiment, "PRJNA000002-homo_sapiens")
  expect_identical(row$active_total, 3L)
  expect_true(is.na(row$active_done))
})

test_that("pipeline_checklist names the active study for a marker-step stage even when NO remaining candidate has per-sample evidence yet", {
  # Built from a real observed case (Håkon, 2026-10-09): pipe_exp_ofst
  # HAS a marker_step ("ofst"), but used to stay "queued"-looking with
  # no active study at all whenever upstream work hadn't fed its next
  # study any ofst samples yet, despite 33/37 studies already being
  # genuinely done for this stage.
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(preset = "RNA-seq", mode = "online", session_dir = tempfile("session_"))
  pipelines <- c(
    fake_pipelines(accession = "PRJNA000001", runs = data.table::data.table(
      Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE", LIBRARYTYPE = "RNA",
      ScientificName = "Homo sapiens")),
    fake_pipelines(accession = "PRJNA000002", runs = data.table::data.table(
      Run = c("SRR003", "SRR004", "SRR005"), LibraryLayout = "SINGLE", LIBRARYTYPE = "RNA",
      ScientificName = "Homo sapiens"))
  )
  # PRJNA000001 fully done for pipe_exp_ofst (both its "exp" and "ofst"
  # flags); PRJNA000002 hasn't had ANY sample's "ofst" marker set yet.
  fake_mark_all_done(config, c("exp", "ofst"), "PRJNA000001-homo_sapiens")

  tab <- suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  row <- tab[stage == "pipe_exp_ofst"]
  expect_identical(row$state, "running") # previously stayed "queued" despite 1/2 studies already done -- the bug this fixes
  expect_identical(row$active_experiment, "PRJNA000002-homo_sapiens")
  expect_identical(row$active_total, 3L)
  expect_true(is.na(row$active_done))
})

test_that("stage_progress_rate_label never seeds generic_progress_rate's baseline with NA active_done (would permanently poison that experiment's rate)", {
  key <- "pipe_na_guard_test PRJNA_na_guard"
  if (exists(key, envir = .rate_state)) rm(list = key, envir = .rate_state)

  label1 <- stage_progress_rate_label("pipe_na_guard_test", "ofst", "running", fake_config(),
                                      active_experiment = "PRJNA_na_guard", active_done = NA_integer_,
                                      n_done = 0L, pipelines = list())
  expect_true(is.na(label1))
  expect_false(exists(key, envir = .rate_state)) # confirms no baseline was seeded by the NA call

  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  testthat::local_mocked_bindings(Sys.time = function() t0, .package = "base")
  label2 <- stage_progress_rate_label("pipe_na_guard_test", "ofst", "running", fake_config(),
                                      active_experiment = "PRJNA_na_guard", active_done = 2L,
                                      n_done = 0L, pipelines = list())
  expect_true(is.na(label2)) # first REAL observation -- seeds the baseline now, nothing to compute a rate from yet

  t1 <- t0 + 3600 # 1 hour later
  testthat::local_mocked_bindings(Sys.time = function() t1, .package = "base")
  label3 <- stage_progress_rate_label("pipe_na_guard_test", "ofst", "running", fake_config(),
                                      active_experiment = "PRJNA_na_guard", active_done = 6L,
                                      n_done = 0L, pipelines = list())
  expect_match(label3, "^4\\.0 samples/hr$") # +4 in 1hr -- NOT poisoned by the earlier NA call
})

test_that("stage_progress_rate_label shows 'waiting for available RAM' instead of the align rate when the active sample is deferred", {
  # See R/pipeline_align_memory.R -- a deferred sample has no STAR
  # progress to report, so this status replaces (never supplements)
  # the usual M reads/hr figure.
  config <- fake_config()
  pipelines <- fake_pipelines(
    runs = data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp <- "PRJNA000001-homo_sapiens"
  mark_waiting_for_ram(config, exp, "SRR001")

  label <- stage_progress_rate_label("pipe_align_clean", "aligned", "running", config,
                                     active_experiment = exp, active_done = 0L, n_done = 0L,
                                     pipelines = pipelines)
  expect_match(label, "^waiting for available RAM, [0-9.]+ min$")
})

test_that("stage_progress_rate_label falls back to the normal align rate once a sample's RAM wait is cleared", {
  config <- fake_config()
  pipelines <- fake_pipelines(
    runs = data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp <- "PRJNA000001-homo_sapiens"
  mark_waiting_for_ram(config, exp, "SRR001")
  clear_waiting_for_ram(config, exp, "SRR001")

  label <- stage_progress_rate_label("pipe_align_clean", "aligned", "running", config,
                                     active_experiment = exp, active_done = 0L, n_done = 0L,
                                     pipelines = pipelines)
  expect_true(is.na(label)) # no Log.progress.out exists in this fixture -> NA, but NOT the waiting text
})

test_that("format_system_usage_line reports the drive cap note only when at/above the cap", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1.2, Memory_Usage_Percent = 3.4,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "50%")
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92,
                                     max_unprocessed_downloads = 30))
  pipelines <- fake_pipelines()
  u <- format_system_usage_line(pipelines, config)
  expect_identical(u[["line"]], "CPU (1.2%), Memory (3.4%), Drive /dev/fake (50%%)")
  expect_identical(u[["drive_cap_note"]], "")
})

test_that("format_system_usage_line: drive cap note fires at exactly the boundary (>=, not >)", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "92%")
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92,
                                     max_unprocessed_downloads = 30))
  pipelines <- fake_pipelines()
  u <- format_system_usage_line(pipelines, config)
  expect_match(u[["drive_cap_note"]], "DRIVE AT/ABOVE 92% CAP")
  expect_match(u[["drive_cap_note"]], "pipe_fetch\\(\\) is pausing new downloads")
})

test_that("format_system_usage_line: unparseable drive percent never errors, drive cap note stays empty", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) NA_character_,
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = NA_character_, Drive_Usage_Percent = NA_character_)
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92,
                                     max_unprocessed_downloads = 30))
  pipelines <- fake_pipelines()
  expect_no_error(u <- format_system_usage_line(pipelines, config))
  expect_identical(u[["drive_cap_note"]], "")
})

test_that("format_system_usage_line: backlog cap note fires when unprocessed downloads exceed the cap", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "10%")
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92,
                                     max_unprocessed_downloads = 30))
  pipelines <- fake_pipelines()
  testthat::local_mocked_bindings(unprocessed_downloads_count = function(...) 31L)

  u <- format_system_usage_line(pipelines, config)
  expect_match(u[["backlog_cap_note"]], "UNPROCESSED DOWNLOADS 31 > 30 CAP")
  expect_match(u[["backlog_cap_note"]], "pipe_fetch\\(\\) is pausing new downloads")
})

test_that("format_system_usage_line: backlog cap note stays empty when under the cap", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "10%"),
    unprocessed_downloads_count = function(...) 5L
  )
  config <- fake_config(extra = list(stop_downloading_new_data_at_drive_usage = 92,
                                     max_unprocessed_downloads = 30))
  pipelines <- fake_pipelines()
  u <- format_system_usage_line(pipelines, config)
  expect_identical(u[["backlog_cap_note"]], "")
})

test_that("unprocessed_downloads_count returns 0 when this config has no fetch step (e.g. local mode)", {
  config <- fake_config(extra = list(flag = c(trim = "trim"))) # no "fetch" name present
  pipelines <- fake_pipelines()
  expect_identical(unprocessed_downloads_count(pipelines, config), 0L)
})

test_that("unprocessed_downloads_count counts experiments whose progress is exactly at the fetch step, from real flags", {
  config <- fake_config(mode = "online") # mode = "online" is what puts "fetch" in config$flag at all
  pipelines <- fake_pipelines(exp_name = "PRJNA000001-homo_sapiens")
  exp <- "PRJNA000001-homo_sapiens"

  # Nothing done yet -- not counted (progress is at 0, not at the fetch index).
  expect_identical(unprocessed_downloads_count(pipelines, config), 0L)

  # start+fetch done, nothing past that -- now counted.
  set_flag(config, "start", exp)
  set_flag(config, "fetch", exp)
  expect_identical(unprocessed_downloads_count(pipelines, config), 1L)
})

test_that("pipeline_checklist writes a fixed-position header (usage/2 cap-notes/2 blanks) before the stage table", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "99%"),
    unprocessed_downloads_count = function(...) 99L
  )
  session_dir <- tempfile("session_")
  config <- fake_config(session_dir = session_dir,
                        extra = list(stop_downloading_new_data_at_drive_usage = 92,
                                    max_unprocessed_downloads = 30))
  pipelines <- fake_pipelines()

  tab <- suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  expect_s3_class(tab, "data.table")

  lines <- readLines(checklist_path(config))
  expect_match(lines[1], "^Pipeline status as of")
  expect_match(lines[2], "^CPU \\(1%\\)")
  expect_match(lines[3], "DRIVE AT/ABOVE 92% CAP")
  expect_match(lines[4], "UNPROCESSED DOWNLOADS 99 > 30 CAP")
  expect_identical(lines[5], "")
  expect_identical(lines[6], "")
  expect_match(lines[7], "^\\[.\\] ") # stage table always starts on line 7
})

test_that("pipeline_checklist: header stays 6 lines with both blank cap-note slots when neither cap fires", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "10%"),
    unprocessed_downloads_count = function(...) 0L
  )
  session_dir <- tempfile("session_")
  config <- fake_config(session_dir = session_dir,
                        extra = list(stop_downloading_new_data_at_drive_usage = 92,
                                    max_unprocessed_downloads = 30))
  pipelines <- fake_pipelines()

  suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  lines <- readLines(checklist_path(config))
  expect_identical(lines[3], "")
  expect_identical(lines[4], "")
  expect_identical(lines[5], "")
  expect_identical(lines[6], "")
  expect_match(lines[7], "^\\[.\\] ") # stage table stays on line 7 regardless
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

test_that("format_elapsed_hours formats to one decimal place", {
  now <- as.POSIXct("2026-10-01 12:00:00", tz = "UTC")
  expect_identical(format_elapsed_hours(now - 3600 * 5.23, now), "5.2 hours")
  expect_identical(format_elapsed_hours(now, now), "0 hours")
})

test_that("pipeline_checklist title line has no fractional seconds", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(session_dir = tempfile("session_"))
  pipelines <- fake_pipelines()
  suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  title <- readLines(checklist_path(config))[1]
  expect_match(title, "^Pipeline status as of \\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}( |$)")
  expect_false(grepl("\\.\\d+", title)) # no decimal/fractional-second remnant
})

test_that("pipeline_checklist: run_status = NULL auto-shows elapsed time from config$init_time", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(session_dir = tempfile("session_"),
                        extra = list(init_time = Sys.time() - 3600 * 5.2))
  pipelines <- fake_pipelines()
  suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  title <- readLines(checklist_path(config))[1]
  expect_match(title, "\\(running for 5\\.2 hours\\)$")
})

test_that("pipeline_checklist: an explicit run_status overrides the auto elapsed-time note", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(session_dir = tempfile("session_"),
                        extra = list(init_time = Sys.time() - 3600 * 5.2))
  pipelines <- fake_pipelines()
  suppressMessages(pipeline_checklist(pipelines, config, print = FALSE, run_status = "aborted after 5.2 hours"))
  title <- readLines(checklist_path(config))[1]
  expect_match(title, "\\(aborted after 5\\.2 hours\\)$")
})

test_that("pipeline_checklist: no parenthetical at all when init_time isn't set and run_status isn't given", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(session_dir = tempfile("session_")) # no init_time
  pipelines <- fake_pipelines()
  suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  title <- readLines(checklist_path(config))[1]
  expect_false(grepl("\\(", title))
})

test_that("pipeline_eta_hours returns 0 once the final stage is already fully done", {
  expect_equal(pipeline_eta_hours("eta_test_done", 5, 5), 0)
})

test_that("pipeline_eta_hours returns NA on the first observation (no rate exists yet)", {
  expect_true(is.na(pipeline_eta_hours("eta_test_first_obs", 1, 10)))
})

test_that("pipeline_eta_hours divides remaining studies by the real observed studies/hr rate", {
  key <- "eta_test_real_rate"
  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  testthat::local_mocked_bindings(Sys.time = function() t0, .package = "base")
  pipeline_eta_hours(key, 2, 10) # seeds the baseline: amount = 2 at t0

  t1 <- t0 + 3600 # 1 hour later
  testthat::local_mocked_bindings(Sys.time = function() t1, .package = "base")
  # +2 studies done in 1 hour -> 2 studies/hr; 10 - 4 = 6 remaining -> 3 hours
  expect_equal(pipeline_eta_hours(key, 4, 10), 3)
})

test_that("pipeline_eta_hours returns NA (not Inf) when the rate is still 0 (no progress since first seen)", {
  key <- "eta_test_zero_rate"
  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  testthat::local_mocked_bindings(Sys.time = function() t0, .package = "base")
  pipeline_eta_hours(key, 2, 10)

  t1 <- t0 + 3600
  testthat::local_mocked_bindings(Sys.time = function() t1, .package = "base")
  expect_true(is.na(pipeline_eta_hours(key, 2, 10))) # no progress -> rate 0 -> NA, never Inf
})

test_that("pipeline_checklist title line shows a real ETA, in hours, once the final stage's own rate is known", {
  # See pipeline_eta_hours() -- built from a real request (Håkon,
  # 2026-10-09) to show a topline ETA alongside "running for N.N hours",
  # computed from the final stage's own real studies/hr throughput
  # rather than a guess.
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  config <- fake_config(session_dir = tempfile("session_"),
                        extra = list(init_time = Sys.time() - 3600))
  pipelines <- c(fake_pipelines(accession = "PRJNA000001"),
                fake_pipelines(accession = "PRJNA000002"),
                fake_pipelines(accession = "PRJNA000003"))

  # .rate_state (R/pipeline_checklist_rates.R) is shared, in-process,
  # global state -- an EARLIER test in this file may already have
  # seeded a "pipe_counts__eta" baseline (e.g. any other
  # pipeline_checklist() call whose fixture also has a "pipe_counts"
  # stage) at the real wall clock time, which would make this test's
  # own "first call seeds the baseline" step below land on an
  # ALREADY-seeded key instead. Clear it first so this test is
  # deterministic regardless of run order.
  if (exists("pipe_counts__eta", envir = .rate_state)) rm("pipe_counts__eta", envir = .rate_state)

  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  testthat::local_mocked_bindings(Sys.time = function() t0, .package = "base")
  fake_mark_all_done(config, "pcounts", "PRJNA000001-homo_sapiens") # 1/3 done -- seeds the rate baseline, no ETA yet
  suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))

  t1 <- t0 + 3600 # 1 hour later
  testthat::local_mocked_bindings(Sys.time = function() t1, .package = "base")
  fake_mark_all_done(config, "pcounts", "PRJNA000002-homo_sapiens") # 2/3 done, +1 study in 1hr -> 1 study/hr; 1 remaining -> ETA 1h
  suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))

  title <- readLines(checklist_path(config))[1]
  expect_match(title, "ETA: 1 hours\\)$")
})

test_that("session_log_dirs / session_checklist_path reflect real tempdir fixtures, newest first", {
  project <- tempfile("mNGSp_test_")
  d1 <- file.path(project, "session_logs", "2020-01-01")
  d2 <- file.path(project, "session_logs", "2025-01-01")
  dir.create(d1, recursive = TRUE); dir.create(d2, recursive = TRUE)
  writeLines("newest", file.path(d2, "checklist.txt"))
  writeLines("oldest", file.path(d1, "checklist.txt"))

  config <- fake_config(project = project)
  dirs <- session_log_dirs(config)
  expect_length(dirs, 2)
  expect_identical(basename(dirs[1]), "2025-01-01") # newest first

  expect_identical(readLines(session_checklist_path(config, index = 1)), "newest")
  expect_identical(readLines(session_checklist_path(config, index = 2)), "oldest")
})

test_that("session_checklist_path errors clearly when the requested index exceeds available sessions", {
  project <- tempfile("mNGSp_test_")
  dir.create(file.path(project, "session_logs", "2025-01-01"), recursive = TRUE)
  config <- fake_config(project = project)
  expect_error(session_checklist_path(config, index = 5), "only 1")
})

test_that("watch_pipeline_checklist(open_in_new_text_window = TRUE) opens an RStudio terminal running watch, instead of looping in this console", {
  project <- tempfile("mNGSp_test_")
  session_dir <- file.path(project, "session_logs", "2025-01-01")
  dir.create(session_dir, recursive = TRUE)
  config <- fake_config(project = project)

  create_calls <- list(); send_calls <- list()
  testthat::local_mocked_bindings(
    isAvailable = function(...) TRUE,
    terminalCreate = function(...) { create_calls[[length(create_calls) + 1]] <<- TRUE; "term-1" },
    terminalSend = function(id, text) send_calls[[length(send_calls) + 1]] <<- list(id = id, text = text),
    .package = "rstudioapi"
  )

  result <- watch_pipeline_checklist(config, interval = 3, open_in_new_text_window = TRUE)

  expect_identical(result, "term-1")
  expect_length(create_calls, 1)
  expect_length(send_calls, 1)
  expect_identical(send_calls[[1]]$id, "term-1")
  expect_match(send_calls[[1]]$text, "watch -n 3 \"cat '")
  expect_match(send_calls[[1]]$text, "checklist.txt")
})

test_that("watch_pipeline_checklist(open_in_new_text_window = TRUE) errors clearly outside RStudio instead of silently doing nothing", {
  project <- tempfile("mNGSp_test_")
  dir.create(file.path(project, "session_logs", "2025-01-01"), recursive = TRUE)
  config <- fake_config(project = project)

  testthat::local_mocked_bindings(isAvailable = function(...) FALSE, .package = "rstudioapi")

  expect_error(
    watch_pipeline_checklist(config, open_in_new_text_window = TRUE),
    "RStudio"
  )
})

test_that("watch_pipeline_checklist() prefers config$session_dir over the newest session_logs entry, when index isn't explicitly given", {
  # The scenario this guards: config$session_dir points at THIS specific
  # session, but a second, newer run_pipeline() call started
  # concurrently elsewhere -- "index = 1 = newest by directory listing"
  # would otherwise silently resolve to that OTHER session instead of
  # the one this config actually belongs to.
  project <- tempfile("mNGSp_test_")
  own_session <- file.path(project, "session_logs", "2020-01-01")
  other_newer_session <- file.path(project, "session_logs", "2025-01-01")
  dir.create(own_session, recursive = TRUE); dir.create(other_newer_session, recursive = TRUE)
  writeLines("own session content", file.path(own_session, "checklist.txt"))
  writeLines("other newer session content", file.path(other_newer_session, "checklist.txt"))

  config <- fake_config(project = project, session_dir = own_session)
  captured <- NULL
  testthat::local_mocked_bindings(
    isAvailable = function(...) TRUE,
    terminalCreate = function(...) "term-1",
    terminalSend = function(id, text) captured <<- text,
    .package = "rstudioapi"
  )

  watch_pipeline_checklist(config, open_in_new_text_window = TRUE)
  expect_match(captured, "2020-01-01", fixed = TRUE)
  expect_false(grepl("2025-01-01", captured, fixed = TRUE))
})

test_that("watch_pipeline_checklist() uses an explicit index even when config$session_dir is set", {
  project <- tempfile("mNGSp_test_")
  own_session <- file.path(project, "session_logs", "2020-01-01")
  other_older_session <- file.path(project, "session_logs", "2019-01-01")
  dir.create(own_session, recursive = TRUE); dir.create(other_older_session, recursive = TRUE)
  writeLines("own", file.path(own_session, "checklist.txt"))
  writeLines("other older", file.path(other_older_session, "checklist.txt"))

  config <- fake_config(project = project, session_dir = own_session)
  captured <- NULL
  testthat::local_mocked_bindings(
    isAvailable = function(...) TRUE,
    terminalCreate = function(...) "term-1",
    terminalSend = function(id, text) captured <<- text,
    .package = "rstudioapi"
  )

  # index = 2 is the OLDER session (2019-01-01), explicitly requested
  # -- must win over config$session_dir (which points at 2020-01-01).
  watch_pipeline_checklist(config, index = 2, open_in_new_text_window = TRUE)
  expect_match(captured, "2019-01-01", fixed = TRUE)
})

test_that("pipeline_checklist() surfaces a live samples/hr rate once an active stage has made progress across two calls", {
  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(CPU_Usage_Percent = 1, Memory_Usage_Percent = 1,
                                          Drive = "/dev/fake", Drive_Usage_Percent = "1%")
  )
  # Ribo-seq specifically: RNA-seq's default stage names don't include
  # "pipe_trim_collapse" at all (no separate fetch/trim stage the same
  # way).
  config <- fake_config(preset = "Ribo-seq", mode = "online", session_dir = tempfile("session_"))
  pipelines <- fake_pipelines(
    runs = data.table::data.table(Run = c("SRR001", "SRR002", "SRR003"), LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp <- "PRJNA000001-homo_sapiens"
  # "trim" (pipe_trim_collapse), not "aligned" (pipe_align_clean) --
  # the align stage routes to the STAR-Log.progress.out-specific rate
  # source instead of the generic samples/hr path (see
  # align_progress_rate(), tested separately), so it would correctly
  # stay NA here with no real log file on disk.
  set_sample_flag(config, "trim", exp, "SRR001")

  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  testthat::local_mocked_bindings(Sys.time = function() t0, .package = "base")
  tab1 <- suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  expect_true(is.na(tab1[stage == "pipe_trim_collapse"]$rate_label)) # first observation, no rate yet

  set_sample_flag(config, "trim", exp, "SRR002")
  testthat::local_mocked_bindings(Sys.time = function() t0 + 3600, .package = "base")
  tab2 <- suppressMessages(pipeline_checklist(pipelines, config, print = FALSE))
  trim_row <- tab2[stage == "pipe_trim_collapse"]
  expect_identical(trim_row$rate_label, "1.0 samples/hr") # +1 sample (1->2 done) in 1 hour

  checklist_txt <- readLines(checklist_path(config))
  expect_true(any(grepl("1.0 samples/hr", checklist_txt, fixed = TRUE)))
})
