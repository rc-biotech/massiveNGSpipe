# Unit tests for R/pipeline_checklist_rates.R -- the live
# processing-speed figures added to checklist.txt (fetch MB/s, STAR's
# own M reads/hr from Log.progress.out, generic samples/studies per
# hour), built from this session's own STAR-index-race incident (a
# ~70-minute silent stall with no visible symptom until investigated by
# hand) showing the cost of a done/total count with no way to tell
# whether it's actually still moving.

test_that("rate_since_first_seen() returns NA on the first observation of a key, then a real rate on the second", {
  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  expect_true(is.na(rate_since_first_seen("k1", 10, now = t0)))
  t1 <- t0 + 3600 # 1 hour later
  expect_equal(rate_since_first_seen("k1", 30, now = t1), 20) # +20 units in 1 hour
})

test_that("rate_since_first_seen() clamps a decrease to 0, never negative", {
  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  rate_since_first_seen("k2", 100, now = t0)
  t1 <- t0 + 1800
  expect_equal(rate_since_first_seen("k2", 40, now = t1), 0)
})

test_that("rate_since_first_seen() returns NA when elapsed time is <= 0 (e.g. a clock that didn't advance)", {
  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  rate_since_first_seen("k3", 1, now = t0)
  expect_true(is.na(rate_since_first_seen("k3", 2, now = t0))) # same instant again
})

test_that("rate_since_first_seen() starts fresh (NA) for a brand new key, independent of other keys", {
  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  rate_since_first_seen("k4a", 50, now = t0)
  expect_true(is.na(rate_since_first_seen("k4b", 999, now = t0))) # different key, never seen before
})

test_that("experiment_run_ids() and experiment_conf() match experiment_sample_counts()'s own iteration", {
  pipelines <- fake_pipelines(
    runs = data.table::data.table(Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  runs <- experiment_run_ids(pipelines)
  expect_identical(runs[["PRJNA000001-homo_sapiens"]], c("SRR001", "SRR002"))

  conf <- experiment_conf(pipelines, "PRJNA000001-homo_sapiens")
  expect_identical(unname(conf["exp"]), "PRJNA000001-homo_sapiens")
  expect_null(experiment_conf(pipelines, "no_such_experiment"))
})

test_that("active_run_id() picks the first not-yet-done run, and NA when every run is already done", {
  config <- fake_config()
  pipelines <- fake_pipelines(
    runs = data.table::data.table(Run = c("SRR001", "SRR002", "SRR003"), LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp <- "PRJNA000001-homo_sapiens"
  set_sample_flag(config, "fetch", exp, "SRR001")

  expect_identical(active_run_id(pipelines, config, "fetch", exp), "SRR002")

  set_sample_flag(config, "fetch", exp, "SRR002")
  set_sample_flag(config, "fetch", exp, "SRR003")
  expect_true(is.na(active_run_id(pipelines, config, "fetch", exp)))

  expect_true(is.na(active_run_id(pipelines, config, "fetch", "no_such_experiment")))
})

test_that("fetch_progress_rate() computes MB/s from real bytes-on-disk growth between two calls", {
  dir <- tempfile("fastq_"); dir.create(dir)
  f <- file.path(dir, "SRR24190359.fastq")
  writeLines(strrep("A", 1e6), f) # ~1MB

  # Must use a fresh key (run id) -- rate_since_first_seen()'s state is
  # shared/global across this whole test file.
  run <- "SRR_fetch_rate_test"
  f2 <- file.path(dir, paste0(run, ".fastq"))
  file.rename(f, f2)

  expect_true(is.na(fetch_progress_rate(dir, run))) # first observation

  writeLines(c(strrep("A", 1e6), strrep("A", 1e6)), f2) # now ~2MB
  rate <- fetch_progress_rate(dir, run)
  expect_true(!is.na(rate) && rate > 0)
})

test_that("fetch_progress_rate() returns NA when nothing matching the run exists yet", {
  dir <- tempfile("fastq_empty_"); dir.create(dir)
  expect_true(is.na(fetch_progress_rate(dir, "SRR_never_started")))
})

test_that("align_progress_rate() parses the last data row's Speed column from a real-shaped Log.progress.out", {
  bam_dir <- tempfile("bam_")
  logs_dir <- file.path(bam_dir, "aligned", "LOGS_SINGLE")
  dir.create(logs_dir, recursive = TRUE)
  run <- "SRR12285169"
  writeLines(c(
    "           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped",
    "                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other",
    "Jun 01 14:39:10   1264.9    21081071       35    53.4%     27.0     0.7%    27.4%     2.8%     0.0%    16.5%     0.0%",
    "Jun 01 14:43:03    978.8    79661848       35    55.4%     27.0     0.8%    23.6%     2.1%     0.0%    18.9%     0.0%"
  ), file.path(logs_dir, paste0("collapsed_trimmed_", run, "_Log.progress.out")))

  expect_equal(align_progress_rate(bam_dir, run), 978.8)
})

test_that("align_progress_rate() strips a trailing 'ALL DONE!' line before reading the last data row", {
  bam_dir <- tempfile("bam_")
  logs_dir <- file.path(bam_dir, "aligned", "LOGS_SINGLE")
  dir.create(logs_dir, recursive = TRUE)
  run <- "SRR_done"
  writeLines(c(
    "           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped",
    "                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other",
    "Jun 01 14:39:10   1264.9    21081071       35    53.4%     27.0     0.7%    27.4%     2.8%     0.0%    16.5%     0.0%",
    "ALL DONE!"
  ), file.path(logs_dir, paste0("collapsed_trimmed_", run, "_Log.progress.out")))

  expect_equal(align_progress_rate(bam_dir, run), 1264.9)
})

test_that("align_progress_rate() returns NA for a header-only log (no data rows yet) or a missing file, without erroring", {
  bam_dir <- tempfile("bam_")
  logs_dir <- file.path(bam_dir, "aligned", "LOGS_SINGLE")
  dir.create(logs_dir, recursive = TRUE)
  writeLines(c(
    "           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped",
    "                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other"
  ), file.path(logs_dir, "collapsed_trimmed_SRR_header_only_Log.progress.out"))

  expect_true(is.na(align_progress_rate(bam_dir, "SRR_header_only")))
  expect_true(is.na(align_progress_rate(bam_dir, "SRR_never_even_started")))
})

test_that("generic_progress_rate() uses a stage+experiment key for the sample case and a stage-only key for the study case", {
  t0 <- as.POSIXct("2026-10-09 10:00:00", tz = "UTC")
  testthat::local_mocked_bindings(Sys.time = function() t0, .package = "base")
  generic_progress_rate("pipe_trim_collapse_test1", "expA", 2)
  generic_progress_rate("pipe_pshift_test1", NA_character_, 1)

  testthat::local_mocked_bindings(Sys.time = function() t0 + 3600, .package = "base")
  rate_sample <- generic_progress_rate("pipe_trim_collapse_test1", "expA", 6)
  rate_study <- generic_progress_rate("pipe_pshift_test1", NA_character_, 3)

  expect_equal(rate_sample, 4)  # +4 samples in 1 hour
  expect_equal(rate_study, 2)   # +2 studies in 1 hour
})
