# Unit tests for R/ai_trim_fix_workflow.R -- the consolidated
# redetect/fix/rerun/verify workflow built to replace the several
# separate hand-run bash/R snippets this session used for each of the
# real barcode fixes (PRJNA926112, PRJNA770650, PRJEB50305).

test_that("sample_currently_processing() reflects whether a process table entry mentions the experiment", {
  testthat::local_mocked_bindings(
    system2 = function(...) c("12345 fastp --in1 x PRJNA999999-homo_sapiens"), .package = "base"
  )
  expect_true(sample_currently_processing("PRJNA999999-homo_sapiens"))
  expect_false(sample_currently_processing("PRJNA000001-homo_sapiens"))
})

test_that("sample_currently_processing() returns FALSE (not an error) when pgrep finds nothing", {
  testthat::local_mocked_bindings(
    system2 = function(...) { warning("no matches"); character() }, .package = "base"
  )
  expect_false(sample_currently_processing("PRJNA000001-homo_sapiens"))
})

test_that("last_real_touch() ignores known diagnostic-noise files and reports the newest real file", {
  config <- fake_config()
  exp <- "PRJNA000001-homo_sapiens"
  bam_dir <- file.path(config$config["bam"], exp)
  dir.create(file.path(bam_dir, "read_length_distribution"), recursive = TRUE)
  dir.create(file.path(bam_dir, "aligned"), recursive = TRUE)
  # Noise: should never influence the result.
  writeLines("x", file.path(bam_dir, "read_length_distribution", "SRR001.csv"))
  Sys.setFileTime(file.path(bam_dir, "read_length_distribution", "SRR001.csv"), Sys.time())
  # Real file, deliberately old.
  real_file <- file.path(bam_dir, "aligned", "SRR001.bam")
  writeLines("x", real_file)
  old_time <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC")
  Sys.setFileTime(real_file, old_time)

  touch <- last_real_touch(exp, config)
  expect_equal(as.numeric(touch), as.numeric(old_time), tolerance = 2)
})

test_that("last_real_touch() returns NA when the experiment has no directories at all", {
  config <- fake_config()
  expect_true(is.na(last_real_touch("PRJNA_NONEXISTENT-homo_sapiens", config)))
})

test_that("sibling_file_snapshot()/changed_files() detect a modified sibling file and ignore an untouched one", {
  bam_dir <- tempfile("bam_")
  dir.create(bam_dir)
  f1 <- file.path(bam_dir, "SRR001.bam")
  f2 <- file.path(bam_dir, "SRR002.bam")
  writeLines("x", f1); writeLines("x", f2)

  before <- sibling_file_snapshot(bam_dir, c("SRR001", "SRR002"))
  expect_length(before, 2)

  Sys.sleep(1.1) # ensure a detectable mtime difference
  writeLines("y", f2) # simulate an accidental sibling touch

  after <- sibling_file_snapshot(bam_dir, c("SRR001", "SRR002"))
  changed <- changed_files(before, after)
  expect_identical(changed, f2)
})

test_that("sibling_file_snapshot() returns empty for no siblings, without touching the filesystem", {
  expect_length(sibling_file_snapshot("/nonexistent", character()), 0)
})

test_that("barcode_fix_candidates() excludes samples already present in manual_trim_fix_log.csv", {
  config <- fake_config()
  outliers <- data.table::data.table(
    study_accession = c("PRJNA001", "PRJNA002"),
    ScientificName = c("Homo sapiens", "Homo sapiens"),
    raw_library = c("SRR001", "SRR002"),
    flag_reason = "barcode_outlier",
    n_study_samples = c(10, 10),
    n_study_barcode_true = c(9, 9),
    trim_mean_length = c(33, 33),
    study_median_len = c(35, 35)
  )
  outliers_path <- tempfile(fileext = ".csv")
  data.table::fwrite(outliers, outliers_path)

  dir.create(dirname(config$complete_metadata), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(data.table::data.table(exp = "PRJNA001-homo_sapiens", run = "SRR001"),
                     file.path(config$project, "manual_trim_fix_log.csv"))

  testthat::local_mocked_bindings(sample_currently_processing = function(exp) FALSE)
  testthat::local_mocked_bindings(last_real_touch = function(exp, config) as.POSIXct("2020-01-01", tz = "UTC"))

  cands <- barcode_fix_candidates(config, outliers_path)
  expect_identical(cands$run, "SRR002")
})

test_that("barcode_fix_candidates() excludes a study that is currently being processed", {
  config <- fake_config()
  outliers <- data.table::data.table(
    study_accession = "PRJNA003", ScientificName = "Homo sapiens", raw_library = "SRR003",
    flag_reason = "barcode_outlier", n_study_samples = 10, n_study_barcode_true = 9,
    trim_mean_length = 33, study_median_len = 35
  )
  outliers_path <- tempfile(fileext = ".csv")
  data.table::fwrite(outliers, outliers_path)

  testthat::local_mocked_bindings(sample_currently_processing = function(exp) TRUE)
  testthat::local_mocked_bindings(last_real_touch = function(exp, config) as.POSIXct("2020-01-01", tz = "UTC"))

  cands <- barcode_fix_candidates(config, outliers_path)
  expect_identical(nrow(cands), 0L)
})

test_that("barcode_fix_candidates() respects min_majority_frac and ranks by it descending", {
  config <- fake_config()
  outliers <- data.table::data.table(
    study_accession = c("PRJNA004", "PRJNA005"), ScientificName = "Homo sapiens",
    raw_library = c("SRR004", "SRR005"), flag_reason = "barcode_outlier",
    n_study_samples = c(10, 10), n_study_barcode_true = c(6, 9), # 0.6 and 0.9
    trim_mean_length = 33, study_median_len = 35
  )
  outliers_path <- tempfile(fileext = ".csv")
  data.table::fwrite(outliers, outliers_path)

  testthat::local_mocked_bindings(sample_currently_processing = function(exp) FALSE)
  testthat::local_mocked_bindings(last_real_touch = function(exp, config) as.POSIXct("2020-01-01", tz = "UTC"))

  cands <- barcode_fix_candidates(config, outliers_path, min_majority_frac = 0.7)
  expect_identical(cands$run, "SRR005")
})

test_that("redetect_barcode_for_sample() calls barcode_detector_single() with check_at_mean_size forced to 0", {
  config <- fake_config()
  dir.create(dirname(config$complete_metadata), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE"), config$complete_metadata)
  stub <- fake_experiment_stub(run_ids = "SRR001", exp_name = "PRJNA001-homo_sapiens")
  bam_dir <- tempfile("bam_")
  captured <- NULL
  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir,
    barcode_detector_single = function(study_sample, fastq_dir, bam_root, trimmed_dir, ...) {
      captured <<- list(...)
      data.table::data.table(id = study_sample$Run)
    }
  )

  res <- redetect_barcode_for_sample("PRJNA001-homo_sapiens", "SRR001", config)

  expect_identical(res$id, "SRR001")
  expect_equal(captured$check_at_mean_size, 0)
  expect_true(captured$redownload_raw_if_needed)
})

test_that("redetect_barcode_for_sample() errors clearly when the run isn't in complete_metadata", {
  config <- fake_config()
  dir.create(dirname(config$complete_metadata), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(data.table::data.table(Run = "SRR999"), config$complete_metadata)
  expect_error(redetect_barcode_for_sample("PRJNA001-homo_sapiens", "SRR001", config), "not found")
})

test_that("apply_trim_fix_and_rerun() refuses to proceed when the study is currently being processed", {
  testthat::local_mocked_bindings(sample_currently_processing = function(exp) TRUE)
  expect_error(
    apply_trim_fix_and_rerun("PRJNA001-homo_sapiens", "SRR001", fake_config(), barcode5p_size = 5, barcode3p_size = 0),
    "currently being processed"
  )
})

test_that("apply_trim_fix_and_rerun() reports success and writes a status file when nothing unexpected changes", {
  config <- fake_config(preset = "Ribo-seq")
  exp <- "PRJNA001-homo_sapiens"
  dir.create(dirname(config$complete_metadata), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(data.table::data.table(study_accession = "PRJNA001", Run = c("SRR001", "SRR002")),
                     config$complete_metadata)
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"), exp_name = exp)
  bam_dir <- tempfile("bam_")
  dir.create(bam_dir, recursive = TRUE)
  writeLines("x", file.path(bam_dir, "SRR002.bam")) # sibling, must stay untouched

  testthat::local_mocked_bindings(sample_currently_processing = function(exp) FALSE)
  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir
  )
  testthat::local_mocked_bindings(apply_trim_fix_to_sample = function(...) invisible(NULL))
  testthat::local_mocked_bindings(pipeline_init_all = function(...) list(PRJNA001 = list()))
  testthat::local_mocked_bindings(run_pipeline = function(...) invisible(NULL))
  fake_mark_all_done(config, names(config$flag), exp)

  result <- apply_trim_fix_and_rerun(exp, "SRR001", config, barcode5p_size = 5, barcode3p_size = 0, note = "test")

  expect_true(result$success)
  expect_length(result$touched_sibling_files, 0)
  status_path <- file.path(config$project, "ai_fix_status", paste0(exp, ".rds"))
  expect_true(file.exists(status_path))
  expect_identical(readRDS(status_path)$exp, exp)
})

test_that("apply_trim_fix_and_rerun() flags failure when a sibling's own file was unexpectedly touched", {
  config <- fake_config(preset = "Ribo-seq")
  exp <- "PRJNA001-homo_sapiens"
  dir.create(dirname(config$complete_metadata), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(data.table::data.table(study_accession = "PRJNA001", Run = c("SRR001", "SRR002")),
                     config$complete_metadata)
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"), exp_name = exp)
  bam_dir <- tempfile("bam_")
  dir.create(bam_dir, recursive = TRUE)
  sibling_file <- file.path(bam_dir, "SRR002.bam")
  writeLines("x", sibling_file)

  testthat::local_mocked_bindings(sample_currently_processing = function(exp) FALSE)
  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir
  )
  testthat::local_mocked_bindings(apply_trim_fix_to_sample = function(...) invisible(NULL))
  # Simulate a bug: run_pipeline() incorrectly touches the sibling's file.
  testthat::local_mocked_bindings(pipeline_init_all = function(...) list(PRJNA001 = list()))
  testthat::local_mocked_bindings(run_pipeline = function(...) {
    Sys.sleep(1.1)
    writeLines("mutated", sibling_file)
  })
  fake_mark_all_done(config, names(config$flag), exp)

  result <- apply_trim_fix_and_rerun(exp, "SRR001", config, barcode5p_size = 5, barcode3p_size = 0)

  expect_false(result$success)
  expect_identical(result$touched_sibling_files, sibling_file)
})
