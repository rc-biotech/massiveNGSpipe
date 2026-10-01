test_that("qc_verdict_report() returns zero rows (not an error) when nothing has been recorded", {
  config <- fake_config()
  pipelines <- fake_pipelines()

  result <- qc_verdict_report(pipelines, config)
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 0)
})

test_that("qc_verdict_report() reports one row per experiment with a recorded diagnostics file", {
  config <- fake_config()
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(bam_dir = bam_dir)

  update_qc_diagnostics(qc_diagnostics_path(bam_dir), too_few_reads = FALSE,
                        raw_reads = 5e6, wrong_organism = TRUE, mapping_rate_pct = 4.2,
                        periodicity_status = "warning", primary_cause = "wrong_organism")

  result <- qc_verdict_report(pipelines, config)
  expect_equal(nrow(result), 1)
  expect_identical(result$experiment, "PRJNA000001-homo_sapiens")
  expect_false(result$too_few_reads)
  expect_true(result$wrong_organism)
  expect_identical(result$periodicity_status, "warning")
  expect_identical(result$primary_cause, "wrong_organism")
})

test_that("qc_verdict_report() fills unrecorded fields with NA, not an error", {
  config <- fake_config()
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(bam_dir = bam_dir)
  # Only the trim-stage signal has been recorded so far.
  update_qc_diagnostics(qc_diagnostics_path(bam_dir), too_few_reads = FALSE, raw_reads = 5e6)

  result <- qc_verdict_report(pipelines, config)
  expect_equal(nrow(result), 1)
  expect_false(result$too_few_reads)
  expect_identical(result$wrong_organism, NA)
  expect_identical(result$periodicity_status, NA_character_)
  expect_identical(result$primary_cause, NA_character_)
})

test_that("qc_verdict_report(only_flagged = TRUE) drops clean studies", {
  config <- fake_config()
  bam_clean <- tempfile("bam_clean_")
  bam_bad <- tempfile("bam_bad_")
  pipelines <- c(
    fake_pipelines(accession = "PRJNA000001", bam_dir = bam_clean),
    fake_pipelines(accession = "PRJNA000002", bam_dir = bam_bad)
  )
  update_qc_diagnostics(qc_diagnostics_path(bam_clean), too_few_reads = FALSE,
                        periodicity_status = "good", primary_cause = "unknown")
  update_qc_diagnostics(qc_diagnostics_path(bam_bad), too_few_reads = TRUE,
                        periodicity_status = "warning", primary_cause = "too_few_reads")

  result <- qc_verdict_report(pipelines, config, only_flagged = TRUE)
  expect_equal(nrow(result), 1)
  expect_identical(result$experiment, "PRJNA000002-homo_sapiens")
})

test_that("qc_verdict_report(only_flagged = TRUE) keeps a study flagged only via periodicity_status", {
  config <- fake_config()
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(bam_dir = bam_dir)
  # no_data with no classified primary_cause -- must still count as flagged.
  update_qc_diagnostics(qc_diagnostics_path(bam_dir), periodicity_status = "no_data")

  result <- qc_verdict_report(pipelines, config, only_flagged = TRUE)
  expect_equal(nrow(result), 1)
})

test_that("pshift_triage() assembles qc/frames/length-distribution data purely from precomputed files", {
  config <- fake_config()
  bam_dir <- tempfile("bam_")
  qc_dir <- file.path(bam_dir, "aligned", "QC_STATS")
  # read_dist_or_null() derives its path as
  # dirname(filepath(df, type))/read_length_distribution/00_aggregated.csv
  # -- so the parent of each *_parent dir below is what the mocked
  # filepath() must return, not the distribution folder itself.
  ofst_parent <- tempfile("ofst_parent_")
  pshifted_parent <- tempfile("pshifted_parent_")
  dir.create(file.path(ofst_parent, "read_length_distribution"), recursive = TRUE)
  dir.create(file.path(pshifted_parent, "read_length_distribution"), recursive = TRUE)

  update_qc_diagnostics(qc_diagnostics_path(bam_dir), primary_cause = "wrong_organism")
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(data.table::data.table(frame = 0, fraction = "RFP", length = 28, score = 10),
                     file.path(qc_dir, "Ribo_frames_all.csv"))
  data.table::fwrite(data.table::data.table(read_length = 28L, count = 100, percent = 100),
                     file.path(ofst_parent, "read_length_distribution", "00_aggregated.csv"))
  data.table::fwrite(data.table::data.table(read_length = 28L, count = 80, percent = 100),
                     file.path(pshifted_parent, "read_length_distribution", "00_aggregated.csv"))

  stub <- fake_experiment_stub(exp_name = "PRJNA000001-homo_sapiens")
  testthat::local_mocked_bindings(
    `read.experiment` = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir,
    QCfolder = function(x) qc_dir,
    filepath = function(df, type, ...) {
      if (type == "ofst") file.path(ofst_parent, "SRR001.ofst")
      else file.path(pshifted_parent, "SRR001_pshifted.ofst")
    }
  )

  result <- pshift_triage("PRJNA000001-homo_sapiens", config)

  expect_identical(result$exp, "PRJNA000001-homo_sapiens")
  expect_identical(result$qc$primary_cause, "wrong_organism")
  expect_equal(result$frames$score, 10)
  expect_equal(result$length_dist_preshift$count, 100)
  expect_equal(result$length_dist_postshift$count, 80)
  expect_identical(result$accepted_lengths, config$accepted_lengths_rpf)
})

test_that("pshift_triage() returns NULL (not an error) for pieces that don't exist yet", {
  config <- fake_config()
  bam_dir <- tempfile("bam_") # nothing written under it at all
  stub <- fake_experiment_stub(exp_name = "PRJNA000002-homo_sapiens")
  testthat::local_mocked_bindings(
    `read.experiment` = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir,
    QCfolder = function(x) file.path(bam_dir, "aligned", "QC_STATS"),
    filepath = function(df, type, ...) file.path(tempfile(), "missing.ofst")
  )

  result <- pshift_triage("PRJNA000002-homo_sapiens", config)

  expect_identical(result$qc, list())
  expect_null(result$frames)
  expect_null(result$length_dist_preshift)
  expect_null(result$length_dist_postshift)
})
