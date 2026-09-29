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
