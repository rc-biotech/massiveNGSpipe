# Note: qc_diagnostics_path()'s equivalence with the real ORFik
# QCfolder(df)-derived path (file.path(QCfolder(df), "qc_diagnostics.rds"))
# was verified live against a real processed experiment
# (PRJNA521496-homo_sapiens) rather than re-checked here: building a real
# ORFik experiment via create.experiment() is too heavy for this fast unit
# suite. What's checked here is qc_diagnostics_path()'s own path shape
# staying exactly what that live check confirmed matches QCfolder()'s
# layout (bam_dir/aligned/QC_STATS/), so a future edit that silently
# changes this shape gets caught.
test_that("qc_diagnostics_path() matches the verified QCfolder() layout", {
  expect_identical(qc_diagnostics_path("/some/bam/dir"),
                   "/some/bam/dir/aligned/QC_STATS/qc_diagnostics.rds")
})

test_that("update_qc_diagnostics()/read_qc_diagnostics() round-trip and merge fields", {
  path <- file.path(tempfile("qcdiag_"), "qc_diagnostics.rds")
  expect_identical(read_qc_diagnostics(path), list())

  update_qc_diagnostics(path, too_few_reads = FALSE, raw_reads = 5e6)
  rec1 <- read_qc_diagnostics(path)
  expect_identical(rec1$too_few_reads, FALSE)
  expect_identical(rec1$raw_reads, 5e6)
  expect_true(!is.null(rec1$last_updated))

  # A second update overwrites only the fields it names, leaving earlier
  # fields from a different stage untouched.
  update_qc_diagnostics(path, wrong_organism = TRUE, mapping_rate_pct = 12.3)
  rec2 <- read_qc_diagnostics(path)
  expect_identical(rec2$too_few_reads, FALSE)
  expect_identical(rec2$raw_reads, 5e6)
  expect_identical(rec2$wrong_organism, TRUE)
  expect_identical(rec2$mapping_rate_pct, 12.3)

  # Re-setting an already-set field overwrites in place, doesn't duplicate.
  update_qc_diagnostics(path, too_few_reads = TRUE)
  rec3 <- read_qc_diagnostics(path)
  expect_identical(rec3$too_few_reads, TRUE)
  expect_length(rec3[names(rec3) == "too_few_reads"], 1)
})

#' Write a minimal fastp-shaped JSON report, as ORFik::trimming.table()
#' actually reads it (see its source: json$summary$before_filtering /
#' after_filtering $total_reads, read1_mean_length).
write_fake_fastp_json <- function(dir, run, raw_reads, trim_reads,
                                  raw_mean_length = 30, trim_mean_length = 28) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  json <- list(summary = list(
    before_filtering = list(total_reads = raw_reads, read1_mean_length = raw_mean_length),
    after_filtering = list(total_reads = trim_reads, read1_mean_length = trim_mean_length)
  ))
  jsonlite::write_json(json, file.path(dir, paste0("report_", run, ".json")),
                       auto_unbox = TRUE)
}

test_that("check_too_few_reads() flags below-threshold raw read counts", {
  bam_dir <- tempfile("bam_")
  trimmed_dir <- tempfile("trim_")
  write_fake_fastp_json(trimmed_dir, "SRR001", raw_reads = 5e4, trim_reads = 4.8e4)

  check_too_few_reads(bam_dir, trimmed_dir, min_raw_reads = 1e5)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_true(diag$too_few_reads)
  # jsonlite parses whole-number totals as integer, not double -- compare
  # by value, not storage type.
  expect_equal(diag$raw_reads, 5e4)
})

test_that("check_too_few_reads() does not flag above-threshold raw read counts", {
  bam_dir <- tempfile("bam_")
  trimmed_dir <- tempfile("trim_")
  write_fake_fastp_json(trimmed_dir, "SRR001", raw_reads = 5e6, trim_reads = 4.8e6)

  check_too_few_reads(bam_dir, trimmed_dir, min_raw_reads = 1e5)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_false(diag$too_few_reads)
})

test_that("check_too_few_reads() sums across multiple runs in one study", {
  bam_dir <- tempfile("bam_")
  trimmed_dir <- tempfile("trim_")
  write_fake_fastp_json(trimmed_dir, "SRR001", raw_reads = 6e4, trim_reads = 5.8e4)
  write_fake_fastp_json(trimmed_dir, "SRR002", raw_reads = 6e4, trim_reads = 5.8e4)

  check_too_few_reads(bam_dir, trimmed_dir, min_raw_reads = 1e5)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_equal(diag$raw_reads, 1.2e5)
  expect_false(diag$too_few_reads)
})

test_that("check_too_few_reads() records NA (not FALSE) when no trim reports exist yet", {
  bam_dir <- tempfile("bam_")
  trimmed_dir <- tempfile("trim_empty_")
  dir.create(trimmed_dir, recursive = TRUE)

  check_too_few_reads(bam_dir, trimmed_dir, min_raw_reads = 1e5)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_identical(diag$too_few_reads, NA)
  expect_identical(diag$raw_reads, NA_real_)
})

test_that("full_process_csv_path() prefers SINGLE over plain, falls back, else NULL", {
  bam_dir <- tempfile("bam_")
  dir.create(bam_dir, recursive = TRUE)
  expect_null(full_process_csv_path(bam_dir))

  write.csv(data.frame(x = 1), file.path(bam_dir, "full_process.csv"), row.names = FALSE)
  expect_identical(full_process_csv_path(bam_dir), file.path(bam_dir, "full_process.csv"))

  write.csv(data.frame(x = 1), file.path(bam_dir, "full_process_SINGLE.csv"), row.names = FALSE)
  expect_identical(full_process_csv_path(bam_dir), file.path(bam_dir, "full_process_SINGLE.csv"))
})

test_that("check_alignment_rate() flags below-threshold unique-mapping rate", {
  bam_dir <- tempfile("bam_")
  dir.create(bam_dir, recursive = TRUE)
  data.table::fwrite(data.table::data.table(
    sample = "SRR001", `Uniquely mapped reads %` = "3.20%"
  ), file.path(bam_dir, "full_process_SINGLE.csv"))

  check_alignment_rate(bam_dir, min_alignment_rate = 10)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_true(diag$wrong_organism)
  expect_equal(diag$mapping_rate_pct, 3.2)
})

test_that("check_alignment_rate() does not flag above-threshold mapping rate", {
  bam_dir <- tempfile("bam_")
  dir.create(bam_dir, recursive = TRUE)
  data.table::fwrite(data.table::data.table(
    sample = "SRR001", `Uniquely mapped reads %` = "87.50%"
  ), file.path(bam_dir, "full_process_SINGLE.csv"))

  check_alignment_rate(bam_dir, min_alignment_rate = 10)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_false(diag$wrong_organism)
})

test_that("check_alignment_rate() records NA when full_process csv doesn't exist yet", {
  bam_dir <- tempfile("bam_")
  check_alignment_rate(bam_dir, min_alignment_rate = 10)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_identical(diag$wrong_organism, NA)
  expect_identical(diag$mapping_rate_pct, NA_real_)
})

test_that("check_adapter_barcode_quality() flags high no-adapter-removed fraction", {
  bam_dir <- tempfile("bam_")
  trimmed_dir <- tempfile("trim_")
  dir.create(trimmed_dir, recursive = TRUE)
  data.table::fwrite(data.table::data.table(
    sample = "SRR001", `reads_no_adapter_removed fastp(%)` = 92.4
  ), file.path(trimmed_dir, "adapter_barcode_table.csv"))

  check_adapter_barcode_quality(bam_dir, trimmed_dir, max_no_adapter_removed_pct = 80)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_true(diag$bad_adapter_barcode)
  expect_equal(diag$no_adapter_removed_pct, 92.4)
})

test_that("check_adapter_barcode_quality() records NA when table doesn't exist yet", {
  bam_dir <- tempfile("bam_")
  trimmed_dir <- tempfile("trim_missing_")
  check_adapter_barcode_quality(bam_dir, trimmed_dir, max_no_adapter_removed_pct = 80)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_identical(diag$bad_adapter_barcode, NA)
})

test_that("classify_bad_shift_cause() follows fixed precedence: too_few_reads first", {
  bam_dir <- tempfile("bam_")
  update_qc_diagnostics(qc_diagnostics_path(bam_dir), too_few_reads = TRUE,
                        wrong_organism = TRUE, bad_adapter_barcode = TRUE)
  expect_identical(classify_bad_shift_cause(bam_dir), "too_few_reads")
})

test_that("classify_bad_shift_cause() falls through to wrong_organism next", {
  bam_dir <- tempfile("bam_")
  update_qc_diagnostics(qc_diagnostics_path(bam_dir), too_few_reads = FALSE,
                        wrong_organism = TRUE, bad_adapter_barcode = TRUE)
  expect_identical(classify_bad_shift_cause(bam_dir), "wrong_organism")
})

test_that("classify_bad_shift_cause() falls through to adapter_barcode_umi last", {
  bam_dir <- tempfile("bam_")
  update_qc_diagnostics(qc_diagnostics_path(bam_dir), too_few_reads = FALSE,
                        wrong_organism = FALSE, bad_adapter_barcode = TRUE)
  expect_identical(classify_bad_shift_cause(bam_dir), "adapter_barcode_umi")
})

test_that("classify_bad_shift_cause() returns unknown when nothing is flagged or recorded", {
  bam_dir <- tempfile("bam_")
  expect_identical(classify_bad_shift_cause(bam_dir), "unknown")

  bam_dir2 <- tempfile("bam_")
  update_qc_diagnostics(qc_diagnostics_path(bam_dir2), too_few_reads = FALSE,
                        wrong_organism = FALSE, bad_adapter_barcode = FALSE)
  expect_identical(classify_bad_shift_cause(bam_dir2), "unknown")
})

test_that("classify_bad_shift_cause() treats NA (not-yet-checked) as not-this-cause", {
  bam_dir <- tempfile("bam_")
  update_qc_diagnostics(qc_diagnostics_path(bam_dir), too_few_reads = NA,
                        wrong_organism = NA, bad_adapter_barcode = TRUE)
  expect_identical(classify_bad_shift_cause(bam_dir), "adapter_barcode_umi")
})

test_that("classify_bad_shift_cause() writes primary_cause back into the record", {
  bam_dir <- tempfile("bam_")
  update_qc_diagnostics(qc_diagnostics_path(bam_dir), wrong_organism = TRUE)
  classify_bad_shift_cause(bam_dir)
  diag <- read_qc_diagnostics(qc_diagnostics_path(bam_dir))
  expect_identical(diag$primary_cause, "wrong_organism")
})
