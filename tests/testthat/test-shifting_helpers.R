test_that("periodicity_check_flag() returns good and writes good.rds when periodicity is fine", {
  qc_dir <- tempfile("qc_")
  dir.create(qc_dir, recursive = TRUE)
  zero_frame <- data.table::data.table(percent_length = c(40, 60))

  status <- periodicity_check_flag(zero_frame, qc_dir, has_cds_data = TRUE)

  expect_identical(status, "good")
  expect_true(file.exists(file.path(qc_dir, "good.rds")))
  expect_false(file.exists(file.path(qc_dir, "warning.rds")))
  expect_false(file.exists(file.path(qc_dir, "no_data.rds")))
})

test_that("periodicity_check_flag() returns warning when a read length is under 25% CDS coverage", {
  qc_dir <- tempfile("qc_")
  dir.create(qc_dir, recursive = TRUE)
  zero_frame <- data.table::data.table(percent_length = c(10, 60))

  expect_warning(periodicity_check_flag(zero_frame, qc_dir, has_cds_data = TRUE))
  status <- suppressWarnings(periodicity_check_flag(zero_frame, qc_dir, has_cds_data = TRUE))

  expect_identical(status, "warning")
  expect_true(file.exists(file.path(qc_dir, "warning.rds")))
  expect_false(file.exists(file.path(qc_dir, "good.rds")))
  expect_false(file.exists(file.path(qc_dir, "no_data.rds")))
})

test_that("periodicity_check_flag() returns no_data (not good) when there is zero CDS frame data", {
  qc_dir <- tempfile("qc_")
  dir.create(qc_dir, recursive = TRUE)
  # Empty zero_frame with has_cds_data = FALSE is the "nothing aligned to
  # real ORFs" case -- this must NOT be reported as a clean "good" pass,
  # which is exactly the bug this fixes (nrow(zero_frame) == 0 used to
  # mean both "perfect periodicity" and "no data" identically).
  zero_frame <- data.table::data.table(percent_length = numeric())

  expect_warning(periodicity_check_flag(zero_frame, qc_dir, has_cds_data = FALSE))
  status <- suppressWarnings(periodicity_check_flag(zero_frame, qc_dir, has_cds_data = FALSE))

  expect_identical(status, "no_data")
  expect_true(file.exists(file.path(qc_dir, "no_data.rds")))
  expect_false(file.exists(file.path(qc_dir, "good.rds")))
  expect_false(file.exists(file.path(qc_dir, "warning.rds")))
})

test_that("periodicity_check_flag() with an empty zero_frame and real CDS data is a genuine good", {
  qc_dir <- tempfile("qc_")
  dir.create(qc_dir, recursive = TRUE)
  # Empty zero_frame with has_cds_data = TRUE means every read length
  # cleared the 25% bar -- a real "good", distinct from the no_data case.
  zero_frame <- data.table::data.table(percent_length = numeric())

  status <- periodicity_check_flag(zero_frame, qc_dir, has_cds_data = TRUE)

  expect_identical(status, "good")
  expect_true(file.exists(file.path(qc_dir, "good.rds")))
})

test_that("periodicity_check_flag()'s default has_cds_data reproduces the pre-fix behavior", {
  qc_dir <- tempfile("qc_")
  dir.create(qc_dir, recursive = TRUE)
  zero_frame <- data.table::data.table(percent_length = numeric())

  # No has_cds_data passed: falls back to nrow(zero_frame) > 0, matching
  # the old (buggy) boolean, for any caller other than shift_qc() itself.
  expect_warning(periodicity_check_flag(zero_frame, qc_dir))
  status <- suppressWarnings(periodicity_check_flag(zero_frame, qc_dir))

  expect_identical(status, "no_data")
})

test_that("periodicity_check_flag() re-run switches status and cleans up the old flag file", {
  qc_dir <- tempfile("qc_")
  dir.create(qc_dir, recursive = TRUE)

  periodicity_check_flag(data.table::data.table(percent_length = c(90)), qc_dir, has_cds_data = TRUE)
  expect_true(file.exists(file.path(qc_dir, "good.rds")))

  expect_warning(periodicity_check_flag(data.table::data.table(percent_length = c(5)), qc_dir, has_cds_data = TRUE))
  expect_true(file.exists(file.path(qc_dir, "warning.rds")))
  expect_false(file.exists(file.path(qc_dir, "good.rds")))
})

test_that("periodicity_check_flag() records periodicity_status into qc_diagnostics.rds", {
  qc_dir <- tempfile("qc_")
  dir.create(qc_dir, recursive = TRUE)
  periodicity_check_flag(data.table::data.table(percent_length = c(90)), qc_dir, has_cds_data = TRUE)

  diag <- read_qc_diagnostics(file.path(qc_dir, "qc_diagnostics.rds"))
  expect_identical(diag$periodicity_status, "good")
})

# regionPerReadLengthPerLib() itself (R/shifting_helpers.R) has the same
# fill=TRUE fix applied as shift_qc_cache.R's own aggregation (2026-10-08,
# same root cause: a library with zero rows never goes through the
# `if (nrow(total) > 0)` column-augmentation, so rbindlist() without
# fill=TRUE would error on the resulting column-count mismatch). No
# dedicated test here: its bare (Depends-resolved, not `BiocParallel::`-
# qualified) `bplapply()` call isn't interceptable via
# testthat::local_mocked_bindings(.package = "BiocParallel") the way an
# explicitly-qualified call is (confirmed by direct attempt -- the real
# dispatch ran regardless of the mock), and a real end-to-end test would
# need ORFik::regionPerReadLength() to actually run against real
# GAlignments/GRanges input. The fix itself is still correct by direct
# inspection -- same one-line change, same proven mechanism.
