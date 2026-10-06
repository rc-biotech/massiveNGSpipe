# Unit tests for R/manual_trim_fix.R -- the manual adapter/barcode
# override writer and the per-sample flag-clearing orchestration.
# Companion to test-manual_reshift.R (same style/scope: massiveNGSpipe's
# own orchestration, not real trim/collapse/align execution).

test_that("upsert_manual_override_row writes a fresh file when none exists", {
  dir <- tempfile(); dir.create(dir)
  path <- file.path(dir, "barcodes_manual.csv")
  upsert_manual_override_row(path, data.table::data.table(Run = "SRR001", barcode5p_size = 5, barcode3p_size = 2))
  out <- data.table::fread(path)
  expect_identical(out$Run, "SRR001")
  expect_equal(out$barcode5p_size, 5)
})

test_that("upsert_manual_override_row replaces only the matching Run's row, preserving the rest", {
  path <- tempfile(fileext = ".csv")
  data.table::fwrite(data.table::data.table(Run = c("SRR001", "SRR002"),
                                            barcode5p_size = c(1, 1), barcode3p_size = c(1, 1)), path)
  upsert_manual_override_row(path, data.table::data.table(Run = "SRR002", barcode5p_size = 5, barcode3p_size = 2))
  out <- data.table::fread(path)[order(Run)]
  expect_identical(out$Run, c("SRR001", "SRR002"))
  expect_equal(out[Run == "SRR001"]$barcode5p_size, 1)
  expect_equal(out[Run == "SRR002"]$barcode5p_size, 5)
})

test_that("log_trim_fix appends, not overwrites, across multiple calls", {
  config <- fake_config()
  log_trim_fix("exp1", "SRR001", 5, 2, NULL, "first", config)
  log_trim_fix("exp2", "SRR002", NA, NA, "AGATCGGAAGAG", "second", config)
  log <- data.table::fread(file.path(config$project, "manual_trim_fix_log.csv"))
  expect_identical(log$exp, c("exp1", "exp2"))
  expect_identical(log$note, c("first", "second"))
})

test_that("apply_trim_fix_to_sample() writes the manual override, clears this sample's per-sample markers, clears downstream flags, and logs", {
  config <- fake_config(preset = "Ribo-seq")
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"), exp_name = "PRJNA000001-homo_sapiens")
  bam_dir <- tempfile("bam_")

  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir
  )

  # Pre-existing per-sample markers for both samples, across every
  # per-sample-tracked step -- only SRR001's should be removed.
  for (step_id in c("trim", "collapsed", "aligned", "ofst", "covrle", "bigwig")) {
    set_sample_flag(config, step_id, "PRJNA000001-homo_sapiens", "SRR001")
    set_sample_flag(config, step_id, "PRJNA000001-homo_sapiens", "SRR002")
  }
  fake_mark_all_done(config, names(config$flag), "PRJNA000001-homo_sapiens")

  apply_trim_fix_to_sample("PRJNA000001-homo_sapiens", "SRR001", config,
                           barcode5p_size = 5, barcode3p_size = 2, note = "undetected barcode")

  # Manual override written.
  barcodes <- data.table::fread(file.path(bam_dir, "trim", "barcodes_manual.csv"))
  expect_identical(barcodes$Run, "SRR001")
  expect_equal(barcodes$barcode5p_size, 5)

  # SRR001's per-sample markers gone, SRR002's untouched.
  for (step_id in c("trim", "collapsed", "aligned", "ofst", "covrle", "bigwig")) {
    expect_identical(samples_done(config, step_id, "PRJNA000001-homo_sapiens"), "SRR002")
  }

  # Experiment-level flags from "trim" onward cleared; "fetch" (before
  # trim) untouched.
  downstream <- names(config$flag)
  downstream <- downstream[which(downstream == "trim"):length(downstream)]
  for (s in downstream) expect_false(step_is_done(config, s, "PRJNA000001-homo_sapiens"))
  if ("fetch" %in% names(config$flag)) expect_true(step_is_done(config, "fetch", "PRJNA000001-homo_sapiens"))

  # Audit log written.
  log <- data.table::fread(file.path(config$project, "manual_trim_fix_log.csv"))
  expect_identical(log$run, "SRR001")
  expect_identical(log$note, "undetected barcode")
})

test_that("apply_trim_fix_to_sample() also writes adapters_manual.csv when an adapter correction is given", {
  config <- fake_config(preset = "Ribo-seq")
  stub <- fake_experiment_stub(run_ids = "SRR001", exp_name = "PRJNA000001-homo_sapiens")
  bam_dir <- tempfile("bam_")
  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir
  )
  fake_mark_all_done(config, names(config$flag), "PRJNA000001-homo_sapiens")

  apply_trim_fix_to_sample("PRJNA000001-homo_sapiens", "SRR001", config, adapter = "AGATCGGAAGAG")

  adapters <- data.table::fread(file.path(bam_dir, "trim", "adapters_manual.csv"))
  expect_identical(adapters$Run, "SRR001")
  expect_identical(adapters$adapter, "AGATCGGAAGAG")
  expect_false(file.exists(file.path(bam_dir, "trim", "barcodes_manual.csv")))
})

test_that("apply_trim_fix_to_sample() requires at least one correction and both barcode sizes together", {
  config <- fake_config()
  expect_error(apply_trim_fix_to_sample("exp1", "SRR001", config))
  expect_error(apply_trim_fix_to_sample("exp1", "SRR001", config, barcode5p_size = 5))
})

test_that("apply_trim_fix_to_sample() backfills sibling markers for a study with NO pre-existing per-sample markers at all", {
  # Reproduces a study fully processed before the per-sample-resume
  # feature existed: experiment-level flags are TRUE, but samples_done()
  # is empty for every sample, not just the one being fixed. Without the
  # backfill, pipeline_trim()/pipeline_collapse() would treat every
  # sibling as not-done too, not just the fixed sample -- confirmed live
  # on PRJNA926112-homo_sapiens, 2026-10-05.
  config <- fake_config(preset = "Ribo-seq")
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002", "SRR003"), exp_name = "PRJNA000001-homo_sapiens")
  bam_dir <- tempfile("bam_")
  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir
  )
  # No set_sample_flag() calls at all -- only the experiment-level flags,
  # matching a pre-feature study's actual on-disk state.
  fake_mark_all_done(config, names(config$flag), "PRJNA000001-homo_sapiens")

  apply_trim_fix_to_sample("PRJNA000001-homo_sapiens", "SRR001", config,
                           barcode5p_size = 7, barcode3p_size = 0, note = "legacy study backfill")

  # Siblings backfilled to "done"; the fixed sample has no marker.
  for (step_id in c("trim", "collapsed", "aligned", "ofst", "covrle", "bigwig")) {
    expect_setequal(samples_done(config, step_id, "PRJNA000001-homo_sapiens"), c("SRR002", "SRR003"))
  }
})

test_that("apply_trim_fix_to_sample()'s backfilled 'trim' markers are rbindlist()-compatible, reusing the existing adapter_barcode_table.csv row when present", {
  # pipeline_trim() reconstructs adapter_barcode_table.csv via
  # rbindlist(sample_flag_values(config, "trim", experiment), fill = TRUE)
  # -- a plain TRUE backfill value breaks that with "Item 1 of input is
  # not a data.frame, data.table or list" (confirmed live, 2026-10-05,
  # PRJNA926112-homo_sapiens, right after the real fastp trim for the
  # fixed sample itself had already succeeded).
  config <- fake_config(preset = "Ribo-seq")
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"), exp_name = "PRJNA000001-homo_sapiens")
  bam_dir <- tempfile("bam_")
  dir.create(file.path(bam_dir, "trim"), recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(id = c("SRR001", "SRR002"), adapter = "passed",
                           barcode_detected = c(FALSE, TRUE), barcode5p_size = c(0, 7)),
    file.path(bam_dir, "trim", "adapter_barcode_table.csv"))
  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    bam_dir_from_df = function(df) bam_dir
  )
  fake_mark_all_done(config, names(config$flag), "PRJNA000001-homo_sapiens")

  apply_trim_fix_to_sample("PRJNA000001-homo_sapiens", "SRR001", config,
                           barcode5p_size = 5, barcode3p_size = 2, note = "rbindlist compatibility")

  rebuilt <- data.table::rbindlist(sample_flag_values(config, "trim", "PRJNA000001-homo_sapiens"), fill = TRUE)
  expect_identical(rebuilt$id, "SRR002")
  expect_equal(rebuilt$barcode5p_size, 7)
})
