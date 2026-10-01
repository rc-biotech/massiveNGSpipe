# Unit tests for manual_reshift()/log_manual_reshift() (R/manual_reshift.R).
# Real ORFik calls (shiftFootprintsByExperiment, shifts_save, shiftPlots,
# read.experiment) and the heavier massiveNGSpipe helpers (shift_qc(),
# save_pshifted_length_distributions()) are all mocked -- this verifies
# only massiveNGSpipe's own orchestration: the right calls happen in the
# right order, the manually_checked_shifts.rds marker and audit log get
# written, and flags are cleared/set correctly.

test_that("manual_reshift() writes the marker, clears downstream flags, sets pshifted, and logs", {
  config <- fake_config(preset = "Ribo-seq")
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"), exp_name = "PRJNA000001-homo_sapiens")
  qc_dir <- tempfile("qc_")
  offsets <- data.table::data.table(fraction = 28:30, offsets_start = -12)

  shift_calls <- list()
  testthat::local_mocked_bindings(
    read.experiment = function(exp, ...) stub,
    filepath = function(df, type, ...) paste0("/fake/", type, "/", runIDs(df), ".ofst"),
    libFolder = function(x, ...) "/fake/aligned",
    QCfolder = function(x) qc_dir
  )
  testthat::local_mocked_bindings(
    shiftFootprintsByExperiment = function(df, output_format, shift.list) {
      shift_calls[["shift"]] <<- shift.list
    },
    shifts_save = function(shifts, folder) shift_calls[["save"]] <<- list(shifts = shifts, folder = folder),
    shiftPlots = function(df, ...) shift_calls[["plots"]] <<- TRUE,
    .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    shift_qc = function(df, ...) shift_calls[["qc"]] <<- TRUE,
    save_pshifted_length_distributions = function(df) shift_calls[["dist"]] <<- TRUE
  )

  manual_reshift("PRJNA000001-homo_sapiens", offsets, config, note = "test fix")

  # The real shift + its persistence + QC regeneration all happened.
  expect_length(shift_calls[["shift"]], 2) # one entry per ofst file
  expect_identical(shift_calls[["save"]]$shifts, shift_calls[["shift"]])
  expect_true(shift_calls[["plots"]])
  expect_true(shift_calls[["qc"]])
  expect_true(shift_calls[["dist"]])

  # Marker file written where QCfolder() points.
  expect_true(file.exists(file.path(qc_dir, "manually_checked_shifts.rds")))
  expect_true(readRDS(file.path(qc_dir, "manually_checked_shifts.rds")))

  # pshifted marked done; everything downstream of it cleared (not done).
  expect_true(step_is_done(config, "pshifted", "PRJNA000001-homo_sapiens"))
  expect_false(step_is_done(config, "valid_pshift", "PRJNA000001-homo_sapiens"))
  expect_false(step_is_done(config, "merged_lib", "PRJNA000001-homo_sapiens"))

  # Audit log row written.
  log_path <- file.path(config$project, "manual_reshift_log.csv")
  expect_true(file.exists(log_path))
  log <- data.table::fread(log_path)
  expect_identical(log$exp, "PRJNA000001-homo_sapiens")
  expect_identical(log$note, "test fix")
  expect_match(log$offsets, "28:-12;29:-12;30:-12")
})

test_that("manual_reshift() rejects a malformed offsets table before touching anything", {
  config <- fake_config(preset = "Ribo-seq")
  bad_offsets <- data.table::data.table(wrong_col = 1:3)
  expect_error(manual_reshift("PRJNA000001-homo_sapiens", bad_offsets, config), "fraction")
})

test_that("log_manual_reshift() appends, not overwrites, across multiple calls", {
  config <- fake_config()
  offsets <- data.table::data.table(fraction = 28, offsets_start = -12)

  log_manual_reshift("exp1", offsets, "first", config)
  log_manual_reshift("exp2", offsets, "second", config)

  log <- data.table::fread(file.path(config$project, "manual_reshift_log.csv"))
  expect_identical(log$exp, c("exp1", "exp2"))
  expect_identical(log$note, c("first", "second"))
})
