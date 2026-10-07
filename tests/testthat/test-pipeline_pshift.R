# Unit tests for pipeline_pshift()'s new read-length-distribution hook
# AND its per-sample resume subsetting (both R/pipeline_preset_steps_sub.R).
# shiftFootprintsByExperimentSafe() (real ORFik shifting) and
# save_pshifted_length_distributions() (real ofst I/O) are both mocked --
# this verifies only massiveNGSpipe's own orchestration: the distribution
# hook fires on a successful shift and not on a failed one, and only
# not-yet-valid samples get passed through to the (mocked) real shift call.
#
# pshift_sample_valid()/pshift_needed_subset() resolve ofst/pshifted
# paths via the ofst_filepath()/pshifted_filepath() wrapper functions
# (not ORFik::filepath() directly) specifically so these can be mocked
# against a fake_exp_stub -- filepath() itself is a plain (non-S4-generic)
# ORFik function with an internal stopifnot(is(df, "experiment")) that
# can never dispatch on a test double.

test_that("pipeline_pshift() saves length distributions after a successful all-mappers shift", {
  config <- fake_config(preset = "Ribo-seq", extra = list(all_mappers = TRUE, split_unique_mappers = FALSE,
                                     reuse_shifts_if_existing = FALSE,
                                     accepted_lengths_rpf = c(20, 21, 25:33),
                                     max_no_adapter_removed_pct = 80))
  stub <- fake_experiment_stub(run_ids = c("SRR001"))
  fake_mark_all_done(config, c("trim", "collapsed", "aligned", "cleanbam", "exp", "ofst"),
                     name(stub))

  called_with <- NULL
  testthat::local_mocked_bindings(
    # No ofst file "exists" for this fake path -- pshift_sample_valid()
    # returns FALSE, so this sample is treated as needing a shift.
    ofst_filepath = function(df) file.path(tempdir(), "does_not_exist.ofst"),
    shiftFootprintsByExperimentSafe = function(df, ...) "ok", # not an "error"-classed object
    save_pshifted_length_distributions = function(df) called_with <<- df
  )

  pipeline_pshift(list(stub), config)

  expect_false(is.null(called_with))
  expect_true(step_is_done(config, "pshifted", name(stub)))
})

test_that("pipeline_pshift() does NOT save length distributions after a failed shift", {
  config <- fake_config(preset = "Ribo-seq", extra = list(all_mappers = TRUE, split_unique_mappers = FALSE,
                                     reuse_shifts_if_existing = FALSE,
                                     accepted_lengths_rpf = c(20, 21, 25:33),
                                     max_no_adapter_removed_pct = 80))
  stub <- fake_experiment_stub(run_ids = c("SRR001"))
  fake_mark_all_done(config, c("trim", "collapsed", "aligned", "cleanbam", "exp", "ofst"),
                     name(stub))

  called <- FALSE
  testthat::local_mocked_bindings(
    ofst_filepath = function(df) file.path(tempdir(), "does_not_exist.ofst"),
    shiftFootprintsByExperimentSafe = function(df, ...) structure(list(), class = "error"),
    save_pshifted_length_distributions = function(df) called <<- TRUE
  )

  pipeline_pshift(list(stub), config)

  expect_false(called)
  expect_false(step_is_done(config, "pshifted", name(stub)))
})

test_that("pipeline_pshift() skips a sample whose pshifted output is already valid, passing only the subset still needing it to the real shift call", {
  config <- fake_config(preset = "Ribo-seq", extra = list(all_mappers = TRUE, split_unique_mappers = FALSE,
                                     reuse_shifts_if_existing = FALSE,
                                     accepted_lengths_rpf = c(20, 21, 25:33),
                                     max_no_adapter_removed_pct = 80))
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"))
  fake_mark_all_done(config, c("trim", "collapsed", "aligned", "cleanbam", "exp", "ofst"),
                     name(stub))

  # SRR001: valid (ofst older than pshifted). SRR002: needs a shift (no
  # pshifted file at all). ofst_filepath()/pshifted_filepath() are keyed
  # off run id so each row gets its own deterministic fake path.
  ofst_dir <- tempfile("ofst_"); dir.create(ofst_dir)
  pshifted_dir <- tempfile("pshifted_"); dir.create(pshifted_dir)
  ofst_paths <- setNames(file.path(ofst_dir, paste0(c("SRR001", "SRR002"), ".ofst")), c("SRR001", "SRR002"))
  for (p in ofst_paths) file.create(p)
  pshifted_srr001 <- file.path(pshifted_dir, "SRR001_pshifted.ofst")
  file.create(pshifted_srr001)
  Sys.setFileTime(pshifted_srr001, Sys.time() + 10) # newer than its ofst -> valid

  subset_runs_seen <- NULL
  testthat::local_mocked_bindings(
    ofst_filepath = function(df) unname(ofst_paths[runIDs(df)]),
    pshifted_filepath = function(df_one_row) {
      if (runIDs(df_one_row) == "SRR001") pshifted_srr001 else file.path(pshifted_dir, "does_not_exist.ofst")
    },
    shiftFootprintsByExperimentSafe = function(df, ...) { subset_runs_seen <<- runIDs(df); "ok" },
    save_pshifted_length_distributions = function(df) NULL
  )

  pipeline_pshift(list(stub), config)

  expect_identical(subset_runs_seen, "SRR002")
})
