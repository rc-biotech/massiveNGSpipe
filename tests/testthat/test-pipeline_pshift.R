# Unit tests for pipeline_pshift()'s new read-length-distribution hook
# (R/pipeline_preset_steps_sub.R). shiftFootprintsByExperimentSafe()
# (real ORFik shifting) and save_pshifted_length_distributions() (real
# ofst I/O) are both mocked -- this verifies only massiveNGSpipe's own
# new orchestration decision: the distribution hook fires on a
# successful shift and not on a failed one.

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
    shiftFootprintsByExperimentSafe = function(df, ...) structure(list(), class = "error"),
    save_pshifted_length_distributions = function(df) called <<- TRUE
  )

  pipeline_pshift(list(stub), config)

  expect_false(called)
  expect_false(step_is_done(config, "pshifted", name(stub)))
})
