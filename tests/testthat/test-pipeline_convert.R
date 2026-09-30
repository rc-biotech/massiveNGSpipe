# Unit tests for convert_per_sample() and the ofst/covRLE/bigwig pipeline
# wrappers in R/pipeline_preset_steps_sub.R -- per-sample resume support
# for stages that are plain serial per-file conversion loops. The real
# ORFik conversion functions (convert_bam_to_ofst/convert_to_covRleList/
# convert_to_bigWig) are out of scope, same as ORFik/STAR elsewhere in
# this suite: tests here use a mock convert_fun (for convert_per_sample())
# or mock convert_per_sample() itself (for the three wrappers), never real
# bam/genomic files.

test_that("convert_per_sample() calls convert_fun once per row and records a flag for each", {
  config <- fake_config()
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002", "SRR003"))
  seen <- character()
  convert_fun <- function(x) seen <<- c(seen, runIDs(x))

  convert_per_sample(stub, config, "ofst", convert_fun)

  expect_identical(seen, c("SRR001", "SRR002", "SRR003"))
  expect_identical(sort(samples_done(config, "ofst", name(stub))),
                   c("SRR001", "SRR002", "SRR003"))
})

test_that("convert_per_sample() skips samples already marked done (resume)", {
  config <- fake_config()
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002", "SRR003"))
  set_sample_flag(config, "ofst", name(stub), "SRR001")
  seen <- character()
  convert_fun <- function(x) seen <<- c(seen, runIDs(x))

  convert_per_sample(stub, config, "ofst", convert_fun)

  # Only the two not-yet-done samples get (re-)converted.
  expect_identical(seen, c("SRR002", "SRR003"))
  expect_identical(sort(samples_done(config, "ofst", name(stub))),
                   c("SRR001", "SRR002", "SRR003"))
})

test_that("convert_per_sample() does nothing when every sample is already done", {
  config <- fake_config()
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"))
  set_sample_flag(config, "ofst", name(stub), "SRR001")
  set_sample_flag(config, "ofst", name(stub), "SRR002")
  called <- FALSE
  convert_fun <- function(x) called <<- TRUE

  convert_per_sample(stub, config, "ofst", convert_fun)
  expect_false(called)
})

test_that("convert_per_sample() keeps all-mappers and _unique passes on independent resume points", {
  config <- fake_config()
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"))
  set_sample_flag(config, "ofst", name(stub), "SRR001")
  set_sample_flag(config, "ofst", name(stub), "SRR002")
  seen <- character()
  convert_fun <- function(x) seen <<- c(seen, runIDs(x))

  # "ofst" is fully done, but "ofst_unique" has never been touched --
  # must still run both samples through the _unique step id.
  convert_per_sample(stub, config, "ofst_unique", convert_fun)
  expect_identical(seen, c("SRR001", "SRR002"))
})

test_that("pipeline_create_ofst() drives convert_per_sample() for the all-mappers pass only when split_unique_mappers is FALSE", {
  config <- fake_config(extra = list(all_mappers = TRUE, split_unique_mappers = FALSE))
  stub <- fake_experiment_stub(run_ids = c("SRR001"))
  # step_is_next_not_done() needs every earlier step done, not just the
  # immediately-preceding one.
  fake_mark_all_done(config, c("aligned", "cleanbam", "exp"), name(stub))
  calls <- list()
  testthat::local_mocked_bindings(
    convert_per_sample = function(df, config, step_id, convert_fun) {
      calls[[length(calls) + 1]] <<- step_id
    }
  )

  pipeline_create_ofst(list(stub), config)

  expect_identical(unlist(calls), "ofst")
  expect_true(step_is_done(config, "ofst", name(stub)))
})

test_that("pipeline_create_ofst() also drives the _unique pass when split_unique_mappers is TRUE", {
  config <- fake_config(extra = list(all_mappers = TRUE, split_unique_mappers = TRUE))
  stub <- fake_experiment_stub(run_ids = c("SRR001"))
  fake_mark_all_done(config, c("aligned", "cleanbam", "exp"), name(stub))
  calls <- list()
  testthat::local_mocked_bindings(
    convert_per_sample = function(df, config, step_id, convert_fun) {
      calls[[length(calls) + 1]] <<- step_id
    }
  )

  pipeline_create_ofst(list(stub), config)

  expect_identical(unlist(calls), c("ofst", "ofst_unique"))
})

test_that("pipeline_create_ofst() skips an experiment whose ofst step is already done", {
  config <- fake_config(extra = list(all_mappers = TRUE, split_unique_mappers = FALSE))
  stub <- fake_experiment_stub(run_ids = c("SRR001"))
  fake_mark_all_done(config, c("aligned", "cleanbam", "exp", "ofst"), name(stub))
  called <- FALSE
  testthat::local_mocked_bindings(
    convert_per_sample = function(df, config, step_id, convert_fun) called <<- TRUE
  )

  pipeline_create_ofst(list(stub), config)
  expect_false(called)
})

test_that("pipeline_convert_covRLE() drives convert_per_sample() with covrle step ids", {
  config <- fake_config(extra = list(all_mappers = TRUE, split_unique_mappers = TRUE))
  stub <- fake_experiment_stub(run_ids = c("SRR001"))
  fake_mark_all_done(config, c("aligned", "cleanbam", "exp", "ofst", "cigar_collapse", "merged_lib"),
                     name(stub))
  calls <- list()
  testthat::local_mocked_bindings(
    convert_per_sample = function(df, config, step_id, convert_fun) {
      calls[[length(calls) + 1]] <<- step_id
    }
  )

  pipeline_convert_covRLE(list(stub), config)

  expect_identical(unlist(calls), c("covrle", "covrle_unique"))
  expect_true(step_is_done(config, "covrle", name(stub)))
})

test_that("pipeline_convert_bigwig() drives convert_per_sample() with bigwig step ids", {
  config <- fake_config(extra = list(all_mappers = TRUE, split_unique_mappers = TRUE))
  stub <- fake_experiment_stub(run_ids = c("SRR001"))
  fake_mark_all_done(config, c("aligned", "cleanbam", "exp", "ofst", "cigar_collapse",
                               "merged_lib", "covrle"), name(stub))
  calls <- list()
  testthat::local_mocked_bindings(
    convert_per_sample = function(df, config, step_id, convert_fun) {
      calls[[length(calls) + 1]] <<- step_id
    }
  )

  pipeline_convert_bigwig(list(stub), config)

  expect_identical(unlist(calls), c("bigwig", "bigwig_unique"))
  expect_true(step_is_done(config, "bigwig", name(stub)))
})
