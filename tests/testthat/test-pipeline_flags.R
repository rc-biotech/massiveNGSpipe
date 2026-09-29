# Unit tests for R/pipeline_flags.R -- massiveNGSpipe's own experiment-level
# flag bookkeeping (which step is done for which study/organism). All
# file-based tests run against tempdir(), no real STAR/fastp/network.

test_that("libtype_flags builds the expected flag recipe per preset", {
  rfp <- libtype_flags("Ribo-seq", "online", contam = FALSE)
  expect_true(all(c("pipe_fetch", "pipe_trim_collapse", "pipe_align_clean",
                    "pipe_exp_ofst", "pipe_pshift_and_validate",
                    "pipe_merge_study", "pipe_convert", "pipe_counts") %in% names(rfp)))
  expect_true("pshifted" %in% rfp)
  expect_true("valid_pshift" %in% rfp)

  rna <- libtype_flags("RNA-seq", "online", contam = FALSE)
  expect_true("cigar_collapse" %in% rna)
  expect_false("trim" %in% rna) # RNA-seq preset clears trim_flags entirely

  empty <- libtype_flags("empty", "online", contam = FALSE)
  expect_true(length(empty) > 0) # still has download/align/exp/merge/format/count flags
  expect_false("pshifted" %in% empty)
})

test_that("libtype_flags: mode='local' drops the fetch/start flags", {
  online <- libtype_flags("RNA-seq", "online")
  local <- libtype_flags("RNA-seq", "local")
  expect_true("start" %in% online)
  expect_false("start" %in% local)
  expect_false("fetch" %in% local)
})

test_that("libtype_flags: contam=TRUE prepends a contam flag to align_clean", {
  no_contam <- libtype_flags("Ribo-seq", contam = FALSE)
  with_contam <- libtype_flags("Ribo-seq", contam = TRUE)
  expect_false("contam" %in% no_contam)
  expect_true("contam" %in% with_contam)
  expect_identical(unname(with_contam[names(with_contam) == "pipe_align_clean"])[1], "contam")
})

test_that("libtype_flags: invalid preset errors with a clean, readable message", {
  expect_error(libtype_flags("not-a-real-preset"),
              regexp = "Currently valid preset pipelines are of types: Ribo-seq, RNA-seq, disome, SSU, empty",
              fixed = TRUE)
})

test_that("pipeline_flags creates real directories and a grouping attribute", {
  project <- tempfile("mNGSp_test_")
  flags <- pipeline_flags(project, mode = "local", preset = "RNA-seq", create_dirs = TRUE)
  expect_true(all(dir.exists(flags)))
  expect_identical(names(attributes(flags))[1], "names")
  expect_true(!is.null(attr(flags, "grouping")))
  expect_identical(length(attr(flags, "grouping")), length(flags))
})

test_that("pipeline_flags(preset='empty') returns character(0) with no dirs made", {
  project <- tempfile("mNGSp_test_")
  flags <- pipeline_flags(project, mode = "local", preset = "empty", create_dirs = TRUE)
  expect_identical(flags, character())
  expect_false(dir.exists(file.path(project, "flags")))
})

test_that("step_is_done reflects real flag files, vectorized over experiment", {
  config <- fake_config()
  exp <- "study1-organism"
  expect_false(step_is_done(config, "aligned", exp))
  set_flag(config, "aligned", exp)
  expect_true(step_is_done(config, "aligned", exp))

  exp2 <- "study2-organism"
  expect_identical(step_is_done(config, "aligned", c(exp, exp2)), c(TRUE, FALSE))
})

test_that("set_flag errors for a step whose directory doesn't exist", {
  config <- fake_config()
  expect_error(set_flag(config, "not_a_real_step", "exp1"))
})

test_that("step_is_next_not_done: fetch must be checked via step_is_done directly", {
  config <- fake_config()
  expect_error(step_is_next_not_done(config, "fetch", "exp1"),
              "fetch should use step_is_done directly", fixed = TRUE)
})

test_that("step_is_next_not_done: first step has no predecessor requirement", {
  config <- fake_config()
  first_step <- names(config$flag)[1]
  expect_true(step_is_next_not_done(config, first_step, "exp1"))
  set_flag(config, first_step, "exp1")
  expect_false(step_is_next_not_done(config, first_step, "exp1"))
})

test_that("step_is_next_not_done: a later step needs all earlier steps done first", {
  config <- fake_config()
  steps <- names(config$flag)
  later_step <- steps[length(steps)]
  expect_false(step_is_next_not_done(config, later_step, "exp1")) # nothing done yet
  for (s in steps[-length(steps)]) set_flag(config, s, "exp1")
  expect_true(step_is_next_not_done(config, later_step, "exp1"))
})

test_that("step_is_next_not_done errors for an unknown step name", {
  config <- fake_config()
  expect_error(step_is_next_not_done(config, "totally_bogus_step", "exp1"))
})

test_that("remove_flag deletes an existing flag and is silent on a missing one by default", {
  config <- fake_config()
  set_flag(config, "aligned", "exp1")
  expect_true(step_is_done(config, "aligned", "exp1"))
  remove_flag(config, "aligned", "exp1")
  expect_false(step_is_done(config, "aligned", "exp1"))
  expect_silent(remove_flag(config, "aligned", "exp1")) # already gone, warning=FALSE default
})

test_that("set_flag_all_exp / remove_flag_all_exp round-trip, including steps='all'", {
  config <- fake_config()
  steps <- names(config$flag)
  set_flag_all_exp(config, steps = "all", exps = "exp1")
  expect_true(all(step_is_done(config, steps, "exp1")))
  remove_flag_all_exp(config, steps = "all", exps = "exp1")
  expect_true(all(!step_is_done(config, steps, "exp1")))
})

test_that("remove_flag_all_exp rejects an invalid step name", {
  config <- fake_config()
  expect_error(remove_flag_all_exp(config, steps = "not_a_step", exps = "exp1"))
})

test_that("remove_flag_all_exp_from resets a step and everything after it, not before", {
  config <- fake_config()
  steps <- names(config$flag)
  skip_if(length(steps) < 3, "preset doesn't have enough steps for this test")
  for (s in steps) set_flag(config, s, "exp1")

  mid <- steps[ceiling(length(steps) / 2)]
  remove_flag_all_exp_from(config, mid, "exp1")

  before_mid <- steps[seq_len(match(mid, steps) - 1)]
  from_mid <- steps[match(mid, steps):length(steps)]
  expect_true(all(step_is_done(config, before_mid, "exp1")))
  expect_true(all(!step_is_done(config, from_mid, "exp1")))
})

test_that("flag_grouping / preset_grouping round-trip a flags vector", {
  # names(flags) are step ids (as pipeline_flags() produces); the
  # separate "grouping" attribute holds the owning function/group name
  # for each step id -- these are two different things, not the same
  # vector reused (a real trap: pipeline_flags() explicitly overwrites
  # names(flags) with the step ids and stashes the original group names
  # in this attribute instead).
  flags <- c(start = "path/start", fetch = "path/fetch", trim = "path/trim")
  attr(flags, "grouping") <- c("pipe_fetch", "pipe_fetch", "pipe_trim_collapse")
  grouped <- flag_grouping(flags)
  expect_identical(grouped$pipe_fetch, c("start", "fetch"))
  expect_identical(grouped$pipe_trim_collapse, "trim")
})

test_that("preset_grouping errors clearly on malformed flags vectors", {
  no_names <- c("start", "fetch")
  expect_error(preset_grouping(no_names), "flags must have names", fixed = TRUE)

  flags <- c(a = "start", a = "fetch")
  expect_error(preset_grouping(flags), regexp = "grouping")
})

test_that("add_step_to_pipeline appends a new step and updates config$flag_steps", {
  config <- fake_config(preset = "empty")
  dummy_step_fun <- function(pipelines, config) invisible(NULL)
  # add_step_to_pipeline resolves the function by name via get(), so it
  # must exist in an environment get() can find -- assign into globalenv
  # for the duration of this test. group_name is passed explicitly here
  # rather than relying on the group_name=name_of_function(FUN) default,
  # which captures the literal argument expression as written at
  # add_step_to_pipeline's OWN call site inside its default -- i.e. it
  # resolves to the string "FUN" (its own parameter name), not whatever
  # symbol the outer caller passed. Not something to depend on in a test.
  assign("dummy_step_fun", dummy_step_fun, envir = .GlobalEnv)
  on.exit(rm("dummy_step_fun", envir = .GlobalEnv), add = TRUE)

  config2 <- add_step_to_pipeline(config, "custom_step", dummy_step_fun,
                                  group_name = "dummy_step_fun")
  expect_true("custom_step" %in% names(config2$flag))
  expect_true(dir.exists(config2$flag["custom_step"]))
  expect_true("dummy_step_fun" %in% names(config2$flag_steps))
})

test_that("add_step_to_pipeline rejects duplicate short_name/group_name", {
  config <- fake_config()
  existing_step <- names(config$flag)[1]
  f <- function(pipelines, config) NULL
  expect_error(add_step_to_pipeline(config, existing_step, f))
})
