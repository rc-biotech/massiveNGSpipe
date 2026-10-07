# Unit tests for R/pipeline_sample_flags.R -- per-sample progress/resume
# markers (fetch/trim/align). Pure tempdir()-based file I/O, no network/
# STAR/fastp. This is the mechanism the interrupt-and-resume feature and
# the progress checklist both depend on -- see md-level notes in the
# source file for the design rationale.

test_that("set_sample_flag / n_samples_done / samples_done round-trip", {
  config <- fake_config()
  exp <- "study1-organism"
  expect_identical(n_samples_done(config, "aligned", exp), 0L)
  expect_identical(samples_done(config, "aligned", exp), character())

  set_sample_flag(config, "aligned", exp, "SRR001")
  set_sample_flag(config, "aligned", exp, "SRR002")
  expect_identical(n_samples_done(config, "aligned", exp), 2L)
  expect_setequal(samples_done(config, "aligned", exp), c("SRR001", "SRR002"))
})

test_that("set_sample_flag stores arbitrary values, not just TRUE", {
  config <- fake_config()
  exp <- "study1-organism"
  row <- data.table::data.table(Run = "SRR001", adapter = "AGATCGGAAGAG")
  set_sample_flag(config, "trim", exp, "SRR001", value = row)

  values <- sample_flag_values(config, "trim", exp)
  expect_length(values, 1)
  expect_identical(values[[1]], row)
})

test_that("sample_flag_values reads back every marker, old and new alike", {
  config <- fake_config()
  exp <- "study1-organism"
  set_sample_flag(config, "trim", exp, "SRR001", value = data.table::data.table(x = 1))
  set_sample_flag(config, "trim", exp, "SRR002", value = data.table::data.table(x = 2))
  combined <- data.table::rbindlist(sample_flag_values(config, "trim", exp))
  expect_identical(sort(combined$x), c(1, 2))
})

test_that("sample_flag_values names its result by run id", {
  # A caller needs to know WHICH run a marker belongs to in order to
  # look up that run's own row elsewhere (e.g. pipeline_trim()'s
  # adapter_barcode_table.csv reconstruction falling back to an
  # existing on-disk row for a legacy plain-TRUE marker) -- confirmed
  # live, 2026-10-07: this used to return an unnamed list, silently
  # breaking exactly that kind of lookup.
  config <- fake_config()
  exp <- "study1-organism"
  set_sample_flag(config, "trim", exp, "SRR001", value = data.table::data.table(x = 1))
  set_sample_flag(config, "trim", exp, "SRR002", value = data.table::data.table(x = 2))

  values <- sample_flag_values(config, "trim", exp)

  expect_setequal(names(values), c("SRR001", "SRR002"))
  expect_identical(values[["SRR002"]]$x, 2)
})

test_that("reset_sample_flags clears markers for one experiment only", {
  config <- fake_config()
  set_sample_flag(config, "aligned", "exp1", "SRR001")
  set_sample_flag(config, "aligned", "exp2", "SRR001")

  reset_sample_flags(config, "aligned", "exp1")
  expect_identical(n_samples_done(config, "aligned", "exp1"), 0L)
  expect_identical(n_samples_done(config, "aligned", "exp2"), 1L) # untouched
})

test_that("reset_sample_flags on a never-created dir is a silent no-op", {
  config <- fake_config()
  expect_silent(reset_sample_flags(config, "aligned", "never_existed"))
})

test_that("n_samples_done only counts .rds files, ignoring stray other files", {
  config <- fake_config()
  exp <- "study1-organism"
  set_sample_flag(config, "aligned", exp, "SRR001")
  d <- file.path(config$project, "sample_flags", "aligned", exp)
  writeLines("not a marker", file.path(d, "notes.txt"))
  expect_identical(n_samples_done(config, "aligned", exp), 1L)
})

test_that("samples_done/n_samples_done/sample_flag_values on a nonexistent dir all return empty, not error", {
  config <- fake_config()
  expect_identical(n_samples_done(config, "aligned", "no_such_exp"), 0L)
  expect_identical(samples_done(config, "aligned", "no_such_exp"), character())
  expect_identical(sample_flag_values(config, "aligned", "no_such_exp"), list())
})
