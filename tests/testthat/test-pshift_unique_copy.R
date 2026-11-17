# Unique-mapper offsets are copied from the all-mappers table, never detected (Håkon 2026-11-17; regression from 26b3bc1).
local_unique_layout <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  for (d in c("aligned", "aligned/pshifted", "aligned/unique_mappers", "aligned/unique_mappers/pshifted"))
    dir.create(file.path(root, d), showWarnings = FALSE)
  stub <- fake_experiment_stub(run_ids = c("SRR001", "SRR002"))
  folder <- function(df) file.path(root, "aligned", if (isTRUE(uniqueMappers(df))) "unique_mappers" else "")
  testthat::local_mocked_bindings(
    ofst_filepath = function(df) file.path(folder(df), paste0(runIDs(df), ".ofst")),
    pshifted_filepath = function(df) file.path(folder(df), "pshifted", paste0(runIDs(df), "_pshifted.ofst")),
    .env = env)
  stub_u <- stub; uniqueMappers(stub_u) <- TRUE
  for (f in c(ofst_filepath(stub), ofst_filepath(stub_u))) writeLines("raw", f)
  off <- function(o) data.frame(fraction = c(28L, 29L), offsets_start = o)
  all_tab <- setNames(list(off(c(-12L, -12L)), off(c(-12L, -13L))), ofst_filepath(stub))
  saveRDS(all_tab, file.path(root, "aligned", "pshifted", "shifting_table.rds"))
  list(root = root, stub = stub, stub_u = stub_u, all_tab = all_tab, off = off,
       unique_table = file.path(root, "aligned", "unique_mappers", "pshifted", "shifting_table.rds"))
}

test_that("all-mappers offsets are staged into the unique table when unique offsets differ", {
  x <- local_unique_layout()
  # unique table detected separately (the regression): SRR002 length 29 differs by 1 nt
  saveRDS(setNames(list(x$off(c(-12L, -12L)), x$off(c(-12L, -12L))), ofst_filepath(x$stub_u)), x$unique_table)
  for (f in pshifted_filepath(x$stub_u)) { Sys.sleep(0.01); writeLines("shifted", f) }
  staged <- pshift_copy_all_shifts_to_unique(x$stub, x$stub_u)
  # SRR001 equal + legacy-valid: untouched; SRR002 copied
  expect_equal(basename(staged), "SRR002.ofst")
  tab <- readRDS(x$unique_table)
  expect_identical(pshift_canonical_offsets(tab[[ofst_filepath(x$stub_u)[2]]]), x$off(c(-12L, -13L)))
  expect_length(tab, 2)
  # the copy is an explicit request: the sample is no longer valid and reuse returns the copied offsets
  expect_false(pshift_sample_valid(x$stub_u[2, ]))
  recipe <- pshift_detection_recipe(28:29, character(), 80)
  expect_identical(pshift_reusable_offsets(x$stub_u[2, ], recipe, reuse = FALSE), x$off(c(-12L, -13L)))
  # running it again does not re-stage SRR001 and keeps SRR002 pending
  expect_equal(basename(pshift_copy_all_shifts_to_unique(x$stub, x$stub_u)), "SRR002.ofst")
})

test_that("a missing unique table is created from the all-mappers table", {
  x <- local_unique_layout()
  staged <- pshift_copy_all_shifts_to_unique(x$stub, x$stub_u)
  expect_length(staged, 2)
  tab <- readRDS(x$unique_table)
  expect_identical(lapply(unname(tab), pshift_canonical_offsets), lapply(unname(x$all_tab), pshift_canonical_offsets))
})

test_that("copy stops when the all-mappers table or a Run entry is missing", {
  x <- local_unique_layout()
  saveRDS(x$all_tab[1], file.path(x$root, "aligned", "pshifted", "shifting_table.rds"))
  expect_error(pshift_copy_all_shifts_to_unique(x$stub, x$stub_u), "no all-mappers offsets for SRR002")
  unlink(file.path(x$root, "aligned", "pshifted", "shifting_table.rds"))
  expect_error(pshift_copy_all_shifts_to_unique(x$stub, x$stub_u), "no all-mappers shifting table")
})

test_that("pipeline_pshift() unique pass copies offsets and never runs the extra-length hook", {
  config <- fake_config(preset = "Ribo-seq", extra = list(all_mappers = TRUE, split_unique_mappers = TRUE,
                                     reuse_shifts_if_existing = TRUE, accepted_lengths_rpf = 28:29,
                                     max_no_adapter_removed_pct = 80, pshift_extra_lengths = 34:35))
  stub <- fake_experiment_stub(run_ids = "SRR001")
  fake_mark_all_done(config, c("trim", "collapsed", "aligned", "cleanbam", "exp", "ofst"), name(stub))
  calls <- character()
  testthat::local_mocked_bindings(
    pshift_needed_subset = function(df, ...) { calls <<- c(calls, paste0("shift:", uniqueMappers(df))); 1 },
    pshift_extra_lengths_hook = function(df, res, ...) { calls <<- c(calls, paste0("extra:", uniqueMappers(df))); res },
    pshift_copy_all_shifts_to_unique = function(df_all, df_unique) {
      calls <<- c(calls, paste0("copy:", uniqueMappers(df_all), uniqueMappers(df_unique))); invisible(character()) })
  pipeline_pshift(list(stub), config)
  expect_equal(calls, c("shift:FALSE", "extra:FALSE", "copy:FALSETRUE", "shift:TRUE"))
  expect_true(step_is_done(config, "pshifted", name(stub)))
})
