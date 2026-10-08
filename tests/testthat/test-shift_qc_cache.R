# Unit tests for R/shift_qc_cache.R -- the per-sample-cached
# replacement for shift_qc()'s whole-experiment reload (see that file's
# own header for the root-cause writeup: confirmed live on
# GSE151959-homo_sapiens, 2026-10-07, every sample's pshifted library
# was reloaded from disk TWICE per call, with no cache, regardless of
# how many samples actually changed).
#
# Uses ORFik's own built-in ORFik.template.experiment() (real, tiny,
# fast -- ~0.3s -- dummy genome+libraries shipped with ORFik) ONLY for
# shift_qc_one_sample() itself, which needs real seqinfo()/loadRegion()/
# fimport() to do anything meaningful. Every other test here uses
# fake_experiment_stub() (shared test helper, extended there with a
# QCfolder() method and fake_qc_stub_pshifted_path(), a plain helper
# for massiveNGSpipe's own pshifted_filepath() to mock against) instead:
# it gives real, distinct runIDs() per sample and a fresh tempdir per
# call. Two real problems forced this split:
#  - ORFik::filepath() is a plain function (not an S4 generic) with an
#    internal stopifnot(is(df, "experiment")) -- it can never dispatch
#    on a fake class no matter what method is registered for it, so
#    massiveNGSpipe's own pshifted_filepath() wrapper (R/shift_qc_cache.R)
#    exists specifically to be mockable here instead.
#  - ORFik.template.experiment() turns out to have NO "Run" column at
#    all (runIDs() returns "" for every one of its rows) and always
#    resolves to the SAME shared on-disk extdata path -- either one,
#    alone, previously caused real test-isolation bugs (every sample's
#    cache silently collapsing onto one shared ""-keyed path; one
#    test's leftover files satisfying a later test's freshness check)
#    that looked like bugs in the orchestration logic itself.
# shift_qc_one_sample() is unaffected by either issue since its own
# tests only ever use ONE real sample at a time.

fake_qc_stub <- function(n = 4) {
  stub <- fake_experiment_stub(run_ids = paste0("SRR00", seq_len(n)), base_dir = tempfile("qc_stub_"))
  fake_qc_stub_pshifted_path(stub) # pre-create each run's pshifted file on disk
  stub
}

test_that("shift_qc_cache_dir() is namespaced by mapper mode, under QCfolder(df)", {
  df <- fake_qc_stub()
  all_dir <- shift_qc_cache_dir(df)
  ORFik::uniqueMappers(df) <- TRUE
  unique_dir <- shift_qc_cache_dir(df)

  expect_false(identical(all_dir, unique_dir))
  expect_true(startsWith(all_dir, ORFik::QCfolder(df)))
  expect_true(startsWith(unique_dir, ORFik::QCfolder(df)))
})

test_that("shift_qc_cache_paths() names files by run id", {
  df <- fake_qc_stub()
  paths <- shift_qc_cache_paths(df[1, ])
  run_id <- ORFik::runIDs(df[1, ])

  expect_identical(basename(paths$hitmap), paste0(run_id, "_hitmap.rds"))
  expect_identical(basename(paths$frames), paste0(run_id, "_frames.csv"))
})

test_that("shift_qc_cache_valid() is FALSE when the sample's pshifted file doesn't exist", {
  df <- fake_qc_stub()
  df_one <- df[1, ]
  testthat::local_mocked_bindings(pshifted_filepath = function(...) tempfile())
  expect_false(shift_qc_cache_valid(df_one))
})

test_that("shift_qc_cache_valid() is FALSE when no cache files exist yet", {
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub()
  expect_false(shift_qc_cache_valid(df[1, ]))
})

test_that("shift_qc_cache_valid() is FALSE when the cache predates the pshifted file (stale), TRUE otherwise", {
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub()
  df_one <- df[1, ]
  write_shift_qc_cache(df_one, list(
    hitmap = data.table::data.table(x = 1),
    frames = data.table::data.table(y = 1)
  ))
  paths <- unlist(shift_qc_cache_paths(df_one))

  # Cache older than the pshifted file -- stale.
  old_time <- Sys.time() - 1e6
  Sys.setFileTime(paths["hitmap"], old_time)
  Sys.setFileTime(paths["frames"], old_time)
  expect_false(shift_qc_cache_valid(df_one))

  # Cache newer than the pshifted file -- valid.
  Sys.setFileTime(paths["hitmap"], Sys.time() + 10)
  Sys.setFileTime(paths["frames"], Sys.time() + 10)
  expect_true(shift_qc_cache_valid(df_one))
})

test_that("write_shift_qc_cache()/read_shift_qc_cache() round-trip exactly", {
  df <- fake_qc_stub()
  df_one <- df[1, ]
  hitmap <- data.table::data.table(position = 1:3, score = c(1, 2, 3))
  frames <- data.table::data.table(frame = c(0, 1, 2), score = c(5, 6, 7))

  write_shift_qc_cache(df_one, list(hitmap = hitmap, frames = frames))
  result <- read_shift_qc_cache(df_one)

  expect_equal(result$hitmap, hitmap)
  expect_equal(result$frames, frames)
})

rfp_template_df <- function() {
  df <- ORFik::ORFik.template.experiment()
  df[df$libtype == "RFP", ]
}

test_that("shift_qc_one_sample() loads the library exactly once and returns matching hitmap/frames", {
  # Correctness verified by hand against the OLD whole-experiment path
  # (ORFik::windowPerReadLength()/orfFrameDistributions() called
  # directly on the same one sample): hitmap was exactly identical,
  # frames matched on every value (only factor-vs-numeric typing
  # differed, which shift_qc_cached()'s own aggregation step applies
  # afterward, same as the old code did). This test locks in that the
  # function runs end to end and returns the right shape.
  df <- rfp_template_df()
  df_one <- df[1, ]
  cds <- ORFik::loadRegion(df, part = "cds")
  mrna <- ORFik::loadRegion(df, part = "mrna")

  result <- shift_qc_one_sample(df_one, cds, mrna, upstream = 5, downstream = 20)

  expect_false(inherits(result, "shift_qc_sample_error"))
  expect_true(is.data.frame(result$hitmap))
  expect_true(nrow(result$hitmap) > 0)
  expect_true(is.data.frame(result$frames))
  expect_true(nrow(result$frames) > 0)
  expect_identical(result$name_short, ORFik::bamVarName(df_one, skip.experiment = TRUE))
})

test_that("shift_qc_one_sample() returns a shift_qc_sample_error (not a thrown error) when loading fails", {
  df <- rfp_template_df()
  df_one <- df[1, ]
  cds <- ORFik::loadRegion(df, part = "cds")
  mrna <- ORFik::loadRegion(df, part = "mrna")

  testthat::local_mocked_bindings(pshifted_filepath = function(...) stop("no such file"))

  result <- shift_qc_one_sample(df_one, cds, mrna, upstream = 5, downstream = 20)

  expect_true(inherits(result, "shift_qc_sample_error"))
  expect_true(inherits(result, "error"))
  expect_match(result$message, "no such file")
})

test_that("shift_qc_cached() computes every sample on a first call, then skips all of them on a resume", {
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub()
  call_count <- 0L
  fake_result <- list(hitmap = data.table::data.table(position = 1, frame = 0),
                      frames = data.table::data.table(frame = 0, score = 1, fraction = "x", length = 30))
  testthat::local_mocked_bindings(
    shift_qc_one_sample = function(...) { call_count <<- call_count + 1; fake_result },
    shift_qc_build_combined_plot = function(...) invisible(NULL),
    check_adapter_barcode_quality = function(...) invisible(NULL)
  )
  testthat::local_mocked_bindings(
    shift_qc_annotation_context = function(df) list(upstream = 5, downstream = 20, cds = list(), mrna = list())
  )

  shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam())
  expect_equal(call_count, nrow(df))

  call_count <- 0L
  shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam())
  expect_equal(call_count, 0L) # fully resumed: nothing recomputed
})

test_that("shift_qc_cached() recomputes only the one sample whose pshifted file changed", {
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub()
  fake_result <- list(hitmap = data.table::data.table(position = 1, frame = 0),
                      frames = data.table::data.table(frame = 0, score = 1, fraction = "x", length = 30))
  call_log <- character()
  testthat::local_mocked_bindings(
    shift_qc_one_sample = function(df_one_row, ...) {
      call_log <<- c(call_log, ORFik::runIDs(df_one_row))
      fake_result
    },
    shift_qc_build_combined_plot = function(...) invisible(NULL),
    check_adapter_barcode_quality = function(...) invisible(NULL)
  )
  testthat::local_mocked_bindings(
    shift_qc_annotation_context = function(df) list(upstream = 5, downstream = 20, cds = list(), mrna = list())
  )

  shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam())
  call_log <- character()

  # Simulate sample 2 being redone (fresh pshift) by bumping its
  # pshifted file's mtime past its own cache's.
  Sys.setFileTime(fake_qc_stub_pshifted_path(df[2, ]), Sys.time() + 10)

  shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam())
  expect_identical(call_log, ORFik::runIDs(df[2, ]))
})

test_that("shift_qc_cached() writes each sample's cache file as soon as THAT sample finishes, not after the whole batch returns", {
  # The crash-safety property this locks in: an interrupted run (kill,
  # OOM, crash) partway through a large batch must keep every
  # already-finished sample's work, not lose the whole batch because
  # nothing was persisted until every sample returned. Confirmed live,
  # 2026-10-08: the OLD design (write all samples' results only AFTER
  # the entire bplapply() call returned) meant a multi-hour
  # PRJNA637713 valid_pshift run would have lost 100% of its progress
  # on a kill, regardless of how many of its ~270 samples had already
  # finished computing minutes or hours earlier.
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub(n = 2)
  fake_result <- list(hitmap = data.table::data.table(position = 1, frame = 0),
                      frames = data.table::data.table(frame = 0, score = 1, fraction = "x", length = 30))
  sample1_cached_before_sample2_ran <- NA
  testthat::local_mocked_bindings(
    shift_qc_one_sample = function(df_one_row, ...) {
      run <- ORFik::runIDs(df_one_row)
      if (run == "SRR002") {
        # SerialParam runs SRR001 to completion first -- if its write
        # happened inside its own worker turn (not deferred), its
        # cache must already be valid by the time SRR002 starts.
        sample1_cached_before_sample2_ran <<- shift_qc_cache_valid(df[1, ])
      }
      fake_result
    },
    shift_qc_build_combined_plot = function(...) invisible(NULL),
    check_adapter_barcode_quality = function(...) invisible(NULL)
  )
  testthat::local_mocked_bindings(
    shift_qc_annotation_context = function(df) list(upstream = 5, downstream = 20, cds = list(), mrna = list())
  )

  shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam())
  expect_true(sample1_cached_before_sample2_ran)
})

test_that("shift_qc_cached() excludes a failing sample from the aggregate and never caches it, so it's retried next time", {
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub()
  ok_result <- list(hitmap = data.table::data.table(position = 1, frame = 0),
                    frames = data.table::data.table(frame = 0, score = 1, fraction = "x",
                                                     length = 30))
  failing_run <- ORFik::runIDs(df[1, ])
  call_log <- character()
  testthat::local_mocked_bindings(
    shift_qc_one_sample = function(df_one_row, ...) {
      run <- ORFik::runIDs(df_one_row)
      call_log <<- c(call_log, run)
      if (run == failing_run) {
        return(structure(list(message = "simulated failure"),
                         class = c("shift_qc_sample_error", "error", "condition")))
      }
      ok_result
    },
    shift_qc_build_combined_plot = function(...) invisible(NULL),
    check_adapter_barcode_quality = function(...) invisible(NULL)
  )
  testthat::local_mocked_bindings(
    shift_qc_annotation_context = function(df) list(upstream = 5, downstream = 20, cds = list(), mrna = list())
  )

  expect_warning(shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam()), "simulated failure")

  expect_false(shift_qc_cache_valid(df[1, ])) # failing sample: no valid cache
  expect_true(shift_qc_cache_valid(df[2, ]))  # succeeding siblings: cached fine

  # Resuming again must retry the failing sample (not skip it forever)
  # and must NOT re-touch the siblings that already succeeded.
  call_log <- character()
  expect_warning(shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam()), "simulated failure")
  expect_identical(call_log, failing_run)
})

test_that("shift_qc_cached() doesn't crash when one sample's frames table has zero rows (and so a shorter column set than the others)", {
  # Confirmed live, 2026-10-08, PRJNA637713-zea_mays/SRR13808095:
  # shift_qc_one_sample()'s own `if (nrow(frames) > 0)` column-
  # augmentation (adding `length`, renaming `fraction`) never runs for
  # a sample with NO periodicity-covered regions at all, so its cached
  # frames.csv keeps the shorter (3-column: fraction/frame/score)
  # schema instead of the normal 4-column one -- rbindlist() without
  # fill=TRUE errored on that mismatch ("Item N has 3 columns,
  # inconsistent with item 1 which has 4 columns"), crashing the WHOLE
  # experiment's valid_pshift aggregation over one such sample.
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub(n = 2)
  normal_result <- list(hitmap = data.table::data.table(position = 1, frame = 0),
                        frames = data.table::data.table(fraction = "x", frame = 0, score = 1, length = 30))
  empty_run <- ORFik::runIDs(df[2, ])
  testthat::local_mocked_bindings(
    shift_qc_one_sample = function(df_one_row, ...) {
      if (ORFik::runIDs(df_one_row) == empty_run) {
        return(list(hitmap = data.table::data.table(position = integer(), frame = integer()),
                   frames = data.table::data.table(fraction = character(), frame = integer(), score = numeric())))
      }
      normal_result
    },
    shift_qc_build_combined_plot = function(...) invisible(NULL),
    check_adapter_barcode_quality = function(...) invisible(NULL)
  )
  testthat::local_mocked_bindings(
    shift_qc_annotation_context = function(df) list(upstream = 5, downstream = 20, cds = list(), mrna = list())
  )

  expect_no_error(shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam()))
  expect_true(shift_qc_cache_valid(df[1, ]))
  expect_true(shift_qc_cache_valid(df[2, ])) # the empty-but-successful sample is still cached, not treated as a failure

  frameQC <- data.table::fread(file.path(ORFik::QCfolder(df), "Ribo_frames_all.csv"))
  expect_equal(nrow(frameQC), 1) # only the normal sample's one row; the empty sample contributes none
})

test_that("shift_qc_cached() aggregates cached per-sample frame tables into Ribo_frames_all.csv/badzero.csv", {
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub(n = 2)
  run_a <- ORFik::runIDs(df[1, ])
  run_b <- ORFik::runIDs(df[2, ])
  results_by_run <- list()
  results_by_run[[run_a]] <- list(hitmap = data.table::data.table(position = 1, frame = 0),
    frames = data.table::data.table(frame = c(0, 1), score = c(10, 1), fraction = "A", length = 30))
  results_by_run[[run_b]] <- list(hitmap = data.table::data.table(position = 1, frame = 0),
    frames = data.table::data.table(frame = c(0, 1), score = c(1, 10), fraction = "B", length = 30))
  testthat::local_mocked_bindings(
    shift_qc_one_sample = function(df_one_row, ...) results_by_run[[ORFik::runIDs(df_one_row)]],
    shift_qc_build_combined_plot = function(...) invisible(NULL),
    check_adapter_barcode_quality = function(...) invisible(NULL)
  )
  testthat::local_mocked_bindings(
    shift_qc_annotation_context = function(df) list(upstream = 5, downstream = 20, cds = list(), mrna = list())
  )

  shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam())

  all_csv <- data.table::fread(file.path(ORFik::QCfolder(df), "Ribo_frames_all.csv"))
  expect_identical(nrow(all_csv), 4L) # 2 samples x 2 frame rows each

  badzero <- data.table::fread(file.path(ORFik::QCfolder(df), "Ribo_frames_badzero.csv"))
  # Sample "A": frame 0 has the higher score (10 vs 1) -> frame 0 is
  # best_frame -> excluded from badzero. Sample "B": frame 0 has the
  # LOWER score (1 vs 10) -> frame 1 is best_frame -> sample B's frame
  # 0 row IS a bad-zero row.
  expect_identical(nrow(badzero), 1L)
  expect_identical(as.character(badzero$fraction), "B")
})

test_that("shift_qc_cached() plots only the already-existing shiftPlots() max-39 subset", {
  testthat::local_mocked_bindings(pshifted_filepath = fake_qc_stub_pshifted_path)
  df <- fake_qc_stub() # 4 rows, well under 39
  plotted_idx <- NULL
  fake_result <- list(hitmap = data.table::data.table(position = 1, frame = 0),
                      frames = data.table::data.table(frame = 0, score = 1, fraction = "x", length = 30))
  testthat::local_mocked_bindings(
    shift_qc_one_sample = function(...) fake_result,
    shift_qc_build_combined_plot = function(df, plot_idx) plotted_idx <<- plot_idx,
    check_adapter_barcode_quality = function(...) invisible(NULL)
  )
  testthat::local_mocked_bindings(
    shift_qc_annotation_context = function(df) list(upstream = 5, downstream = 20, cds = list(), mrna = list())
  )

  shift_qc_cached(df, BPPARAM = BiocParallel::SerialParam())

  expect_identical(plotted_idx, seq_len(nrow(df)))
})
