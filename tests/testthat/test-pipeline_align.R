# Unit tests for pipeline_align() and pipeline_cleanup()
# (R/pipeline_preset_steps_sub.R) -- both previously had NO
# per-sample-resume-safe file resolution: run_files_organizer()/
# match_bam_to_metadata() were called against the FULL study before any
# not-done filtering, which fails the moment an already-finished
# sibling's input file is gone (delete_trimmed_files/delete_collapsed_files),
# and pipeline_cleanup() additionally had no resume filtering at all,
# letting fs::file_move() try to move an already-renamed sibling's BAM
# onto itself -- which DELETES it instead of no-op'ing. Confirmed live
# and fixed on PRJNA926112-homo_sapiens, 2026-10-05/06; these tests lock
# that fix in so it can't silently regress.

test_that("pipeline_align() only resolves/aligns samples not yet marked done, and only deletes THEIR input files", {
  config <- fake_config(preset = "Ribo-seq", extra = list(delete_collapsed_files = TRUE))
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(
    bam_dir = bam_dir,
    runs = data.table::data.table(Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp_name <- "PRJNA000001-homo_sapiens"
  fake_mark_all_done(config, c("trim", "collapsed"), exp_name)
  trimmed_single_dir <- file.path(bam_dir, "trim", "SINGLE")
  dir.create(trimmed_single_dir, recursive = TRUE)
  # Only SRR002's collapsed input still exists -- SRR001 is already
  # aligned and its own input was already cleaned up earlier.
  srr002_input <- file.path(trimmed_single_dir, "collapsed_trimmed_SRR002.fasta.gz")
  file.create(srr002_input)
  set_sample_flag(config, "aligned", exp_name, "SRR001")

  resolved_for <- character()
  testthat::local_mocked_bindings(
    run_files_organizer = function(runs, ...) {
      resolved_for <<- runs$Run
      stats::setNames(as.list(file.path(trimmed_single_dir, paste0("collapsed_trimmed_", runs$Run, ".fasta.gz"))), NULL)
    }
  )
  testthat::local_mocked_bindings(
    STAR.install = function(...) "star", install.fastp = function(...) "fastp",
    STAR.align.single = function(...) invisible(NULL), .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    run_experiment_subprocess = function(func, args, ...) do.call(func, args)
  )
  checked_runs <- NULL
  checked_pairs <- NULL
  testthat::local_mocked_bindings(
    alignment_final_checks = function(input_dir, output_dir, runs, pairs, config, steps) {
      checked_runs <<- runs$Run
      checked_pairs <<- pairs
    }
  )

  pipeline_align(pipelines[["PRJNA000001"]], config)

  expect_identical(resolved_for, "SRR002")
  expect_identical(checked_runs, "SRR002")
  expect_setequal(samples_done(config, "aligned", exp_name), c("SRR001", "SRR002"))
  # What alignment_final_checks() receives as `pairs` is exactly what its
  # own delete_collapsed_files branch later deletes (unlist(pairs)) --
  # confirming this is SRR002-only is what makes that deletion safe.
  expect_identical(unlist(checked_pairs, use.names = FALSE), srr002_input)
})

test_that("pipeline_align()'s delete_collapsed_files only deletes the input files it actually used, never the whole directory", {
  bam_dir <- tempfile("bam_")
  trimmed_single_dir <- file.path(bam_dir, "trim", "SINGLE")
  aligned_dir <- file.path(bam_dir, "aligned")
  dir.create(trimmed_single_dir, recursive = TRUE)
  dir.create(aligned_dir, recursive = TRUE)
  used_file <- file.path(trimmed_single_dir, "used.fasta.gz")
  sibling_file <- file.path(trimmed_single_dir, "sibling_untouched.fasta.gz")
  writeLines("x", used_file); writeLines("x", sibling_file)
  writeLines("x", file.path(aligned_dir, "SRR001.bam")) # non-empty, matches the one run being checked

  testthat::local_mocked_bindings(
    system2 = function(...) invisible(NULL), .package = "base"
  )
  testthat::local_mocked_bindings(
    check_alignment_rate = function(...) invisible(NULL)
  )

  alignment_final_checks(
    input_dir = trimmed_single_dir, output_dir = bam_dir,
    runs = data.table::data.table(Run = "SRR001"),
    pairs = list(used_file),
    config = list(min_alignment_rate_pshift = 0, delete_collapsed_files = TRUE),
    steps = NULL # ORFik::STAR.allsteps.multiQC()'s own first line short-circuits on NULL steps; irrelevant to what this test checks (the deletion)
  )

  # The only file deleted must be the one actually resolved/used for
  # THIS sample -- never the sibling's, even though both sat in the
  # same directory (fs::dir_delete(input_dir) used to delete both).
  expect_false(file.exists(used_file))
  expect_true(file.exists(sibling_file))
})

test_that("pipeline_cleanup() skips renaming a sample whose BAM already has its final <Run>.bam name", {
  config <- fake_config(preset = "Ribo-seq")
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(
    bam_dir = bam_dir,
    runs = data.table::data.table(Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp_name <- "PRJNA000001-homo_sapiens"
  fake_mark_all_done(config, c("trim", "collapsed", "aligned"), exp_name)
  aligned_dir <- file.path(bam_dir, "aligned")
  dir.create(aligned_dir, recursive = TRUE)
  # SRR001 already has its final name (already cleaned earlier).
  # SRR002 still has STAR's native output name (needs cleanup now).
  file.create(file.path(aligned_dir, "SRR001.bam"))
  file.create(file.path(aligned_dir, "SRR002_Aligned.sortedByCoord.out.bam"))

  moved <- list()
  testthat::local_mocked_bindings(
    file_move = function(from, to) { moved[[length(moved) + 1]] <<- list(from = from, to = to) },
    .package = "fs"
  )

  pipeline_cleanup(pipelines[["PRJNA000001"]], config)

  expect_length(moved, 1)
  expect_identical(basename(moved[[1]]$from), "SRR002_Aligned.sortedByCoord.out.bam")
  expect_identical(basename(moved[[1]]$to), "SRR002.bam")
  # SRR001's already-correct file must never even be considered for a move.
  expect_true(file.exists(file.path(aligned_dir, "SRR001.bam")))
})
