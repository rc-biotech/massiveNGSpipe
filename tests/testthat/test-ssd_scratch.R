# Unit tests for resolve_align_scratch() (R/pipeline_preset_steps_sub.R)
# and its wiring into pipeline_align_one_organism()'s per-sample loop --
# the opportunistic SSD scratch copy / copy-back for STAR alignment,
# built from the 2026-10-09 benchmark findings (see
# ~/shared_workspace/hakon_playground/server_and_software_benchmarks/
# 2026-10-06_star_alignment_ssd_vs_main_drive/FINDINGS.md). Purely
# opportunistic: every failure mode here must fall back to the main
# drive, never fail the sample itself.

test_that("resolve_align_scratch() is a no-op when ssd_scratch_dir is NULL", {
  config <- fake_config() # ssd_scratch_dir = NULL by default
  file1 <- tempfile("fastq_"); file.create(file1)
  output_dir <- tempfile("bam_")

  result <- resolve_align_scratch(file1, output_dir, "SRR001", config)

  expect_false(result$used_ssd)
  expect_identical(result$output.dir, output_dir)
  expect_identical(result$file1, file1)
})

test_that("resolve_align_scratch() uses SSD scratch when free space is sufficient, and cleanup() removes it", {
  ssd_dir <- tempfile("ssd_scratch_"); dir.create(ssd_dir)
  config <- fake_config(extra = list(ssd_scratch_dir = ssd_dir, ssd_min_free_gb = 1))
  file1 <- tempfile("fastq_"); writeLines("x", file1) # tiny -- needed_gb is tiny too
  output_dir <- file.path(tempfile("bam_"))

  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake", .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Drive_Free = "999G"), .package = "ORFik"
  )

  result <- resolve_align_scratch(file1, output_dir, "SRR001", config)

  expect_true(result$used_ssd)
  expect_true(file.exists(result$file1))
  expect_identical(dirname(result$file1), result$output.dir)
  expect_true(startsWith(result$output.dir, ssd_dir))

  result$cleanup()
  expect_false(dir.exists(result$output.dir))
})

test_that("resolve_align_scratch() falls back to the main drive when free space is insufficient, without attempting a copy", {
  ssd_dir <- tempfile("ssd_scratch_"); dir.create(ssd_dir)
  config <- fake_config(extra = list(ssd_scratch_dir = ssd_dir, ssd_min_free_gb = 1e6)) # impossible margin
  file1 <- tempfile("fastq_"); writeLines("x", file1)
  output_dir <- tempfile("bam_")

  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake", .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Drive_Free = "10G"), .package = "ORFik"
  )

  result <- resolve_align_scratch(file1, output_dir, "SRR001", config)

  expect_false(result$used_ssd)
  expect_identical(result$output.dir, output_dir)
  expect_length(list.files(ssd_dir, recursive = TRUE), 0) # nothing was copied
})

test_that("resolve_align_scratch() falls back and cleans up a partial scratch dir when the copy itself errors", {
  ssd_dir <- tempfile("ssd_scratch_"); dir.create(ssd_dir)
  config <- fake_config(extra = list(ssd_scratch_dir = ssd_dir, ssd_min_free_gb = 1))
  file1 <- tempfile("fastq_"); writeLines("x", file1)
  output_dir <- tempfile("bam_")

  testthat::local_mocked_bindings(
    detect_drive = function(...) "/dev/fake", .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Drive_Free = "999G"), .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    file_copy = function(...) stop("simulated ENOSPC: no space left on device"), .package = "fs"
  )

  # Not expect_warning(): in this testthat version it captures and
  # returns the warning CONDITION itself, not resolve_align_scratch()'s
  # own return value -- withCallingHandlers() + muffleWarning gets both.
  warned <- FALSE
  result <- withCallingHandlers(
    resolve_align_scratch(file1, output_dir, "SRR001", config),
    warning = function(w) {
      expect_match(conditionMessage(w), "falling back to main drive")
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )

  expect_true(warned)
  expect_false(result$used_ssd)
  expect_identical(result$output.dir, output_dir)
  expect_length(list.files(ssd_dir, recursive = TRUE), 0) # partial scratch dir was cleaned up
})

test_that("resolve_align_scratch() treats an unresolvable drive (detect_drive error) the same as insufficient space", {
  ssd_dir <- tempfile("ssd_scratch_"); dir.create(ssd_dir)
  config <- fake_config(extra = list(ssd_scratch_dir = ssd_dir, ssd_min_free_gb = 1))
  file1 <- tempfile("fastq_"); writeLines("x", file1)
  output_dir <- tempfile("bam_")

  testthat::local_mocked_bindings(
    detect_drive = function(...) stop("simulated: could not determine drive"), .package = "ORFik"
  )

  result <- resolve_align_scratch(file1, output_dir, "SRR001", config)
  expect_false(result$used_ssd)
})

test_that("pipeline_align_one_organism() retries on the main drive when STAR fails on the SSD scratch copy, and propagates a genuine failure when SSD is disabled", {
  config <- fake_config(preset = "Ribo-seq", extra = list(ssd_scratch_dir = tempfile("ssd_"), ssd_min_free_gb = 1))
  dir.create(config$ssd_scratch_dir)
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(
    bam_dir = bam_dir,
    runs = data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp_name <- "PRJNA000001-homo_sapiens"
  fake_mark_all_done(config, c("trim", "collapsed"), exp_name)
  trimmed_single_dir <- file.path(bam_dir, "trim", "SINGLE")
  dir.create(trimmed_single_dir, recursive = TRUE)
  input_file <- file.path(trimmed_single_dir, "collapsed_trimmed_SRR001.fasta.gz")
  writeLines("x", input_file)

  testthat::local_mocked_bindings(
    run_files_organizer = function(runs, ...) list(input_file)
  )
  testthat::local_mocked_bindings(
    STAR.install = function(...) "star", install.fastp = function(...) "fastp",
    detect_drive = function(...) "/dev/fake",
    get_system_usage = function(...) list(Drive_Free = "999G"),
    .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    run_experiment_subprocess = function(func, args, ...) do.call(func, args)
  )
  testthat::local_mocked_bindings(
    alignment_final_checks = function(...) invisible(NULL)
  )

  calls <- list()
  testthat::local_mocked_bindings(
    STAR.align.single = function(file1, file2, output.dir, ...) {
      calls[[length(calls) + 1]] <<- output.dir
      # Fail only the FIRST call (the SSD-scratch attempt); succeed on retry.
      if (length(calls) == 1) stop("simulated STAR crash on SSD scratch")
      invisible(NULL)
    },
    .package = "ORFik"
  )

  expect_warning(pipeline_align_one_organism(pipelines[["PRJNA000001"]], "Homo sapiens", config),
                 "retrying on main drive")

  expect_length(calls, 2)
  expect_true(startsWith(calls[[1]], config$ssd_scratch_dir)) # first attempt was on SSD scratch
  expect_identical(unname(calls[[2]]), bam_dir) # retry used the real main-drive output dir (output_dir = conf["bam"] carries a "bam" name)
  expect_setequal(samples_done(config, "aligned", exp_name), "SRR001") # sample still ends up marked done
})

test_that("pipeline_align_one_organism() propagates a genuine STAR failure unchanged when SSD scratch is disabled", {
  config <- fake_config(preset = "Ribo-seq") # ssd_scratch_dir = NULL
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(
    bam_dir = bam_dir,
    runs = data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                                  LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens")
  )
  exp_name <- "PRJNA000001-homo_sapiens"
  fake_mark_all_done(config, c("trim", "collapsed"), exp_name)
  trimmed_single_dir <- file.path(bam_dir, "trim", "SINGLE")
  dir.create(trimmed_single_dir, recursive = TRUE)
  input_file <- file.path(trimmed_single_dir, "collapsed_trimmed_SRR001.fasta.gz")
  writeLines("x", input_file)

  testthat::local_mocked_bindings(
    run_files_organizer = function(runs, ...) list(input_file)
  )
  testthat::local_mocked_bindings(
    STAR.install = function(...) "star", install.fastp = function(...) "fastp",
    STAR.align.single = function(...) stop("simulated genuine STAR crash"),
    .package = "ORFik"
  )
  testthat::local_mocked_bindings(
    run_experiment_subprocess = function(func, args, ...) do.call(func, args)
  )

  expect_error(pipeline_align_one_organism(pipelines[["PRJNA000001"]], "Homo sapiens", config),
              "simulated genuine STAR crash")
  expect_setequal(samples_done(config, "aligned", exp_name), character()) # never marked done
})
