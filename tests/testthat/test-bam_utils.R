# Unit tests for R/bam_utils.R's true (uncollapsed) alignment-metrics
# orchestration: save_expanded_alignment_metrics() and
# expanded_alignment_metrics_plot(). get_expanded_alignment_metrics()
# itself (real BAM/fasta parsing) is mocked out here -- these tests cover
# massiveNGSpipe's own orchestration (SINGLE/PAIRED filtering, missing-
# file handling, never-fails-the-caller error containment), not
# Rsamtools/GenomicAlignments parsing correctness.

fake_metrics_row <- function(input = 1000, aligned = 950, unique = 900, multi = 50) {
  data.table::data.table(
    total_input_seqs = input, total_input_seqs_collapsed = input,
    total_aligned_reads = aligned, total_aligned_reads_collapsed = aligned,
    total_alignments = aligned, total_alignments_collapsed = aligned,
    unique_aligned_reads = unique, unique_aligned_reads_collapsed = unique,
    multimapping_aligned_reads = multi, multimapping_aligned_reads_collapsed = multi,
    `total_input_seqs(million reads)` = input / 1e6,
    `total_aligned_reads%` = round(100 * aligned / input, 2),
    `unique_alignments%` = round(100 * unique / input, 2),
    `multimapping_aligned_reads%` = round(100 * multi / input, 2)
  )
}

test_that("get_expanded_alignment_metrics() reproduces the user's own worked example, with no duplicate columns", {
  # The motivating scenario, verbatim: one collapsed sequence represents
  # 1,000,000 original reads and maps; a second, singleton sequence (1
  # original read) fails to map. Collapsed-level reporting would show a
  # misleading 50% (1 of 2 unique sequences mapped); the true rate is
  # ~99.9999%. get_qname() (the only real-BAM-touching part) is mocked
  # out -- everything else here is the function's own real dt_summary
  # construction logic, which is what's actually being verified.
  fasta_file <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(big_x1000000 = "ACGTACGTACGT", singleton_x1 = "TTTTGGGGCCCC")),
    fasta_file, format = "fasta")

  testthat::local_mocked_bindings(get_qname = function(bam, ...) "big_x1000000")

  res <- get_expanded_alignment_metrics(fasta_file, "fake.bam")

  expect_identical(anyDuplicated(colnames(res)), 0L)
  expect_equal(res$total_input_seqs, 1000001)
  expect_equal(res$total_aligned_reads, 1000000)
  expect_equal(res$unique_aligned_reads, 1000000)
  expect_equal(res$unique_aligned_reads_collapsed, 1)
  expect_equal(res[["total_aligned_reads%"]], 100) # rounds to 100%, not the misleading 50%
})

test_that("save_expanded_alignment_metrics() writes a CSV + PNG for SINGLE-end samples", {
  bam_dir <- tempfile("bam_"); dir.create(bam_dir)
  collapsed_dir <- tempfile("collapsed_"); dir.create(collapsed_dir)
  out_dir <- tempfile("qc_")

  study_org <- data.table::data.table(Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE")
  for (run in study_org$Run) {
    file.create(file.path(bam_dir, paste0(run, ".bam")))
    # Real naming convention written by pipeline_collapse() -- see the
    # regression test below for why this prefix matters.
    file.create(file.path(collapsed_dir, paste0("collapsed_trimmed_", run, ".fasta.gz")))
  }

  testthat::local_mocked_bindings(
    get_expanded_alignment_metrics = function(fasta_file, bam_file) fake_metrics_row()
  )

  save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org, out_dir = out_dir)

  csv_path <- file.path(out_dir, "expanded_alignment_metrics.csv")
  expect_true(file.exists(csv_path))
  dt <- data.table::fread(csv_path)
  expect_identical(dt$Run, c("SRR001", "SRR002"))
  expect_true(file.exists(file.path(out_dir, "expanded_alignment_metrics.png")))
})

test_that("save_expanded_alignment_metrics() skips PAIRED samples with a message, still processes SINGLE", {
  bam_dir <- tempfile("bam_"); dir.create(bam_dir)
  collapsed_dir <- tempfile("collapsed_"); dir.create(collapsed_dir)
  out_dir <- tempfile("qc_")

  study_org <- data.table::data.table(Run = c("SRR001", "SRR002"),
                                      LibraryLayout = c("SINGLE", "PAIRED"))
  for (run in study_org$Run) {
    file.create(file.path(bam_dir, paste0(run, ".bam")))
    file.create(file.path(collapsed_dir, paste0("collapsed_trimmed_", run, ".fasta.gz")))
  }
  testthat::local_mocked_bindings(
    get_expanded_alignment_metrics = function(fasta_file, bam_file) fake_metrics_row()
  )

  expect_message(
    save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org, out_dir = out_dir),
    "Skipping expanded alignment metrics for 1 PAIRED-end"
  )

  dt <- data.table::fread(file.path(out_dir, "expanded_alignment_metrics.csv"))
  expect_identical(dt$Run, "SRR001")
})

test_that("save_expanded_alignment_metrics() skips samples missing a bam or collapsed fasta", {
  bam_dir <- tempfile("bam_"); dir.create(bam_dir)
  collapsed_dir <- tempfile("collapsed_"); dir.create(collapsed_dir)
  out_dir <- tempfile("qc_")

  study_org <- data.table::data.table(Run = c("SRR001", "SRR002"), LibraryLayout = "SINGLE")
  # Only SRR001 gets both files; SRR002 is missing its collapsed fasta.
  file.create(file.path(bam_dir, "SRR001.bam"))
  file.create(file.path(collapsed_dir, "collapsed_trimmed_SRR001.fasta.gz"))
  file.create(file.path(bam_dir, "SRR002.bam"))

  testthat::local_mocked_bindings(
    get_expanded_alignment_metrics = function(fasta_file, bam_file) fake_metrics_row()
  )

  expect_message(
    save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org, out_dir = out_dir),
    "Missing bam or collapsed fasta for 1 sample"
  )
  dt <- data.table::fread(file.path(out_dir, "expanded_alignment_metrics.csv"))
  expect_identical(dt$Run, "SRR001")
})

test_that("save_expanded_alignment_metrics() falls back to a plain .fasta when .fasta.gz doesn't exist", {
  bam_dir <- tempfile("bam_"); dir.create(bam_dir)
  collapsed_dir <- tempfile("collapsed_"); dir.create(collapsed_dir)
  out_dir <- tempfile("qc_")

  study_org <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE")
  file.create(file.path(bam_dir, "SRR001.bam"))
  file.create(file.path(collapsed_dir, "collapsed_trimmed_SRR001.fasta")) # no .gz

  seen_fasta <- NULL
  testthat::local_mocked_bindings(
    get_expanded_alignment_metrics = function(fasta_file, bam_file) {
      seen_fasta <<- fasta_file
      fake_metrics_row()
    }
  )
  save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org, out_dir = out_dir)
  expect_identical(seen_fasta, file.path(collapsed_dir, "collapsed_trimmed_SRR001.fasta"))
})

test_that("save_expanded_alignment_metrics() finds the real collapsed_trimmed_<Run> naming, not a plain <Run>.fasta.gz guess", {
  # Regression test: pipeline_collapse() actually writes
  # collapsed_trimmed_<Run>.fasta.gz (verified live against production
  # data), not <Run>.fasta.gz -- a hardcoded paste0(Run, ".fasta.gz")
  # guess silently matches nothing for every real study. Also covers the
  # plain "<Run>.fasta.gz" (no prefix) and "trimmed_<Run>.fasta.gz"
  # variants run_files_organizer() supports, so any of the three
  # naming eras this pipeline has used resolve correctly.
  bam_dir <- tempfile("bam_"); dir.create(bam_dir)
  collapsed_dir <- tempfile("collapsed_"); dir.create(collapsed_dir)
  out_dir <- tempfile("qc_")

  study_org <- data.table::data.table(Run = c("SRR001", "SRR002", "SRR003"),
                                      LibraryLayout = "SINGLE")
  for (run in study_org$Run) file.create(file.path(bam_dir, paste0(run, ".bam")))
  file.create(file.path(collapsed_dir, "collapsed_trimmed_SRR001.fasta.gz"))
  file.create(file.path(collapsed_dir, "trimmed_SRR002.fasta.gz"))
  file.create(file.path(collapsed_dir, "SRR003.fasta.gz"))

  testthat::local_mocked_bindings(
    get_expanded_alignment_metrics = function(fasta_file, bam_file) fake_metrics_row()
  )

  save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org, out_dir = out_dir)

  dt <- data.table::fread(file.path(out_dir, "expanded_alignment_metrics.csv"))
  expect_setequal(dt$Run, c("SRR001", "SRR002", "SRR003"))
})

test_that("save_expanded_alignment_metrics() does nothing (no error) when there are no SINGLE-end samples", {
  bam_dir <- tempfile("bam_"); dir.create(bam_dir)
  collapsed_dir <- tempfile("collapsed_"); dir.create(collapsed_dir)
  out_dir <- tempfile("qc_")
  study_org <- data.table::data.table(Run = "SRR001", LibraryLayout = "PAIRED")

  expect_message(
    save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org, out_dir = out_dir),
    "No SINGLE-end samples"
  )
  expect_false(dir.exists(out_dir))
})

test_that("save_expanded_alignment_metrics() never throws, even if the underlying computation errors", {
  # This is the core safety property: a QC-only failure here must not
  # propagate up and permanently skip the whole study for the rest of
  # the session via pipe_align_clean()'s own outer try()/report_failed_pipe().
  bam_dir <- tempfile("bam_"); dir.create(bam_dir)
  collapsed_dir <- tempfile("collapsed_"); dir.create(collapsed_dir)
  out_dir <- tempfile("qc_")
  study_org <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE")
  file.create(file.path(bam_dir, "SRR001.bam"))
  file.create(file.path(collapsed_dir, "collapsed_trimmed_SRR001.fasta.gz"))

  testthat::local_mocked_bindings(
    get_expanded_alignment_metrics = function(fasta_file, bam_file) stop("boom: corrupt bam")
  )

  expect_warning(
    expect_no_error(save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org, out_dir = out_dir)),
    "QC-only, not fatal"
  )
})

test_that("expanded_alignment_metrics_plot() returns a ggplot built from the unique/multi/unmapped percentages", {
  dt <- data.table::rbindlist(list(
    cbind(Run = "SRR001", fake_metrics_row(input = 1000, aligned = 950, unique = 900, multi = 50)),
    cbind(Run = "SRR002", fake_metrics_row(input = 2000, aligned = 100, unique = 60, multi = 40))
  ))
  p <- expanded_alignment_metrics_plot(dt)
  expect_s3_class(p, "ggplot")
  expect_identical(sort(unique(as.character(p$data$Run))), c("SRR001", "SRR002"))
  expect_setequal(levels(p$data$category), c("Unmapped", "Multimapping", "Unique"))
  # SRR002: 100 - total_aligned_reads% (95%) = 5% unmapped for SRR001;
  # for SRR002 total_aligned_reads% = 5%, so unmapped = 95%.
  unmapped_srr002 <- p$data$percent[p$data$Run == "SRR002" & p$data$category == "Unmapped"]
  expect_equal(unmapped_srr002, 95)
})
