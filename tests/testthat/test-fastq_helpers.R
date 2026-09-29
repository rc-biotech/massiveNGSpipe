# Unit tests for R/fastq_helpers.R -- the pure/tempdir()-friendly pieces.
# The real STAR/fastqc/adapter-detection subprocess chain
# (detect_adapter_and_trim, run_fastqc, barcode_detector_single/pipeline,
# run_barcode_detection_and_trim) is out of scope here.

test_that("run_files_organizer_internal matches a SINGLE run's plain filename", {
  runs <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE")
  files <- c("/d/SRR001.fastq.gz", "/d/SRR002.fastq.gz")
  result <- run_files_organizer_internal(1, runs, files)
  expect_identical(result[1], "/d/SRR001.fastq.gz")
})

test_that("run_files_organizer_internal matches PAIRED _1/_2 suffixes", {
  runs <- data.table::data.table(Run = "SRR001", LibraryLayout = "PAIRED")
  files <- c("/d/SRR001_1.fastq.gz", "/d/SRR001_2.fastq.gz")
  result <- run_files_organizer_internal(1, runs, files)
  expect_identical(unname(result), c("/d/SRR001_1.fastq.gz", "/d/SRR001_2.fastq.gz"))
})

test_that("run_files_organizer_internal falls back to the second paired_end_suffixes variant", {
  runs <- data.table::data.table(Run = "SRR001", LibraryLayout = "PAIRED")
  files <- c("/d/SRR001_R1_001.fastq.gz", "/d/SRR001_R2_001.fastq.gz")
  result <- run_files_organizer_internal(1, runs, files)
  expect_identical(unname(result), c("/d/SRR001_R1_001.fastq.gz", "/d/SRR001_R2_001.fastq.gz"))
})

test_that("run_files_organizer_internal errors clearly (rather than guessing) on a genuinely ambiguous match", {
  runs <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE")
  files <- c("/d/SRR001.fastq.gz", "/d/trimmed_SRR001.fastq.gz")
  # grep("SRR001", files) matches BOTH; this doesn't resolve to a single
  # exact prefix+format+compression combination, so it errors instead of
  # silently picking one -- verified this is the real, current behavior
  # rather than assumed.
  expect_error(run_files_organizer_internal(1, runs, files), "File format could not be detected")
})

test_that("run_files_organizer_internal errors when a run has no matching file at all", {
  runs <- data.table::data.table(Run = "SRR999", LibraryLayout = "SINGLE")
  files <- c("/d/SRR001.fastq.gz")
  expect_error(run_files_organizer_internal(1, runs, files), "does not exist")
})

test_that("run_files_organizer: end-to-end over a real tempdir with mixed SINGLE/PAIRED runs", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, c("SRR001.fastq.gz", "SRR002_1.fastq.gz", "SRR002_2.fastq.gz")))
  runs <- data.table::data.table(
    Run = c("SRR001", "SRR002"),
    LibraryLayout = c("SINGLE", "PAIRED")
  )
  result <- run_files_organizer(runs, d)
  expect_length(result, 2)
  expect_identical(basename(result[[1]]), "SRR001.fastq.gz")
  expect_setequal(basename(result[[2]]), c("SRR002_1.fastq.gz", "SRR002_2.fastq.gz"))
})

test_that("run_files_organizer errors if source_dir doesn't exist, or has no matching-format files", {
  runs <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE")
  expect_error(run_files_organizer(runs, tempfile()), "does not exist")

  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "not_a_fastq.json"))
  expect_error(run_files_organizer(runs, d), "no files of specified formats")
})

test_that("run_files_organizer warns about extra unmatched files when extra_files_warning=TRUE", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, c("SRR001.fastq.gz", "SRR002.fastq.gz")))
  runs <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE") # only 1 run, 2 files present
  expect_warning(run_files_organizer(runs, d), "additional files")
  expect_no_warning(run_files_organizer(runs, d, extra_files_warning = FALSE))
})

test_that("adapter_list writes and returns the built-in candidate table on first call", {
  f <- tempfile()
  expect_false(file.exists(f))
  candidates <- adapter_list(f)
  expect_true(file.exists(f))
  expect_setequal(colnames(candidates), c("name", "value"))
  expect_true(nrow(candidates) > 0)
})

test_that("adapter_list re-reads an existing candidates file on subsequent calls", {
  f <- tempfile()
  first <- adapter_list(f)
  second <- adapter_list(f)
  expect_identical(first, second)
})

test_that("move_trimmed_files moves the trimmed fastq plus its json/html siblings, stripping the trimmed_ prefix", {
  trimmed_dir <- tempfile(); dir.create(trimmed_dir)
  barcode_dir <- tempfile(); dir.create(barcode_dir)
  writeLines("x", file.path(trimmed_dir, "trimmed_SRR001.fastq.gz"))
  writeLines("{}", file.path(trimmed_dir, "trimmed_SRR001.json"))
  writeLines("<html></html>", file.path(trimmed_dir, "trimmed_SRR001.html"))
  study_sample <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE")

  move_trimmed_files(study_sample, trimmed_dir, barcode_dir)

  expect_setequal(list.files(barcode_dir), c("SRR001.fastq.gz", "SRR001.json", "SRR001.html"))
  expect_length(list.files(trimmed_dir), 0)
})

test_that("barcodes_manual_assign_table writes the expected CSV columns/content", {
  d <- tempfile(); dir.create(d)
  barcodes_manual_assign_table(d, run_ids = c("SRR001", "SRR002"),
                               barcode5p_size = 3, barcode3p_size = 4)
  out <- data.table::fread(file.path(d, "barcodes_manual.csv"))
  expect_identical(out$Run, c("SRR001", "SRR002"))
  expect_true(all(out$barcode5p_size == 3))
  expect_true(all(out$barcode3p_size == 4))
})

test_that("adapters_manual_assign_table writes the expected CSV columns/content", {
  d <- tempfile(); dir.create(d)
  adapters_manual_assign_table(d, run_ids = c("SRR001", "SRR002"), adapters = "AGATCGGAAGAG")
  out <- data.table::fread(file.path(d, "adapters_manual.csv"))
  expect_identical(out$Run, c("SRR001", "SRR002"))
  expect_true(all(out$adapters == "AGATCGGAAGAG"))
})

test_that("trim_flanks_ORFik trims fixed-length 5'/3' flanks and records consensus attributes", {
  reads <- Biostrings::DNAStringSet(c("AAAACCCCGGGG", "AAAATTTTGGGG"))
  trimmed <- trim_flanks_ORFik(reads, left = 4, right = 4)
  expect_identical(as.character(trimmed), c("CCCC", "TTTT"))
  expect_true(!is.null(attr(trimmed, "flank_5p")) || !is.null(attr(trimmed, "flank_left")) ||
             length(attributes(trimmed)) > length(attributes(reads)))
})

test_that("trim_flanks_ORFik: left=0 and right=0 leaves reads untouched", {
  reads <- Biostrings::DNAStringSet(c("AAAACCCCGGGG"))
  trimmed <- trim_flanks_ORFik(reads, left = 0, right = 0)
  expect_identical(as.character(trimmed), as.character(reads))
})

test_that("trim_flanks_ORFik rejects negative left/right", {
  reads <- Biostrings::DNAStringSet(c("AAAACCCCGGGG"))
  expect_error(trim_flanks_ORFik(reads, left = -1, right = 0))
})

test_that("subseqSafe passes through unchanged when start and end are both NA", {
  reads <- Biostrings::DNAStringSet(c("ACGT", "TTTT"))
  expect_identical(as.character(subseqSafe(reads)), as.character(reads))
})

test_that("subseqSafe blanks out-of-range elements when include_invalid=TRUE, drops them when FALSE", {
  reads <- Biostrings::DNAStringSet(c("ACGT", "TT"))
  # end=3 is out of range for the second (2nt) read
  blanked <- subseqSafe(reads, start = 1, end = 3, include_invalid = TRUE)
  expect_length(blanked, 2)
  expect_identical(as.character(blanked[2]), "")

  dropped <- subseqSafe(reads, start = 1, end = 3, include_invalid = FALSE)
  expect_length(dropped, 1)
})

test_that("subseqSafeNonEmptyTrim3p trims each read at its first adapter-hit position", {
  reads <- Biostrings::DNAStringSet(c("ACGTACGT", "TTTTTTTT"))
  hits <- IRanges::IntegerList(list(5L, integer(0))) # hit at pos 5 in read 1, none in read 2
  trimmed <- subseqSafeNonEmptyTrim3p(reads, hits)
  expect_identical(as.character(trimmed[1]), "ACGT")
  expect_identical(as.character(trimmed[2]), "TTTTTTTT") # no hit -> untouched
  expect_identical(attr(trimmed, "empty"), c(FALSE, TRUE))
})

test_that("remove_adapter_ORFik trims reads at the adapter and can attach statistics", {
  reads <- Biostrings::DNAStringSet(c("ACGTAGATCGGAAGAG", "TTTTTTTTTTTTTTTT"))
  invisible(capture.output(suppressMessages(
    trimmed <- remove_adapter_ORFik(reads, adapter = "AGATCGGAAGAG", add_statistics = TRUE)
  )))
  expect_identical(as.character(trimmed[1]), "ACGT")
  expect_true(!is.null(attr(trimmed, "statistics")) || !is.null(attributes(trimmed)$statistics))
})

test_that("remove_adapter_ORFik: add_statistics=FALSE skips the stats attribute", {
  reads <- Biostrings::DNAStringSet(c("ACGTAGATCGGAAGAG"))
  invisible(capture.output(suppressMessages(
    trimmed <- remove_adapter_ORFik(reads, adapter = "AGATCGGAAGAG", add_statistics = FALSE)
  )))
  expect_null(attr(trimmed, "statistics"))
})

test_that("barcode_change_point infers a barcode size on a synthetic step-change curve", {
  skip_if_not_installed("changepoint")
  # curves is a list of "curve type" groups, each a list of individual
  # curve vectors -- NOT a flat numeric vector (a real structural trap:
  # getting this nesting wrong sends single numbers through cpt.mean()
  # one at a time, which errors with "Minimum segment length is too
  # large"). A clear step: low values (barcode region) then a jump up.
  curve <- c(rep(1, 20), rep(20, 30))
  curves <- list(list(curve))
  result <- barcode_change_point(curves, max_barcode_left_size = 15,
                                 max_size_before = 40, max_size_after = 50,
                                 minimum_size = 10, Q = 1)
  expect_named(result, c("barcode5p_size", "barcode3p_size"))
  expect_true(all(result >= 0))
})
