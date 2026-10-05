# Unit tests for R/download_sra.R -- the pure/tempdir()-friendly pieces
# only (existence/validity checks, filename reconciliation, size
# estimation). The real network/subprocess download chain
# (download_sra/download_raw_srr/download_sra_aws/ascp/ebi/sra_to_fastq)
# is out of scope: needs real network access or real sratoolkit/aws/ascp
# binaries, not a "quick unit test."

test_that("fastq_output_exists_and_valid: SINGLE layout, compressed", {
  d <- tempfile(); dir.create(d)
  expect_false(fastq_output_exists_and_valid("SRR001", d, FALSE, TRUE))
  writeLines("x", file.path(d, "SRR001.fastq.gz"))
  expect_true(fastq_output_exists_and_valid("SRR001", d, FALSE, TRUE))
})

test_that("fastq_output_exists_and_valid: PAIRED layout needs both files, non-empty", {
  d <- tempfile(); dir.create(d)
  writeLines("x", file.path(d, "SRR001_1.fastq.gz"))
  expect_false(fastq_output_exists_and_valid("SRR001", d, TRUE, TRUE)) # only _1 present
  writeLines("x", file.path(d, "SRR001_2.fastq.gz"))
  expect_true(fastq_output_exists_and_valid("SRR001", d, TRUE, TRUE))
})

test_that("fastq_output_exists_and_valid: a zero-byte file doesn't count as valid", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001.fastq.gz"))
  expect_identical(file.info(file.path(d, "SRR001.fastq.gz"))$size, 0)
  expect_false(fastq_output_exists_and_valid("SRR001", d, FALSE, TRUE))
})

test_that("fastq_output_exists_and_valid: compress=FALSE expects the uncompressed suffix", {
  d <- tempfile(); dir.create(d)
  writeLines("x", file.path(d, "SRR001.fastq"))
  expect_true(fastq_output_exists_and_valid("SRR001", d, FALSE, FALSE))
  expect_false(fastq_output_exists_and_valid("SRR001", d, FALSE, TRUE)) # wrong suffix
})

test_that("sra_or_direct_fastq_format defaults to TRUE (attempt fast path) when info is missing/ambiguous", {
  expect_true(sra_or_direct_fastq_format(NULL, "SRR001"))
  expect_true(sra_or_direct_fastq_format(data.table::data.table(Run = "SRR001"), "SRR001")) # no size_MB col
  dup <- data.table::data.table(Run = c("SRR001", "SRR001"), size_MB = c(50, 200))
  expect_true(sra_or_direct_fastq_format(dup, "SRR001")) # matches 2 rows, not exactly 1
})

test_that("sra_or_direct_fastq_format: uses the >100MB threshold (strict >) when info is unambiguous", {
  big <- data.table::data.table(Run = "SRR001", size_MB = 150)
  small <- data.table::data.table(Run = "SRR001", size_MB = 50)
  boundary <- data.table::data.table(Run = "SRR001", size_MB = 100)
  expect_true(sra_or_direct_fastq_format(big, "SRR001"))
  expect_false(sra_or_direct_fastq_format(small, "SRR001"))
  expect_false(sra_or_direct_fastq_format(boundary, "SRR001")) # strict >, not >=
})

test_that("sra_or_direct_fastq_format: NA size_MB falls back to TRUE", {
  na_size <- data.table::data.table(Run = "SRR001", size_MB = NA_real_)
  expect_true(sra_or_direct_fastq_format(na_size, "SRR001"))
})

test_that("delete_existing_preformat_files removes only the matching preformat files", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001"), file.path(d, "SRR002"), file.path(d, "SRR001.fastq.gz"))
  delete_existing_preformat_files(d, c("SRR001"))
  expect_false(file.exists(file.path(d, "SRR001")))
  expect_true(file.exists(file.path(d, "SRR002"))) # not in accessions, untouched
  expect_true(file.exists(file.path(d, "SRR001.fastq.gz"))) # different file, untouched
})

test_that("delete_existing_preformat_files: delete_srr_preformat=FALSE is a no-op", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001"))
  delete_existing_preformat_files(d, "SRR001", delete_srr_preformat = FALSE)
  expect_true(file.exists(file.path(d, "SRR001")))
})

test_that("delete_existing_preformat_files: none of the accessions present is a silent no-op", {
  d <- tempfile(); dir.create(d)
  expect_no_error(delete_existing_preformat_files(d, c("SRR999")))
})

test_that("validate_fastq_download: single fastq present, PAIRED_END=FALSE, matches expectation", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001.fastq"))
  result <- validate_fastq_download("SRR001.fastq", "SRR001", PAIRED_END = FALSE,
                                    is_compressed = FALSE, outdir = d)
  expect_identical(result, "SRR001.fastq")
})

test_that("validate_fastq_download: renames a lone <accession>_1.fastq to the plain name for SINGLE-end", {
  # Regression test: the old fastq-dump fallback (used when both the
  # AWS/SRA-toolkit path and the EBI fallback fail -- verified live for
  # SRR27790697, PRJNA1071171-homo_sapiens, in a network-restricted
  # environment) names its SINGLE-end output <accession>_1.fastq, not the
  # plain name every other path expects. Before this fix, a fully and
  # correctly downloaded file was deleted and the whole download treated
  # as failed purely over this naming difference.
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001_1.fastq"))
  result <- validate_fastq_download("SRR001_1.fastq", "SRR001", PAIRED_END = FALSE,
                                    is_compressed = FALSE, outdir = d)
  expect_identical(result, "SRR001.fastq")
  expect_true(file.exists(file.path(d, "SRR001.fastq")))
  expect_false(file.exists(file.path(d, "SRR001_1.fastq")))
})

test_that("validate_fastq_download: a lone <accession>_1.fastq still fails when PAIRED_END=TRUE (missing _2)", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001_1.fastq"))
  expect_error(
    validate_fastq_download("SRR001_1.fastq", "SRR001", PAIRED_END = TRUE,
                            is_compressed = FALSE, outdir = d),
    "invalid filenames"
  )
})

test_that("validate_fastq_download: single fastq present but PAIRED_END=TRUE fails", {
  d <- tempfile(); dir.create(d)
  expect_error(
    validate_fastq_download("SRR001.fastq", "SRR001", PAIRED_END = TRUE,
                            is_compressed = FALSE, outdir = d)
  )
})

test_that("validate_fastq_download: two fastqs present, PAIRED_END=TRUE, matches expectation", {
  d <- tempfile(); dir.create(d)
  files <- c("SRR001_1.fastq", "SRR001_2.fastq")
  result <- validate_fastq_download(files, "SRR001", PAIRED_END = TRUE,
                                    is_compressed = FALSE, outdir = d)
  expect_setequal(result, files)
})

test_that("validate_fastq_download: deletes a wrong-compression duplicate with a warning", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001.fastq"))
  expect_warning(
    result <- validate_fastq_download(c("SRR001.fastq.gz", "SRR001.fastq"), "SRR001",
                                      PAIRED_END = FALSE, is_compressed = TRUE, outdir = d),
    "Non compressed version"
  )
  expect_false(file.exists(file.path(d, "SRR001.fastq"))) # the uncompressed dup got deleted
  expect_identical(result, "SRR001.fastq.gz")
})

test_that("validate_fastq_download: totally unrecognized filenames get deleted and it stops", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001_weird_name.fastq"))
  expect_error(
    validate_fastq_download("SRR001_weird_name.fastq", "SRR001", PAIRED_END = FALSE,
                            is_compressed = FALSE, outdir = d),
    "invalid filenames"
  )
  expect_false(file.exists(file.path(d, "SRR001_weird_name.fastq")))
})

test_that("cleanup_and_validate_fastq_download errors when nothing matches the accession", {
  d <- tempfile(); dir.create(d)
  expect_error(
    cleanup_and_validate_fastq_download("SRR001", d, PAIRED_END = FALSE, compress = FALSE),
    "no files with pattern"
  )
})

test_that("cleanup_and_validate_fastq_download errors when files match accession but not the suffix", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001.json"))
  expect_error(
    cleanup_and_validate_fastq_download("SRR001", d, PAIRED_END = FALSE, compress = FALSE),
    "no file with valid suffix"
  )
})

test_that("cleanup_and_validate_fastq_download: happy path for a single-end run", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "SRR001.fastq"))
  result <- cleanup_and_validate_fastq_download("SRR001", d, PAIRED_END = FALSE, compress = FALSE)
  expect_identical(result, "SRR001.fastq")
})

test_that("estimate_fastq_tmp_gb: falls back to bases/spots when avgLength is NA", {
  row <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                                spots = 1000, bases = 100000, avgLength = NA_real_)
  est <- estimate_fastq_tmp_gb(row)
  expect_identical(est$read_length, 100L) # 100000/1000
  expect_true(est$total_bytes > 0)
})

test_that("estimate_fastq_tmp_gb: PAIRED doubles the read count relative to spots", {
  se <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                               spots = 1000, bases = 100000, avgLength = 100)
  pe <- data.table::data.table(Run = "SRR001", LibraryLayout = "PAIRED",
                               spots = 1000, bases = 100000, avgLength = 100)
  expect_identical(estimate_fastq_tmp_gb(pe)$reads, 2 * estimate_fastq_tmp_gb(se)$reads)
})

test_that("estimate_fastq_tmp_gb: read_len override takes precedence over avgLength", {
  row <- data.table::data.table(Run = "SRR001", LibraryLayout = "SINGLE",
                                spots = 100, bases = 10000, avgLength = 100)
  est <- estimate_fastq_tmp_gb(row, read_len = 50)
  expect_identical(est$read_length, 50L)
})
