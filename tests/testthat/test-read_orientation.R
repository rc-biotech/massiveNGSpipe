mk <- function(d, run, seqs = c(seq1_x10 = "ACGTTT", seq2_x3 = "GGGCAA")) {
  f <- file.path(d, "SINGLE", paste0("collapsed_trimmed_", run, ".fasta.gz")); dir.create(dirname(f), FALSE, TRUE)
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), f, compress = TRUE); f
}
marks <- function(d, runs, rev = "TRUE") fwrite(data.table(Run = runs, reverse = rev, source = "test", note = ""), file.path(d, "reads_reverse_manual.csv"))

test_that("orient_collapsed_fasta: opt-in, keeps headers/weights/lengths, idempotent, re-collapse flips again", {
  d <- tempfile("trim"); f <- mk(d, "R1")
  expect_identical(orient_collapsed_fasta(f, "R1", d), "not_marked")
  marks(d, c("R1", "R2"), c("TRUE", "FALSE")); expect_identical(reads_reverse_runs(d), "R1")
  expect_identical(orient_collapsed_fasta(f, "R1", d), "flipped")
  x <- Biostrings::readDNAStringSet(f)
  expect_identical(as.character(x), c(seq1_x10 = "AAACGT", seq2_x3 = "TTGCCC"))
  expect_identical(orient_collapsed_fasta(f, "R1", d), "already_flipped")
  f <- mk(d, "R1")                                                                   # re-collapse: native again
  expect_identical(orient_collapsed_fasta(f, "R1", d), "flipped")
  expect_identical(as.character(Biostrings::reverseComplement(Biostrings::readDNAStringSet(f))), c(seq1_x10 = "ACGTTT", seq2_x3 = "GGGCAA"))
  expect_error(orient_collapsed_fasta(f, "R1", d, paired = TRUE), "paired-end")
  expect_false(any(grepl("orient_tmp", list.files(file.path(d, "SINGLE"), all.files = TRUE))))
})

test_that("crash after replace but before the final record: retry repairs state, no double flip", {
  d <- tempfile("trim"); f <- mk(d, "R1"); marks(d, "R1")
  orient_collapsed_fasta(f, "R1", d)
  sp <- orientation_state_path(d, "R1"); st <- jsonlite::read_json(sp); st$state <- "prepared"
  jsonlite::write_json(st, sp, auto_unbox = TRUE)                                    # simulate crash before 'done'
  expect_identical(orient_collapsed_fasta(f, "R1", d), "already_flipped")
  expect_identical(jsonlite::read_json(sp)$state, "done")
  expect_identical(as.character(Biostrings::readDNAStringSet(f))[[1]], "AAACGT")
})

test_that("crash before replace (prepared, file still native): retry flips once", {
  d <- tempfile("trim"); f <- mk(d, "R1"); marks(d, "R1")
  write_json_atomic(list(Run = "R1", input_md5 = file_md5(f), output_md5 = "planned-but-not-written", state = "prepared"), orientation_state_path(d, "R1"))
  expect_identical(orient_collapsed_fasta(f, "R1", d), "flipped")
  expect_identical(as.character(Biostrings::readDNAStringSet(f))[[1]], "AAACGT")
})

test_that("failed rename leaves the original untouched and no temp file", {
  d <- tempfile("trim"); f <- mk(d, "R1"); marks(d, "R1")
  real_rename <- base::file.rename
  testthat::local_mocked_bindings(file.rename = function(from, to) if (grepl("orient_tmp", from)) FALSE else real_rename(from, to), .package = "base")
  expect_error(orient_collapsed_fasta(f, "R1", d), "cannot replace")
  expect_identical(as.character(Biostrings::readDNAStringSet(f))[[1]], "ACGTTT")
  expect_false(any(grepl("orient_tmp", list.files(file.path(d, "SINGLE"), all.files = TRUE))))
})

test_that("concurrent runs keep separate state", {
  d <- tempfile("trim"); fs <- vapply(paste0("R", 1:4), function(r) mk(d, r), character(1)); marks(d, paste0("R", 1:4))
  res <- parallel::mclapply(paste0("R", 1:4), function(r) orient_collapsed_fasta(fs[[r]], r, d), mc.cores = 4)
  expect_true(all(unlist(res) == "flipped"))
  s <- orientation_summary(d); expect_identical(sort(s$Run), paste0("R", 1:4)); expect_true(all(s$state == "done"))
  expect_true(all(unlist(parallel::mclapply(paste0("R", 1:4), function(r) orient_collapsed_fasta(fs[[r]], r, d), mc.cores = 4)) == "already_flipped"))
})

test_that("manual table: malformed bool and duplicated runs are rejected", {
  d <- tempfile("trim"); dir.create(d)
  fwrite(data.table(Run = "R1", reverse = "yes"), file.path(d, "reads_reverse_manual.csv")); expect_error(reads_reverse_runs(d), "TRUE or FALSE")
  fwrite(data.table(Run = c("R1", "R1"), reverse = c("TRUE", "FALSE")), file.path(d, "reads_reverse_manual.csv")); expect_error(reads_reverse_runs(d), "duplicated")
})
