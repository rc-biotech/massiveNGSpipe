test_that("orient_collapsed_fasta: opt-in, reverse complement keeps headers/weights/lengths, idempotent, re-collapse flips again", {
  d <- tempfile("trim"); dir.create(file.path(d, "SINGLE"), recursive = TRUE)
  f <- file.path(d, "SINGLE", "collapsed_trimmed_R1.fasta.gz")
  orig <- Biostrings::DNAStringSet(c(seq1_x10 = "ACGTTT", seq2_x3 = "GGGCAA"))
  Biostrings::writeXStringSet(orig, f, compress = TRUE)
  expect_identical(orient_collapsed_fasta(f, "R1", d), "not_marked")                 # default off
  fwrite(data.table(Run = c("R1", "R2"), reverse = c("TRUE", "FALSE"), source = "test", note = ""), file.path(d, "reads_reverse_manual.csv"))
  expect_identical(reads_reverse_runs(d), "R1")
  expect_identical(orient_collapsed_fasta(f, "R1", d), "flipped")
  x <- Biostrings::readDNAStringSet(f)
  expect_identical(names(x), names(orig)); expect_identical(width(x), width(orig))      # headers = weights, lengths
  expect_identical(as.character(x), c(seq1_x10 = "AAACGT", seq2_x3 = "TTGCCC"))
  expect_identical(orient_collapsed_fasta(f, "R1", d), "already_flipped")              # no double flip on resume
  expect_identical(as.character(Biostrings::readDNAStringSet(f)), as.character(x))
  Biostrings::writeXStringSet(orig, f, compress = TRUE)                               # re-collapse writes native again
  expect_identical(orient_collapsed_fasta(f, "R1", d), "flipped")
  expect_identical(as.character(Biostrings::reverseComplement(Biostrings::readDNAStringSet(f))), as.character(orig))  # involution
  expect_error(orient_collapsed_fasta(f, "R1", d, paired = TRUE), "paired-end")
  expect_false(any(grepl("orient_tmp", list.files(file.path(d, "SINGLE"), all.files = TRUE))))   # no temp leftovers
  expect_identical(nrow(fread(file.path(d, ".orientation_manifest.csv"))), 2L)
})
