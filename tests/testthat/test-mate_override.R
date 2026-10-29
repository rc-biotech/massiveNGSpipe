test_that("mate override trims and collapses read 2, reverse complemented, named by Run", {
  skip_if(inherits(try(massiveNGSpipe:::trimcollapse_binary(), silent = TRUE), "try-error"), "trimcollapse binary not available")
  dir <- file.path(tempdir(), "mate_override_test"); unlink(dir, recursive = TRUE)
  raw <- file.path(dir, "raw"); trimmed_dir <- file.path(dir, "trim"); dir.create(raw, recursive = TRUE); dir.create(trimmed_dir)
  set.seed(1)
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  fp <- vapply(1:200, function(i) paste(sample(c("A", "C", "G", "T"), 30, TRUE), collapse = ""), "")
  fp <- rep(fp, each = 2)                                                           # every footprint twice -> collapses
  umi <- vapply(seq_along(fp), function(i) paste(sample(c("A", "C", "G", "T"), 5, TRUE), collapse = ""), "")
  r2 <- paste0("GAT", umi, rc(fp), "TC", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
  r1 <- rep(paste0("ACGTCTGAACTCCAGTCACGATC", strrep("G", 40)), length(fp))
  wfq <- function(s, f) writeLines(paste0("@r", seq_along(s), "\n", s, "\n+\n", strrep("I", nchar(s))), f)
  wfq(r1, file.path(raw, "SRRTEST_1.fastq")); wfq(r2, file.path(raw, "SRRTEST_2.fastq"))
  data.table::fwrite(data.table::data.table(Run = "SRRTEST", mate = 2, reverse = TRUE, adapter = "AGATCGGAAGAGCGT",
                                            trim_front = 8, trim_tail = 2, source = "test"), file.path(trimmed_dir, "mate_override_manual.csv"))
  expect_equal(massiveNGSpipe:::mate_override_runs(trimmed_dir), "SRRTEST")
  runs <- data.table::data.table(Run = c("SRRTEST", "SRROTHER"), LibraryLayout = "PAIRED")
  expect_equal(massiveNGSpipe:::apply_mate_override_layout(runs, trimmed_dir)$LibraryLayout, c("SINGLE", "PAIRED"))
  expect_equal(runs$LibraryLayout, c("PAIRED", "PAIRED"))                             # input table not modified
  bt <- massiveNGSpipe:::trim_with_mate_override("SRRTEST", file.path(raw, c("SRRTEST_1.fastq", "SRRTEST_2.fastq")), dir, trimmed_dir)
  out <- file.path(trimmed_dir, "SINGLE", "collapsed_trimmed_SRRTEST.fasta.gz")
  expect_true(file.exists(out)); expect_true(file.exists(file.path(trimmed_dir, "report_SRRTEST.json")))
  expect_true(file.exists(file.path(trimmed_dir, ".mate_override", "SRRTEST.json")))
  s <- Biostrings::readDNAStringSet(out)
  expect_setequal(as.character(s), unique(fp))                                       # sense footprints recovered exactly
  expect_true(all(grepl("_x2$", names(s))))                                          # duplicates collapsed with count 2
  expect_equal(bt$mate_override, 2L)
})

test_that("mate override table is parsed strictly", {
  d <- file.path(tempdir(), "mate_override_bad"); unlink(d, recursive = TRUE); dir.create(d)
  expect_equal(nrow(massiveNGSpipe:::mate_overrides(d)), 0)                         # no table: nothing changes
  data.table::fwrite(data.table::data.table(Run = "A", mate = 3, reverse = TRUE, adapter = "AGATC", trim_front = 0,
                                            trim_tail = 0, source = "x"), file.path(d, "mate_override_manual.csv"))
  expect_error(massiveNGSpipe:::mate_overrides(d), "mate must be 1 or 2")
})
