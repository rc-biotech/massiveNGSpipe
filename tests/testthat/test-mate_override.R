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

test_that("paired_footprint_mate picks read 2 only when read 1 has no insert", {
  dir <- file.path(tempdir(), "footprint_mate"); unlink(dir, recursive = TRUE); dir.create(dir)
  set.seed(2)
  rnd <- function(n, k) vapply(seq_len(n), function(i) paste(sample(c("A", "C", "G", "T"), k, TRUE), collapse = ""), "")
  wfq <- function(s, f) writeLines(paste0("@r", seq_along(s), "\n", s, "\n+\n", strrep("I", nchar(s))), f)
  n <- 2000; ad <- "AGATCGGAAGAGCGTCGTGTAGGG"
  junk_r1 <- paste0("ACGTCTGAACTCCAGTCACGATC", strrep("G", 40))
  good_r1 <- paste0(rnd(n, 30), "AGATCGGAAGAGCACACGTCTG")
  good_r2 <- paste0("GAT", rnd(n, 5), rnd(n, 30), "TC", ad)
  wfq(rep(junk_r1, n), file.path(dir, "junk_1.fastq")); wfq(good_r2, file.path(dir, "junk_2.fastq"))
  wfq(good_r1, file.path(dir, "ok_1.fastq")); wfq(good_r2, file.path(dir, "ok_2.fastq"))
  d1 <- massiveNGSpipe:::paired_footprint_mate(file.path(dir, "junk_1.fastq"), file.path(dir, "junk_2.fastq"))
  expect_equal(d1$mate, 2L); expect_true(d1$reverse)
  d2 <- massiveNGSpipe:::paired_footprint_mate(file.path(dir, "ok_1.fastq"), file.path(dir, "ok_2.fastq"))
  expect_equal(d2$mate, 1L); expect_false(d2$reverse)                           # usable read 1: default path unchanged
  wfq(rep(junk_r1, n), file.path(dir, "bad_2.fastq"))
  d3 <- massiveNGSpipe:::paired_footprint_mate(file.path(dir, "junk_1.fastq"), file.path(dir, "bad_2.fastq"))
  expect_equal(d3$mate, 1L)                                                     # read 2 without inserts: no switch
})

test_that("paired footprint runs are trimmed and collapsed from one mate only, as single-end (read 1 or read 2)", {
  skip_if(inherits(try(massiveNGSpipe:::trimcollapse_binary(), silent = TRUE), "try-error"), "trimcollapse binary not available")
  skip_if(!nzchar(Sys.which("fastp")), "fastp not available")
  set.seed(3)
  rnd <- function(n, k) vapply(seq_len(n), function(i) paste(sample(c("A", "C", "G", "T"), k, TRUE), collapse = ""), "")
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  wfq <- function(s, f) writeLines(paste0("@r", seq_along(s), "\n", s, "\n+\n", strrep("I", nchar(s))), f)
  fp <- rep(rnd(300, 30), each = 3)                                                # each footprint 3x -> collapses
  run_case <- function(name, r1, r2) {
    dir <- file.path(tempdir(), paste0("footprint_mate_", name)); unlink(dir, recursive = TRUE)
    raw <- file.path(dir, "raw"); trimmed_dir <- file.path(dir, "trim"); dir.create(raw, recursive = TRUE); dir.create(trimmed_dir)
    f1 <- file.path(raw, paste0("SRR", name, "_1.fastq")); f2 <- file.path(raw, paste0("SRR", name, "_2.fastq"))
    wfq(r1, f1); wfq(r2, f2)
    d <- massiveNGSpipe:::paired_footprint_mate(f1, f2)
    ss <- data.table::data.table(Run = paste0("SRR", name), LibraryLayout = "PAIRED", LIBRARYTYPE = "RFP")
    bt <- massiveNGSpipe:::trim_from_footprint_mate(ss, if (d$mate == 2L) f2 else f1, dir, trimmed_dir, "online", d)
    list(d = d, bt = bt, trimmed_dir = trimmed_dir, run = ss$Run,
         out = file.path(trimmed_dir, "SINGLE", paste0("collapsed_trimmed_", ss$Run, ".fasta.gz")))
  }
  # standard design: read 1 = footprint + linker + UMI + TruSeq, read 2 = UMI + reverse complement + junk
  r1 <- paste0(fp, "CTGTAGGCACCATCAAT", rnd(length(fp), 10), "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC")
  r2 <- paste0(rnd(length(fp), 10), "ATTGATGGTGCCTACAG", rc(fp), "AGATCGGAAGAGCGTCGTGTAGGG")
  a <- run_case("MATE1", r1, r2)
  expect_equal(a$d$mate, 1L); expect_false(a$d$reverse)
  expect_true(file.exists(a$out))
  s <- Biostrings::readDNAStringSet(a$out)
  expect_gt(mean(as.character(s) %in% unique(fp)), 0.95)                           # sense footprints from read 1
  auto <- data.table::fread(file.path(a$trimmed_dir, "mate_override_auto.csv"), colClasses = "character")
  expect_equal(auto$Run, a$run); expect_equal(auto$mate, "1"); expect_equal(auto$reverse, "FALSE")
  runs <- data.table::data.table(Run = a$run, LibraryLayout = "PAIRED")
  expect_equal(massiveNGSpipe:::apply_mate_override_layout(runs, a$trimmed_dir)$LibraryLayout, "SINGLE")
  # read 1 without insert (adapter / index read-through + poly-G): read 2, reverse complemented. Fixed flanks are
  # detected automatically; random UMIs inside read 2 are not (use mate_override_manual.csv trim_front for those)
  junk1 <- rep(paste0("ACGTCTGAACTCCAGTCACGATC", strrep("G", 40)), length(fp))
  r2b <- paste0("GAT", rc(fp), "TC", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
  b <- run_case("MATE2", junk1, r2b)
  expect_equal(b$d$mate, 2L); expect_true(b$d$reverse)
  s2 <- Biostrings::readDNAStringSet(b$out)
  expect_gt(mean(as.character(s2) %in% unique(fp)), 0.95)                          # back in the RNA sense
})
