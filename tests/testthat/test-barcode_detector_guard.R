# Guarded barcode detector (2026-11-23): z-score retry calls need support, too-large calls are rejected instead of
# clipped, a length guard for Ribo-seq and a poly(A) adapter fallback. options(massiveNGSpipe.barcode_detector =
# "legacy") keeps the detector of 1914c90 (2026-09-29) .. this change.

rnd_seqs <- function(n, k) {
  if (n == 0) return(character())
  s <- paste(sample(c("A", "C", "G", "T"), n * k, replace = TRUE), collapse = "")
  substring(s, seq(1, n * k, k), seq(k, n * k, k))
}

rnd_var <- function(lengths) {   # random sequences of the given lengths
  s <- rnd_seqs(1, sum(lengths))
  ends <- cumsum(lengths)
  substring(s, ends - lengths + 1, ends)
}

write_fastq_lines <- function(seqs, file) {
  n <- length(seqs)
  writeLines(c(rbind(paste0("@r", seq_len(n)), seqs, "+", strrep("I", nchar(seqs)))), file)
  file
}

illumina <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAT"

# Run barcode_detector_single() the way the pipeline does (cheap detection sample, no download) on one raw fastq.
run_detector <- function(raw_file, run, mode = "guarded", libtype = "RFP") {
  d <- tempfile("bcdet_"); dir.create(file.path(d, "fq"), recursive = TRUE)
  on.exit(unlink(d, recursive = TRUE))
  ext <- if (grepl("\\.gz$", raw_file)) ".fastq.gz" else ".fastq"
  file.copy(raw_file, file.path(d, "fq", paste0(run, ext)))
  ss <- data.table::data.table(Run = run, LibraryLayout = "SINGLE", LIBRARYTYPE = libtype)
  out <- NULL
  withr::with_options(list(massiveNGSpipe.barcode_detector = mode), {
    invisible(capture.output(suppressMessages(suppressWarnings(
      out <- barcode_detector_single(ss, file.path(d, "fq"), d, file.path(d, "trim"),
                                     redownload_raw_if_needed = FALSE, detect_reads_to_process = 1e6)))))
  })
  out
}

skip_without_trim <- function() {
  skip_if(inherits(try(massiveNGSpipe:::trimcollapse_binary(), silent = TRUE), "try-error"),
          "trimcollapse binary not available")
}

fixture <- function(run) testthat::test_path("fixtures", "barcode_regression", paste0(run, ".fastq.gz"))

# --- helpers --------------------------------------------------------------------------------------------------------

test_that("barcode_detector_mode defaults to guarded, accepts legacy and rejects anything else", {
  withr::with_options(list(massiveNGSpipe.barcode_detector = NULL), expect_identical(barcode_detector_mode(), "guarded"))
  withr::with_options(list(massiveNGSpipe.barcode_detector = "legacy"), expect_identical(barcode_detector_mode(), "legacy"))
  withr::with_options(list(massiveNGSpipe.barcode_detector = "new"), expect_error(barcode_detector_mode(), "guarded"))
})

test_that("constrain_barcode_sizes rejects (does not clip) a call that leaves fewer than minimum_size nt", {
  # SRR5766466: 2 + 11 nt on 34 nt reads; legacy clipped the 3' call to 6
  expect_identical(unname(suppressMessages(constrain_barcode_sizes(2, 11, 34, 26))), c(2, 0))
  expect_identical(unname(constrain_barcode_sizes(2, 11, 34, 26, clip = TRUE)), c(2, 6))
  # 5' alone still too large: both rejected
  expect_identical(unname(suppressMessages(constrain_barcode_sizes(10, 3, 30, 26))), c(0, 0))
  # within limits: untouched in both modes
  expect_identical(unname(constrain_barcode_sizes(4, 4, 36, 26)), c(4, 4))
  expect_identical(unname(constrain_barcode_sizes(4, 4, 36, 26, clip = TRUE)), c(4, 4))
})

read_end_fixture <- function(run, side)
  readLines(testthat::test_path("fixtures", "barcode_regression", paste0(run, "_", side, ".txt.gz")))
read_fixture_reads <- function(run)
  readLines(testthat::test_path("fixtures", "barcode_regression", paste0(run, "_reads.txt.gz")))
end_call <- function(seqs, side, ...) {
  x <- suppressMessages(end_barcode_call(seqs, side, ...))
  c(fixed = x[["fixed"]], umi = x[["umi"]])
}

test_that("end_templating_profile: kappa near 0 at random positions, near 1 at templated ones", {
  set.seed(31)
  fp <- rnd_var(sample(28:34, 3000, replace = TRUE))
  ins <- sample(fp, 20000, replace = TRUE)
  ins <- substr(ins, sample(0:2, 20000, replace = TRUE) + 1, 40)            # ragged 5' ends: distinct reads
  p <- end_templating_profile(paste0(rnd_seqs(20000, 3), ins), max_k = 8)
  expect_true(all(p$kappa[3] < 0.2))                                          # innermost random position
  expect_true(all(p$kappa[4:8] > 0.9))                                        # insert
  q <- end_templating_profile(paste0(rnd_seqs(20000, 3), ins), max_k = 3, anchor_start = 4L)
  expect_true(all(q$kappa < 0.2))                                             # outer positions against the insert anchor
  expect_true(all(is.na(end_templating_profile(character(), 5)$kappa)))
})

test_that("end_barcode_call on real reads (12k-read samples of 1M-read ENA prefixes): calls match the soft-clip truth", {
  # GSE73553: master cut 5' 3 (fastp quality curve) + 3' 4 (clipped z-score call); nothing at either end
  expect_equal(end_call(read_end_fixture("SRR2535267", "5p"), "5p"), c(fixed = 0L, umi = 0L))
  expect_equal(end_call(read_end_fixture("SRR2535267", "3p"), "3p"), c(fixed = 0L, umi = 0L))
  # McGlincy & Ingolia 2017: 5' NN from the RT primer; 3' fixed AGATC + 5 nt UMI (master: 5' 5 from a curve tie)
  expect_equal(end_call(read_end_fixture("SRR38141702", "5p"), "5p"), c(fixed = 0L, umi = 2L))
  expect_equal(end_call(read_end_fixture("SRR38141702", "3p"), "3p"), c(fixed = 5L, umi = 5L))
  expect_equal(end_call(read_end_fixture("SRR16039918", "3p"), "3p"), c(fixed = 5L, umi = 5L))
  expect_equal(end_call(read_end_fixture("SRR16039918", "5p"), "5p"), c(fixed = 0L, umi = 2L))
  # RiboLace: random 5 / 4 nt (master 3': 17 clipped to 11)
  expect_equal(end_call(read_end_fixture("SRR34280150", "5p"), "5p"), c(fixed = 0L, umi = 5L))
  expect_equal(end_call(read_end_fixture("SRR34280150", "3p"), "3p"), c(fixed = 0L, umi = 4L))
  # PRJNA1301971 read 2: GAT + 5 nt UMI
  expect_equal(end_call(read_end_fixture("SRR34874063", "5p"), "5p"), c(fixed = 3L, umi = 5L))
  # SRR28976196: 4 random 5' nt (position 5 random in ~30% of reads only: reported, not cut), 6 random 3' nt
  expect_message(end_barcode_call(read_end_fixture("SRR28976196", "5p"), "5p", sample = "t"), "part of the reads")
  expect_equal(end_call(read_end_fixture("SRR28976196", "5p"), "5p"), c(fixed = 0L, umi = 4L))
  expect_equal(end_call(read_end_fixture("SRR28976196", "3p"), "3p"), c(fixed = 0L, umi = 6L))
  # Random flanks without many PCR copies, missed by the duplication (merge) test of cbbd61a:
  # OTTR 5' NNNNN + Y + C (7 nt), 3' none
  ottr <- read_fixture_reads("SRR23242339")
  expect_equal(sum(end_call(ottr, "5p")), 7L)
  expect_equal(sum(end_call(ottr, "3p")), 0L)
  # 7 random nt at both ends (TruSeq-based adapters with N7): outer positions have inflated kappa against their own
  # (random) anchor and are re-tested against the insert anchor
  n7 <- read_fixture_reads("ERR4467418")
  expect_equal(end_call(n7, "5p"), c(fixed = 0L, umi = 7L))
  expect_equal(end_call(n7, "3p"), c(fixed = 0L, umi = 7L))
  # NEXTflex 4N + 4N; and a library described as Ingolia 2012 whose reads carry 4N + 4N as well
  expect_equal(sum(end_call(read_fixture_reads("SRR24083061"), "5p")), 4L)
  expect_equal(sum(end_call(read_fixture_reads("SRR24083061"), "3p")), 4L)
  expect_equal(sum(end_call(read_fixture_reads("SRR12407385"), "5p")), 4L)
  expect_equal(sum(end_call(read_fixture_reads("SRR12407385"), "3p")), 4L)
  # one non-templated 5' base (RealSeq-AC), nothing at 3'
  expect_equal(sum(end_call(read_fixture_reads("SRR22990340"), "5p")), 1L)
  expect_equal(sum(end_call(read_fixture_reads("SRR22990340"), "3p")), 0L)
})

test_that("end_barcode_call finds random flanks in shallow samples (3,000 reads)", {
  set.seed(5)
  s <- read_end_fixture("SRR34280150", "5p"); s <- s[sample(length(s), 3000)]
  expect_equal(sum(end_call(s, "5p")), 5L)
  s <- read_fixture_reads("SRR23242339"); s <- s[sample(length(s), 3000)]
  expect_equal(sum(end_call(s, "5p")), 7L)
})

test_that("end_barcode_call: synthetic fixed + random, deduplicated random flanks, tails, partial additions, few reads", {
  set.seed(21)
  fp <- rnd_var(sample(28:32, 4000, replace = TRUE))
  copies <- sample(fp, 40000, replace = TRUE)
  s <- paste0(copies, rnd_seqs(40000, 4), "CTAGA")                     # insert + 4 nt UMI + fixed CTAGA
  x <- suppressMessages(end_barcode_call(s, "3p"))
  expect_equal(c(x[["fixed"]], x[["umi"]]), c(5L, 4L))
  expect_length(attr(x, "review"), 0)
  # one read per insert (no copy of any molecule; ragged footprint ends) with 3 random nt on each end: still found
  ins <- unique(substr(copies, sample(0:3, 40000, replace = TRUE) + 1, nchar(copies) - sample(0:3, 40000, replace = TRUE)))
  d <- paste0(rnd_seqs(length(ins), 3), ins, rnd_seqs(length(ins), 3))
  expect_equal(end_call(d, "5p"), c(fixed = 0L, umi = 3L))
  expect_equal(end_call(d, "3p"), c(fixed = 0L, umi = 3L))
  # a random base added to 30% of the reads of a real library without flanks: not cut (soft clipping handles it)
  real <- read_end_fixture("SRR2535267", "5p")
  part <- paste0(ifelse(runif(length(real)) < 0.3, rnd_seqs(length(real), 1), ""), real)
  expect_equal(sum(end_call(part, "5p")), 0L)
  # ... and to 50% of the reads (kappa 0.15-0.3: either the single-position guard or the partial band keeps it uncut)
  part <- paste0(ifelse(runif(length(real)) < 0.5, rnd_seqs(length(real), 1), ""), real)
  expect_message(p <- end_barcode_call(part, "5p", sample = "t"), "MANUAL REVIEW")
  expect_equal(sum(p), 0L)
  expect_equal(sum(suppressMessages(end_barcode_call(part, "5p", kappa_cut = 0.35))), 0L)
  # a 6 nt A run is a tail, not a barcode: no cut and a review line
  expect_message(y <- end_barcode_call(paste0(copies, "AAAAAA"), "3p", sample = "t"), "MANUAL REVIEW")
  expect_equal(sum(y), 0L)
  # no flank: nothing
  expect_equal(sum(end_call(ins, "5p")), 0L)
  expect_equal(sum(end_call(ins, "3p")), 0L)
  # random reads (nothing testable): only a fixed run is cut, with a review line
  w <- suppressMessages(end_barcode_call(paste0("CC", rnd_var(sample(36:40, 3000, replace = TRUE))), "5p"))
  expect_equal(c(w[["fixed"]], w[["umi"]]), c(2L, 0L))
  expect_match(attr(w, "review"), "too few read pairs")
  # untrimmed reads (full cycle length, no molecule 3' end) are left out of the 3' test
  long <- substr(paste0(sample(fp, 30000, replace = TRUE), strrep("G", 20)), 1, 36)
  mix <- c(paste0(sample(ins, 20000, replace = TRUE), rnd_seqs(20000, 4)), long[1:8000])
  expect_equal(sum(end_call(unique(mix), "3p", full_length = 36)), 4L)
})

test_that("post_trim_length_guard rejects cuts that move the footprint mode out of the accepted lengths", {
  acc <- c(20, 21, 25:33)
  reads <- vapply(rep(c(28, 29, 29, 30), 50), function(l) strrep("A", l), "")
  bs <- c(barcode5p_size = 0, barcode3p_size = 6)
  expect_warning(out <- post_trim_length_guard(reads, bs, acc, "SRRtest"), "outside accepted_lengths_rpf")
  expect_identical(unname(out), c(0, 0))
  # reject = FALSE (guarded detector): warn, keep the cut
  expect_warning(kept <- post_trim_length_guard(reads, bs, acc, "SRRtest", reject = FALSE), "MANUAL REVIEW")
  expect_identical(kept, bs)
  # a cut that keeps the mode inside is kept
  expect_identical(post_trim_length_guard(reads, c(barcode5p_size = 2, barcode3p_size = 1), acc), c(barcode5p_size = 2, barcode3p_size = 1))
  # untrimmed mode outside (real 4 + 4 UMIs around 29 nt footprints): the lengths cannot judge, call kept
  umi_reads <- vapply(rep(c(36, 37, 37, 38), 50), function(l) strrep("A", l), "")
  expect_identical(post_trim_length_guard(umi_reads, c(barcode5p_size = 4, barcode3p_size = 4), acc),
                   c(barcode5p_size = 4, barcode3p_size = 4))
  expect_identical(post_trim_length_guard(reads, c(barcode5p_size = 0, barcode3p_size = 0), acc),
                   c(barcode5p_size = 0, barcode3p_size = 0))
})

test_that("is_ribo_library: RFP or no LIBRARYTYPE column", {
  expect_true(is_ribo_library(data.table::data.table(Run = "a")))
  expect_true(is_ribo_library(data.table::data.table(Run = "a", LIBRARYTYPE = "RFP")))
  expect_false(is_ribo_library(data.table::data.table(Run = "a", LIBRARYTYPE = "disome")))
})

test_that("polyA_tail_fraction counts footprint + A tails, not pure A or adapter read-through", {
  set.seed(2)
  fp <- paste0(rnd_seqs(700, 26), "C")
  tails <- substr(paste0(fp, strrep("A", 20)), 1, 36)
  pureA <- rep(strrep("A", 36), 200)
  p7 <- rep(paste0("ATCTCGTATGCCGTCTTCTGCTTG", strrep("A", 12)), 200)   # flow-cell read-through (PRJNA530200)
  plain <- rnd_seqs(1000, 30)
  expect_equal(polyA_tail_fraction(c(tails, plain)), 700 / 1700)
  expect_equal(polyA_tail_fraction(c(pureA, p7, plain)), 0)
  expect_equal(polyA_tail_fraction(character()), 0)
  # too short an insert before the A run does not count
  expect_equal(polyA_tail_fraction(paste0(rnd_seqs(10, 10), "C", strrep("A", 20))), 0)
})

# --- whole detector on synthetic libraries ---------------------------------------------------------------------------

test_that("synthetic 5' 4 nt random UMI + 3' fixed 4 nt barcode are found by the guarded detector", {
  # (random UMIs at both ends hide each other from the duplication test: no copies become identical)
  skip_without_trim()
  set.seed(11)
  fp <- rnd_var(sample(28:32, 3000, replace = TRUE))
  copies <- sample(fp, 30000, replace = TRUE)
  reads <- substr(paste0(rnd_seqs(30000, 4), copies, "TGCA", illumina), 1, 60)
  f <- write_fastq_lines(reads, tempfile(fileext = ".fastq"))
  new <- run_detector(f, "SRR0000001")
  expect_equal(c(new$barcode5p_size, new$barcode3p_size), c(4, 4), ignore_attr = TRUE)
})

test_that("synthetic library without barcode: no cut", {
  skip_without_trim()
  set.seed(12)
  ins <- rnd_var(sample(28:36, 30000, replace = TRUE))
  reads <- substr(paste0(ins, illumina), 1, 60)
  f <- write_fastq_lines(reads, tempfile(fileext = ".fastq"))
  new <- run_detector(f, "SRR0000002")
  expect_equal(c(new$barcode5p_size, new$barcode3p_size), c(0, 0), ignore_attr = TRUE)
})

test_that("synthetic poly(A)-tailed library: poly(A) adapter, no cut", {
  skip_without_trim()
  set.seed(13)
  ins <- paste0(rnd_var(sample(25:29, 30000, replace = TRUE)), "G")
  reads <- substr(paste0(ins, strrep("A", 30)), 1, 35)
  f <- write_fastq_lines(reads, tempfile(fileext = ".fastq"))
  new <- run_detector(f, "SRR0000003")
  expect_identical(new$adapter, adapter_list()[name == "polyA", value])
  expect_equal(c(new$barcode5p_size, new$barcode3p_size), c(0, 0), ignore_attr = TRUE)
  old <- run_detector(f, "SRR0000003", mode = "legacy")
  expect_identical(old$adapter, "passed")
})

# --- regression: 20,000-read prefixes of the libraries of 2026-11-22 -------------------------------------------------

test_that("REGRESSION SRR5766466 (GSE100626 yeast): no spurious 3' cut", {
  skip_without_trim()
  set.seed(1)
  new <- run_detector(fixture("SRR5766466"), "SRR5766466")
  expect_identical(new$adapter, "TGGAATTCTCGG")
  expect_equal(new$barcode3p_size, 0, ignore_attr = TRUE)
  expect_equal(new$barcode5p_size, 2, ignore_attr = TRUE)     # fixed CC, as before the regression
  set.seed(1)
  old <- run_detector(fixture("SRR5766466"), "SRR5766466", mode = "legacy")
  expect_gt(old$barcode3p_size, 0)                             # the regression: 6 nt on the pipeline sample
})

test_that("REGRESSION SRR449467 (GSE36892 poly(A)-tailed): poly(A) adapter instead of a fixed 3' cut", {
  skip_without_trim()
  set.seed(1)
  new <- run_detector(fixture("SRR449467"), "SRR449467")
  expect_identical(new$adapter, adapter_list()[name == "polyA", value])
  expect_equal(c(new$barcode5p_size, new$barcode3p_size), c(0, 0), ignore_attr = TRUE)
  set.seed(1)
  old <- run_detector(fixture("SRR449467"), "SRR449467", mode = "legacy")
  expect_identical(old$adapter, "passed")
  expect_gt(old$barcode3p_size, 0)
})
