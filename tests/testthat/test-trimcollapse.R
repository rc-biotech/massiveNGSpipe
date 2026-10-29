library(data.table)

# Synthetic single-end reads: [UMI5 2][insert][UMI3+barcode 19][TruSeq adapter][polyG], 150 nt, NextSeq-style header
make_reads <- function(file, n = 3000, seed = 1) {
  set.seed(seed)
  rnd <- function(k) paste(sample(c("A", "C", "G", "T"), k, replace = TRUE), collapse = "")
  adapter <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  lines <- character(4 * n)
  for (i in seq_len(n)) {
    seq <- substr(paste0(rnd(2), rnd(sample(5:45, 1)), rnd(19), adapter, strrep("G", 150)), 1, 150)
    if (i %% 17 == 0) seq <- rnd(150)                                   # no adapter
    if (i %% 29 == 0) seq <- substr(paste0("AGATCGTGAT", rnd(140)), 1, 150)  # adapter-like start (fastp indel fallback)
    q <- paste(sample(c("F", ":", ","), 150, replace = TRUE, prob = c(0.85, 0.1, 0.05)), collapse = "")
    lines[(4 * i - 3):(4 * i)] <- c(sprintf("@VH00270:1:FC:1:1101:%d:1000 1:N:0:AAAA", i), seq, "+", q)
  }
  con <- gzfile(file, "w"); writeLines(lines, con); close(con)
  file
}
read_fq <- function(f) { x <- readLines(f); setNames(x[seq(2, length(x), 4)], sub(" .*", "", x[seq(1, length(x), 4)])) }

test_that("trimcollapse builds and writes the pipeline's file names and locations", {
  skip_if(Sys.which("cc") == "", "no C compiler")
  d <- file.path(tempdir(), "tc_names"); dir.create(d, showWarnings = FALSE)
  f <- make_reads(file.path(d, "SRR000001.fastq.gz"))
  out <- trimcollapse_run(f, file.path(d, "proc"), "AGATCGGAAGAG", threads = 2)
  expect_equal(out, file.path(d, "proc", "trim", "SINGLE", "collapsed_trimmed_SRR000001.fasta.gz"))
  js <- jsonlite::fromJSON(file.path(d, "proc", "trim", "report_SRR000001.json"))
  expect_equal(js$summary$before_filtering$total_reads, 3000)
  expect_equal(js$adapter_cutting$read1_adapter_sequence, "AGATCGGAAGAG")
  tt <- ORFik::trimming.table(file.path(d, "proc", "trim"))                 # ORFik's own reader of fastp reports
  expect_equal(tt$raw_reads, 3000)
  fa <- Biostrings::readDNAStringSet(out)
  expect_equal(sum(as.integer(sub(".*_x", "", names(fa)))), js$summary$after_filtering$total_reads)
})

test_that("fast path is used only where it reproduces fastp", {
  skip_if(Sys.which("cc") == "", "no C compiler")
  d <- file.path(tempdir(), "tc_rules"); dir.create(d, showWarnings = FALSE)
  fq <- file.path(d, "a.fastq.gz"); fq2 <- file.path(d, "a.fq"); fa <- file.path(d, "a.fasta.gz")
  for (f in c(fq, fq2, fa)) file.create(f)
  expect_true(fast_trim_enabled(fq, NULL, "AGATCGGAAGAG"))
  expect_true(fast_trim_enabled(fq2, NA, "disable"))
  expect_false(fast_trim_enabled(fq, file.path(d, "a_2.fastq.gz"), "AGATCGGAAGAG"))   # paired end
  expect_false(fast_trim_enabled(fq, NULL, "auto"))                                  # fastp auto detection
  expect_false(fast_trim_enabled(fa, NULL, "AGATCGGAAGAG"))                          # fasta input
  expect_false(fast_trim_enabled(file.path(d, "missing.fastq.gz"), NULL, "AGATCGGAAGAG"))
  withr::with_options(list(massiveNGSpipe.fast_trim = FALSE), expect_false(fast_trim_enabled(fq, NULL, "AGATCGGAAGAG")))
})

fastp_bin <- path.expand("~/bin/fastp")
test_that("trimmed reads and report equal fastp 1.3.7 (main pass, post-adapter tail, detection sample)", {
  skip_if(Sys.which("cc") == "" || !file.exists(fastp_bin), "needs a C compiler and ~/bin/fastp")
  skip_if(!grepl("1\\.3\\.7", paste(system2(fastp_bin, "--version", stdout = TRUE, stderr = TRUE), collapse = " ")), "fastp 1.3.7 required")
  d <- file.path(tempdir(), "tc_fastp"); dir.create(d, showWarnings = FALSE)
  f <- make_reads(file.path(d, "S1.fastq.gz"), n = 4000, seed = 3)
  cmp <- function(fp_fq, fp_js, tc_fq, tc_js) {
    expect_identical(read_fq(tc_fq)[names(read_fq(fp_fq))], read_fq(fp_fq))
    a <- jsonlite::fromJSON(fp_js); b <- jsonlite::fromJSON(tc_js)
    expect_equal(b$filtering_result[names(a$filtering_result)], a$filtering_result)   # same values; key order may differ
    expect_equal(b$summary$after_filtering, a$summary$after_filtering, tolerance = 1e-4)
    expect_equal(b$adapter_cutting$adapter_trimmed_reads, a$adapter_cutting$adapter_trimmed_reads)
    expect_equal(b$read1_after_filtering$content_curves, a$read1_after_filtering$content_curves, tolerance = 1e-3)
  }
  # main pass as ORFik::STAR.align.single(steps = "tr")
  system2(fastp_bin, c("--in1", f, "--out1", file.path(d, "fp.fq"), "--json", file.path(d, "fp.json"), "--html", file.path(d, "fp.html"),
                       "--trim_front1", 2, "--trim_tail1", 0, "--length_required", 20, "--thread", 4, "--disable_quality_filtering",
                       "--adapter_sequence", "AGATCGGAAGAG"), stdout = FALSE, stderr = FALSE)
  trimcollapse_run(f, d, "AGATCGGAAGAG", trim.front = 2, fastq_out = file.path(d, "tc.fq"), json = file.path(d, "tc.json"), threads = 2)
  cmp(file.path(d, "fp.fq"), file.path(d, "fp.json"), file.path(d, "tc.fq"), file.path(d, "tc.json"))
  # 3' barcode/UMI after the adapter, as fastp_post_adapter_tail()
  fastp_post_adapter_tail(f, NULL, output.dir = file.path(d, "pat"), adapter = "AGATCGGAAGAG", trim.front = 2, trim.tail = 19, threads = 2)
  trimcollapse_run(f, d, "AGATCGGAAGAG", trim.front = 2, tail_after_adapter = 19, fastq_out = file.path(d, "tc_pat.fq"),
                   json = file.path(d, "tc_pat.json"), threads = 2)
  cmp(file.path(d, "pat", "trim", "trimmed_S1.fastq"), file.path(d, "pat", "trim", "report_S1.json"), file.path(d, "tc_pat.fq"), file.path(d, "tc_pat.json"))
  # barcode detection sample (fastp_detect_sample_args: quality filter on)
  args <- fastp_detect_sample_args(f, file.path(d, "fs.fq"), file.path(d, "fs.json"), file.path(d, "fs.html"), "AGATCGGAAGAG", 2000)
  system2(fastp_bin, args, stdout = FALSE, stderr = FALSE)
  trimcollapse_run(f, d, "AGATCGGAAGAG", quality_filtering = TRUE, reads_to_process = 2000, fastq_out = file.path(d, "ts.fq"),
                   json = file.path(d, "ts.json"), threads = 2)
  cmp(file.path(d, "fs.fq"), file.path(d, "fs.json"), file.path(d, "ts.fq"), file.path(d, "ts.json"))
})

test_that("fast_collapse_fastq equals ORFik::collapse.fastq (sequences, counts, names, file name)", {
  skip_if(inherits(try(massiveNGSpipe:::trimcollapse_binary(), silent = TRUE), "try-error"), "trimcollapse binary not available")
  d <- file.path(tempdir(), "fast_collapse_test"); unlink(d, recursive = TRUE); dir.create(file.path(d, "a"), recursive = TRUE); dir.create(file.path(d, "b"))
  set.seed(3)
  s <- vapply(sample(1:300, 5000, TRUE), function(i) { set.seed(i); paste(sample(c("A", "C", "G", "T"), sample(18:40, 1), TRUE), collapse = "") }, "")
  f <- file.path(d, "trimmed_SRRX.fastq")
  writeLines(paste0("@r", seq_along(s), "\n", s, "\n+\n", strrep("I", nchar(s))), f)
  ORFik::collapse.fastq(f, file.path(d, "a"), compress = TRUE)
  out <- massiveNGSpipe:::fast_collapse_fastq(f, file.path(d, "b"), threads = 2)
  a <- list.files(file.path(d, "a"), full.names = TRUE)
  expect_equal(basename(out), basename(a))
  ra <- Biostrings::readDNAStringSet(a); rb <- Biostrings::readDNAStringSet(out)
  o <- order(as.character(ra)); p <- order(as.character(rb))
  expect_identical(unname(as.character(ra))[o], unname(as.character(rb))[p])          # seq<i> index may differ among tied counts
  expect_identical(sub(".*_x", "", names(ra))[o], sub(".*_x", "", names(rb))[p])
})

test_that("adapter_scan_reads picks the insert/adapter junction (known, poly-A tail, linker, barcodes, de novo, none)", {
  d <- file.path(tempdir(), "adapter_scan"); unlink(d, recursive = TRUE); dir.create(d)
  set.seed(5)
  rnd <- function(n, k) vapply(seq_len(n), function(i) paste(sample(c("A", "C", "G", "T"), k, TRUE), collapse = ""), "")
  wfq <- function(s, f) { writeLines(paste0("@r", seq_along(s), "\n", s, "\n+\n", strrep("I", nchar(s))), f); f }
  n <- 3000; ins <- rnd(n, 30); univ <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  scan <- function(s, nm) unname(massiveNGSpipe:::adapter_scan_reads(wfq(substr(s, 1, 75), file.path(d, nm))))
  expect_equal(scan(paste0(ins, univ, strrep("G", 20)), "u.fq"), "AGATCGGAAGAG")                       # standard
  expect_equal(scan(paste0(ins, strrep("A", 15), univ), "pa.fq"), "AAAAAAAAAA")                        # poly(A) tail
  expect_equal(scan(paste0(ins, strrep("A", 14), "TGGTCA", univ), "pal.fq"), "AAAAAAAAAA")             # tail + linker
  expect_equal(scan(paste0(ins, "TCTCCTTGCATAATCACCAACC", univ), "rl.fq"), "TCTCCTTGCATA")              # RiboLace linker
  bc <- sample(c("ATCGT", "CTAGA", "CGTAA", "GATCA", "GCATA"), n, TRUE)
  expect_equal(scan(paste0(ins, bc, univ), "bc.fq"), "AGATCGGAAGAG")                                   # 5 nt barcodes
  expect_equal(scan(paste0(ins, "GTTCAGAGTTCTACAGTCCGACGATC"), "novo.fq"), "GTTCAGAGTTCT")             # not in the list
  expect_equal(scan(rnd(n, 60), "none.fq"), "disable")                                                  # no adapter
})

test_that("collapse survives a dominant sequence (shard buffer growth, Valen PRJNA563539 crash 2026-10-29)", {
  bin <- tryCatch(massiveNGSpipe:::trimcollapse_binary(), error = function(e) NULL)
  skip_if(is.null(bin), "trimcollapse binary not available")
  dir <- file.path(tempdir(), "tc_dominant"); unlink(dir, recursive = TRUE); dir.create(dir)
  set.seed(4); n <- 400000; dom <- "ACGTACGTACGTACGTACGTACGTACGTAC"
  rnd <- vapply(seq_len(n * 0.4), function(i) paste(sample(c("A", "C", "G", "T"), 30, TRUE), collapse = ""), "")
  s <- sample(c(rep(dom, n - length(rnd)), rnd))
  f <- file.path(dir, "dom.fastq"); writeLines(paste0("@r", seq_along(s), "\n", s, "\n+\n", strrep("I", 30)), f)
  out <- massiveNGSpipe:::fast_collapse_fastq(f, dir, threads = 4)
  x <- Biostrings::readDNAStringSet(out)
  counts <- as.integer(sub(".*_x", "", names(x)))
  expect_equal(sum(counts), n)
  expect_equal(counts[as.character(x) == dom], n - length(rnd))
})
