# Single-read shift QC helpers (R/shift_qc_single_pass.R): the frame table, the read-length distribution and the
# provenance table index must equal what the ORFik / legacy code returns.

make_qc_fixture <- function() {
  cds <- GenomicRanges::GRangesList(
    txP = GenomicRanges::GRanges("1", IRanges::IRanges(c(100, 200), c(129, 229)), "+"),
    txM = GenomicRanges::GRanges("1", IRanges::IRanges(c(500, 400), c(530, 420)), "-"),
    txQ = GenomicRanges::GRanges("2", IRanges::IRanges(10, 49), "+"))
  set.seed(1)
  pos <- c(sample(c(100:129, 200:229), 60, TRUE), sample(c(400:420, 500:530), 60, TRUE), sample(10:49, 30, TRUE), 5000, 130)
  chr <- c(rep("1", 120), rep("2", 30), "1", "1")
  str <- c(rep("+", 60), rep("-", 60), rep("+", 30), "+", "-")
  reads <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1), str)
  reads$size <- sample(c(26L, 27L, 28L, 30L), length(reads), TRUE)
  reads$score <- sample(1:3, length(reads), TRUE)
  reads
}

test_that("qc_frame_sums() equals ORFik::regionPerReadLength(frameSumPerL) incl. row order", {
  skip_if_not_installed("fst")
  cds <- GenomicRanges::GRangesList(
    txP = GenomicRanges::GRanges("1", IRanges::IRanges(c(100, 200), c(129, 229)), "+"),
    txM = GenomicRanges::GRanges("1", IRanges::IRanges(c(500, 400), c(530, 420)), "-"),
    txQ = GenomicRanges::GRanges("2", IRanges::IRanges(10, 49), "+"))
  reads <- make_qc_fixture()
  f <- tempfile(fileext = ".ofst")
  fst::write_fst(data.frame(seqnames = as.character(GenomicRanges::seqnames(reads)), start = BiocGenerics::start(reads),
                            width = 1L, strand = as.character(BiocGenerics::strand(reads)), size = reads$size,
                            score = reads$score), f)
  dt <- qc_read_ofst(f)
  new <- qc_frame_sums(dt, qc_cds_exon_table(cds))
  old <- ORFik::regionPerReadLength(cds, reads, withFrames = TRUE, scoring = "frameSumPerL", weight = "score",
                                    drop.zero.dt = TRUE, exclude.zero.cov.grl = FALSE,
                                    BPPARAM = BiocParallel::SerialParam())
  expect_identical(new$fraction, old$fraction)
  expect_identical(new$frame, old$frame)
  expect_equal(new$score, old$score)
})

test_that("qc_length_distribution() equals read_length_distribution() and handles an empty library", {
  skip_if_not_installed("fst")
  reads <- make_qc_fixture()
  f <- tempfile(fileext = ".ofst")
  fst::write_fst(data.frame(seqnames = as.character(GenomicRanges::seqnames(reads)), start = BiocGenerics::start(reads),
                            width = 1L, strand = as.character(BiocGenerics::strand(reads)), size = reads$size,
                            score = reads$score), f)
  ref <- read_length_distribution(f)
  expect_equal(qc_length_distribution(qc_read_ofst(f)), ref)
  g <- tempfile(fileext = ".ofst")
  fst::write_fst(data.frame(seqnames = character(), start = integer(), width = integer(), strand = character(),
                            size = integer(), score = integer()), g)
  expect_equal(nrow(qc_length_distribution(qc_read_ofst(g))), 0L)
  expect_equal(nrow(qc_frame_sums(qc_read_ofst(g), qc_cds_exon_table(GenomicRanges::GRangesList(a = GenomicRanges::GRanges("1", IRanges::IRanges(1, 9), "+"))))), 0L)
})

test_that("pshift_table_entry() with a prebuilt index equals the per-call read", {
  raw <- tempfile(fileext = ".ofst"); file.create(raw)
  tab <- list(data.frame(fraction = 28:29, offsets_start = c(-12L, -13L)))
  names(tab) <- raw
  tp <- tempfile(fileext = ".rds"); saveRDS(tab, tp)
  idx <- pshift_table_index(tp)
  expect_identical(pshift_table_entry(raw, tp), pshift_table_entry(raw, tp, index = idx))
  expect_null(pshift_table_entry(tempfile(), tp, index = idx))
  expect_null(pshift_table_index(tempfile()))
})

test_that("remove_sample_flags(runs =) deletes only those samples' markers; default removes the whole folder", {
  config <- fake_config(preset = "Ribo-seq")
  step <- names(config$flag)[2]
  dir <- sample_flag_dir(config, step, "expA")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (r in c("SRR1", "SRR2", "SRR3")) saveRDS(TRUE, file.path(dir, paste0(r, ".rds")))
  remove_sample_flags(config, step, "expA", runs = c("SRR1", "SRR3"))
  expect_identical(list.files(dir), "SRR2.rds")
  remove_sample_flags(config, step, "expA")
  expect_false(dir.exists(dir))
})

test_that("qc_hitmap() equals ORFik::windowPerReadLength() (identical)", {
  skip_if_not_installed("fst")
  cds <- GenomicRanges::GRangesList(
    txP = GenomicRanges::GRanges("1", IRanges::IRanges(c(100, 200), c(129, 229)), "+"),
    txM = GenomicRanges::GRanges("1", IRanges::IRanges(c(500, 400), c(530, 420)), "-"))
  mrna <- GenomicRanges::GRangesList(
    txP = GenomicRanges::GRanges("1", IRanges::IRanges(c(90, 200), c(129, 260)), "+"),
    txM = GenomicRanges::GRanges("1", IRanges::IRanges(c(500, 380), c(540, 420)), "-"))
  set.seed(2)
  pos <- c(sample(c(95:125, 200:215), 80, TRUE), sample(c(405:420, 500:535), 80, TRUE))
  reads <- GenomicRanges::GRanges("1", IRanges::IRanges(pos, width = 1), rep(c("+", "-"), each = 80))
  reads$size <- sample(c(26L, 27L, 28L), length(reads), TRUE); reads$score <- sample(1:3, length(reads), TRUE)
  f <- tempfile(fileext = ".ofst")
  fst::write_fst(data.frame(seqnames = "1", start = BiocGenerics::start(reads), width = 1L,
                            strand = as.character(BiocGenerics::strand(reads)), size = reads$size, score = reads$score), f)
  old <- ORFik::windowPerReadLength(cds, mrna, reads, upstream = 5, downstream = 20)
  wt <- qc_window_table(cds, mrna, 5, 20)
  expect_identical(qc_hitmap(qc_read_ofst(f), wt, 5), old)
})
