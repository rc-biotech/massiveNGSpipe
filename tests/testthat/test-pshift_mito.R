# Mitochondria-aware P-shifting (R/pshift_mito.R)

test_that("is_mito_sample recognises mito FRACTION terms, also when ORFik joins columns", {
  df <- data.frame(fraction = c("mito", "Mitochondria__chl", "mitoribosomes", "MitoticArrest", "cyto", NA))
  expect_identical(is_mito_sample(df), c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(is_mito_sample(data.frame(condition = "WT")), FALSE)
})

test_that("offsets come from the preferred frame, within 11..13 nt", {
  ft <- data.table::data.table(len = c(29L, 33L, 34L, 30L, 31L), reads = c(5000, 5000, 5000, 500, 5000),
                               best_frame = c(0L, 2L, 1L, 0L, 0L), frame_frac = c(0.5, 0.55, 0.5, 0.6, 0.35))
  tab <- mito_offsets_from_frames(ft, min_reads = 1000, min_frame_frac = 0.42)
  # 30 nt: too few reads; 31 nt: no frame preference
  expect_identical(tab$fraction, c(29L, 33L, 34L))
  expect_identical(tab$offsets_start, c(-12L, -13L, -11L))
  # P-site lands on codon position 1: (frame + offset) %% 3 == 0
  expect_true(all((ft$best_frame[1:3] - tab$offsets_start) %% 3L == 0L))
  expect_identical(nrow(mito_offsets_from_frames(ft[0])), 0L)
})

test_that("mito_frame_table finds the frame of 5' ends on a mitochondrial CDS body", {
  cds <- GenomicRanges::GRangesList(tx1 = GenomicRanges::GRanges("MT", IRanges::IRanges(1001, 1600), "+"))
  # 30-nt reads with 5' ends at CDS positions p, (p - 1) %% 3 == 2 (frame 2), plus a few off-frame reads
  p <- seq(33, 500, by = 3)
  starts <- c(1000 + p, 1000 + p[1:20] + 1)
  reads <- GenomicRanges::GRanges("MT", IRanges::IRanges(starts, width = 30), "+")
  S4Vectors::mcols(reads)$score <- c(rep(20L, length(p)), rep(1L, 20))
  ft <- mito_frame_table(reads, cds, lengths = 30)
  expect_identical(ft$best_frame, 2L)
  expect_gt(ft$frame_frac, 0.9)
  expect_equal(ft$reads, 20 * length(p) + 20)
  expect_identical(mito_offsets_from_frames(ft, min_reads = 100)$offsets_start, -13L)
  # reads outside the CDS body / other lengths are ignored
  expect_identical(nrow(mito_frame_table(reads, cds, lengths = 28)), 0L)
})
