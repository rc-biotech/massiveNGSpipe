# frame_concentration_table (R/pshift_diagnostics.R)

test_that("frame_concentration_table separates broad phase from a single-locus artefact", {
  # 20 genes of 600 nt on + strand; broad phased 30-nt reads (frame 0) on all genes, plus a 26-nt pile at one site
  starts <- seq(1, by = 1000, length.out = 20)
  cds <- GenomicRanges::GRangesList(lapply(starts, function(s) GenomicRanges::GRanges("1", IRanges::IRanges(s, width = 600), "+")))
  names(cds) <- paste0("g", 1:20)
  pos30 <- unlist(lapply(starts, function(s) s - 1 + seq(40, 500, by = 3)))
  r30 <- GenomicRanges::GRanges("1", IRanges::IRanges(pos30, width = 30), "+")
  pile <- GenomicRanges::GRanges("1", IRanges::IRanges(rep(starts[1] - 1 + 101, 1), width = 26), "+")
  noise26 <- GenomicRanges::GRanges("1", IRanges::IRanges(unlist(lapply(starts, function(s) s - 1 + 40:60)), width = 26), "+")
  reads <- c(r30, pile, noise26)
  S4Vectors::mcols(reads)$score <- c(rep(1L, length(r30)), 5000L, rep(1L, length(noise26)))
  ft <- frame_concentration_table(reads, cds, lengths = c(26, 30))
  l30 <- ft[len == 30]; l26 <- ft[len == 26]
  expect_identical(l30$best_frame, 0L); expect_gt(l30$frame_frac, 0.95); expect_false(l30$concentrated)
  expect_equal(l30$genes, 20L)
  # 26 nt: looks phased overall, but one site carries almost all reads and the rest has no preference
  expect_gt(l26$frame_frac, 0.9); expect_true(l26$concentrated)
  expect_gt(l26$top_site_frac, 0.9); expect_lt(l26$frame_frac_without_top, 0.5)
  expect_identical(nrow(frame_concentration_table(reads[0], cds)), 0L)
})
