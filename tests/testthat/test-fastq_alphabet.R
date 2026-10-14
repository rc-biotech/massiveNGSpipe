test_that("fastq_is_colorspace detects SOLiD colour calls, not base space", {
  w <- function(s) { f <- tempfile(fileext = ".fastq.gz"); con <- gzfile(f, "w")
    writeLines(unlist(lapply(seq_along(s), function(i) c(paste0("@r", i), s[i], "+", strrep("I", nchar(s[i]))))), con); close(con); f }
  cs <- w(c("T0120.3301222", "T3301200112", "G0012303")); bs <- w(c("ACGTTGCA", "TTGACGNNA", "GGGCCCAT"))
  expect_true(fastq_is_colorspace(cs)); expect_false(fastq_is_colorspace(bs))
  expect_false(fastq_is_colorspace(tempfile()))
  e <- tempfile(fileext = ".fastq"); file.create(e); expect_false(fastq_is_colorspace(e))
  h <- tempfile(fileext = ".fastq"); writeLines("@r1", h); expect_false(fastq_is_colorspace(h))
  expect_error(assert_base_space_fastq(c(bs, cs), "SRR1"), "colour-space")
  expect_true(assert_base_space_fastq(bs, "SRR1"))
})
