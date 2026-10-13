test_that("umap_from_counts is reproducible, thread independent and validates input", {
  skip_if_not_installed("uwot")
  set.seed(1)
  counts <- matrix(rpois(300 * 60, rep(c(2, 20), each = 300 * 30)), 300, 60,
                   dimnames = list(paste0("g", 1:300), paste0("s", 1:60)))
  a <- umap_from_counts(counts, ntop = 100, seed = 1, threads = 1)
  b <- umap_from_counts(counts, ntop = 100, seed = 1, threads = 4)
  expect_identical(names(a), c("UMAP 1", "UMAP 2"))
  expect_equal(nrow(a), 60)
  expect_identical(a, b)
  # chunked log-normalisation/variance gives the same result as one chunk
  expect_identical(a, umap_from_counts(counts, ntop = 100, seed = 1, threads = 1, chunk_rows = 7))
  counts[, 3] <- 0L
  expect_error(umap_from_counts(counts, ntop = 100), "zero counts")
})

test_that("umap_sample_metadata matches the curated metadata layout", {
  meta <- data.table(Run = c("R1", "R2", "R3"), AUTHOR = c("Makar", "", "Smith"),
                     CELL_LINE = c("HeLa", "", NA), TISSUE = c("Cervix", "Liver", ""),
                     INHIBITOR = c("CHX", "", "harr"), BioProject = c("P1", "P2", "P3"))
  m <- umap_sample_metadata(c("R3", "R1", "R2"), meta)
  expect_identical(m$sample, c("R3", "R1", "R2"))
  expect_identical(m$author, c("Smith", "Unknown", "Unknown"))
  expect_identical(m$cell_line, c("unknown", "hela", "unknown"))
  expect_identical(m$tissues_cell_lines, c("unknown | unknown", "cervix | hela", "liver | unknown"))
  expect_identical(m$inhibitors, c("harr", "chx", "unknown"))
  expect_identical(m$studies, m$BioProject)
  expect_error(umap_sample_metadata("nope", meta), "No samples")
})
