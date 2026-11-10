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
  six <- umap_from_counts(counts, ntop = 100, seed = 1, threads = 1, n_components = 6)
  expect_identical(names(six), paste("UMAP", 1:6))
  expect_identical(six, umap_from_counts(counts, ntop = 100, seed = 1, threads = 1, n_components = 6))
  expect_identical(six, umap_from_counts(counts, ntop = 100, seed = 1, threads = 4, n_components = 6))
  expect_error(umap_from_counts(counts, n_components=1), "n_components")
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

test_that("UMAP axis annotation scores modules and labels axes by within-study correlation", {
  modules <- umap_marker_modules()
  base <- attr(modules, "baseline")
  prolif_genes <- modules$up[[match("proliferation", modules$module)]]
  genes <- c(base, prolif_genes, "VIM", "FN1", "CDH1", "OTHER")
  set.seed(2)
  n <- 40
  prolif <- rep(c(1, 8), each = n / 2)          # proliferation varies with the axis
  counts <- matrix(rpois(length(genes) * n, 50), length(genes), n,
                   dimnames = list(paste0("tx", seq_along(genes)), paste0("R", seq_len(n))))
  p <- genes %in% prolif_genes
  counts[p, ] <- rpois(sum(p) * n, rep(prolif * 20, each = sum(p)))
  scores <- umap_module_scores(counts, genes, modules)
  expect_identical(scores$sample, colnames(counts))
  expect_true("proliferation" %in% names(scores))
  # EMT has 3 of 15 genes present: skipped at the default fraction, scored when allowed
  expect_false("EMT" %in% names(scores))
  expect_true("EMT" %in% names(umap_module_scores(counts, genes, modules, min_gene_fraction = 0.2)))
  expect_false("interferon" %in% names(scores))   # no interferon genes present
  expect_gt(mean(scores$proliferation[21:40]) - mean(scores$proliferation[1:20]), 1)

  u <- data.table(sample = colnames(counts), BioProject = rep(c("P1", "P2"), n / 2),
                  `UMAP 1` = prolif + rnorm(n, sd = 0.1), `UMAP 2` = rnorm(n))
  assoc <- umap_axis_associations(u, scores, 2)
  expect_identical(nrow(assoc), 2L * (ncol(scores) - 1L))
  expect_identical(assoc[axis == 1][which.max(abs(within_study)), module], "proliferation")
  expect_identical(umap_axis_label(c("a", "b", "c"), c(0.1, -0.6, 0.3)), "b (-), c (+)")
  expect_identical(umap_axis_label(c("a", "b"), c(0.1, -0.15)), "b (-; weak), a (+; weak)")
  expect_error(umap_module_scores(counts[-(1:12), ], genes[-(1:12)], modules), "housekeeping")
})
