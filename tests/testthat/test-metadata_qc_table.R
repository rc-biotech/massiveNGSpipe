test_that("qc helpers: metadata prep, file signature and column order/rounding", {
  m <- data.table(Run = c("SRR1", "SRR2"), study_accession = "PRJ1", ScientificName = "Homo sapiens", KEEP = TRUE, CT = "x")
  p <- qc_prepare_metadata(m)
  expect_identical(p$study, rep("PRJ1-homo_sapiens", 2)); expect_false(any(c("KEEP", "CT") %in% names(p)))
  expect_true(all(p$mode == "online"))
  f <- tempfile(); writeLines("a", f); s1 <- qc_file_sig(f); Sys.sleep(1.1); writeLines("ab", f)
  expect_false(identical(s1, qc_file_sig(f))); expect_identical(qc_file_sig(NA_character_), "NA NA")
  prev <- tempfile(fileext = ".csv"); fwrite(data.table(Run = "x", b = 1, a = 2), prev)
  out <- qc_order_and_round(data.table(a = 1.23456, Run = "SRR1", b = 2L, z = "k"), prev)
  expect_identical(names(out), c("Run", "b", "a", "z")); expect_identical(out$a, 1.23); expect_identical(out$b, 2L)
})
test_that("qc_cov_metrics returns an NA row for a missing coverage file", {
  r <- qc_cov_metrics(tempfile(), list())
  expect_true(is.na(r$Total_Cov_Signal))
})

test_that("qc_merge_subset keeps other organisms' rows and replaces the subset", {
  f <- tempfile(fileext = ".csv")
  data.table::fwrite(data.table::data.table(Run = c("a", "b", "c", "gone"), ScientificName = c("Homo sapiens", "Mus musculus", "Homo sapiens", "Mus musculus"),
                                            x = c(1, 2, 3, 4)), f)
  out <- data.table::data.table(Run = c("a", "c"), ScientificName = "Homo sapiens", x = c(10, 30), y = c("p", "q"))
  res <- qc_merge_subset(out, f, "Homo sapiens", all_runs = c("a", "b", "c"))
  expect_equal(res$Run, c("a", "b", "c"))
  expect_equal(res$x, c(10, 2, 30))
  expect_true(is.na(res[Run == "b", y]))
  expect_identical(qc_merge_subset(out, tempfile(), "Homo sapiens", c("a", "c")), out)
})
