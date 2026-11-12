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

test_that("QC cache reads keep positions (chunks contiguous, n = 1 works)", {
  for (n in c(1, 2, 5, 17, 1000)) {
    ch <- massiveNGSpipe:::qc_chunks(n, 64)
    expect_identical(unlist(ch, use.names = FALSE), seq_len(n)); expect_lte(length(ch), 64)
  }
  f <- tempfile(fileext = rep(".rds", 7)); for (i in seq_along(f)) saveRDS(list(sig = i), f[i])
  r <- massiveNGSpipe:::qc_read_caches(c(f, file.path(tempdir(), "missing.rds")), 3)
  expect_identical(vapply(r[1:7], `[[`, 0, "sig"), as.numeric(1:7)); expect_null(r[[8]])
  expect_length(massiveNGSpipe:::qc_read_caches(f[1], 16), 1)
})

test_that("prefetched cache lookup returns entries by path, missing ones read", {
  f <- tempfile(fileext = rep(".rds", 3)); for (i in seq_along(f)) saveRDS(list(sig = i), f[i])
  pre <- setNames(massiveNGSpipe:::qc_read_caches(f[1:2], 2), f[1:2])
  r <- massiveNGSpipe:::qc_cached(f[c(3, 1, 2)], pre, 2)
  expect_identical(vapply(r, `[[`, 0, "sig"), c(3, 1, 2))
})

test_that("terminal bias is cached per file and equals the uncached result", {
  d <- file.path(tempdir(), "tb"); unlink(d, recursive = TRUE); dir.create(file.path(d, "s1"), recursive = TRUE)
  w <- function(run, side, x) data.table::fwrite(data.table::data.table(x, count = c(5, 3, 2)), file.path(d, "s1", paste0(run, "_", side, "_bias.tsv")), sep = "\t")
  for (r in c("SRR1", "SRR2")) { w(r, "5prime", data.table::data.table(read_5prime_3nt = c("ACG", "TTT", "GAA")))
    w(r, "3prime", data.table::data.table(read_3prime_3nt = c("CCA", "AGT", "AAA"))) }
  cd <- file.path(tempdir(), "tbc"); unlink(cd, recursive = TRUE); dir.create(cd)
  ref <- massiveNGSpipe:::qc_terminal_bias(d, c("SRR1", "SRR2"))
  a <- massiveNGSpipe:::qc_terminal_bias(d, c("SRR1", "SRR2"), cd, 1)
  expect_true(file.exists(file.path(cd, "_terminal_bias.rds")))
  b <- massiveNGSpipe:::qc_terminal_bias(d, c("SRR1", "SRR2"), cd, 1)
  expect_equal(a, ref); expect_equal(b, ref); expect_equal(ref[Run == "SRR1"][["A5'"]], 0.5)
})
