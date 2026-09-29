# Unit tests for the pure helpers in R/pipeline_preset_steps.R
# (libtype_long_to_short/short_to_long, library_strategies_allowed,
# safe_se_cbind). The pipe_*()/pipeline_collection_*()/pipeline_merge_*()
# orchestrators are out of scope -- they need real ORFik experiments,
# BAMs, and STAR/BiocParallel infra.

test_that("libtype_long_to_short / libtype_short_to_long recode the known pairs and pass through unknowns", {
  expect_identical(libtype_long_to_short("Ribo-seq"), "RFP")
  expect_identical(libtype_long_to_short("RNA-seq"), "RNA")
  expect_identical(libtype_long_to_short("disome"), "disome") # passthrough
  expect_identical(libtype_long_to_short(c("Ribo-seq", "RNA-seq")), c("RFP", "RNA"))

  expect_identical(libtype_short_to_long("RFP"), "Ribo-seq")
  expect_identical(libtype_short_to_long("RNA"), "RNA-seq")
  expect_identical(libtype_short_to_long("disome"), "disome")
})

test_that("libtype_long_to_short / libtype_short_to_long round-trip for the two known pairs", {
  expect_identical(libtype_short_to_long(libtype_long_to_short("Ribo-seq")), "Ribo-seq")
  expect_identical(libtype_short_to_long(libtype_long_to_short("RNA-seq")), "RNA-seq")
})

test_that("library_strategies_allowed returns the expected fixed set", {
  expect_setequal(library_strategies_allowed(),
                  c("RNA-Seq", "Ribo-seq", "ssRNA-seq", "ncRNA-Seq", "miRNA-Seq", "OTHER"))
})

test_that("safe_se_cbind cbinds matching-nrow SummarizedExperiments on the happy path", {
  se1 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2)))
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(5:8, nrow = 2)))
  result <- safe_se_cbind(list(a = se1, b = se2))
  expect_identical(ncol(result), 4L)
})

test_that("safe_se_cbind reports a clear diagnostic error for mismatched row counts", {
  se_2row <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2)))
  se_3row <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:6, nrow = 3)))
  err <- tryCatch(
    safe_se_cbind(list(good1 = se_2row, good2 = se_2row, bad = se_3row)),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "bad", fixed = TRUE)
  expect_match(conditionMessage(err), "nrow=3", fixed = TRUE)
})
