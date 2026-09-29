# Unit tests for the pure helpers in R/pipeline_utils.R. is_config() is
# covered in test-pipeline_config.R (it's really a config-validation
# concept, tested alongside get_fun_name()). The RStudio-only navigate*()
# wrappers and the real-conda/real-`sed`-dependent functions are out of
# scope.

test_that("name_of_function captures the literal argument name(s) as written at its OWN call site", {
  expect_identical(name_of_function(mean), "mean")
  expect_identical(name_of_function(mean, sd), c("mean", "sd"))
})

test_that("dirExp builds the expected path for character exp names, by rel_dir and suffix", {
  # Note the doubled "//": dirExp() does file.path(config()["bam"], sub_dir)
  # first, and for an exp name with no _SSU/_RNA-seq/_disome suffix
  # sub_dir is "" -- file.path(x, "") really does produce a trailing "/",
  # confirmed directly rather than assumed. Not collapsed away downstream.
  testthat::local_mocked_bindings(
    config = function(...) c(bam = "/data/bam", fastq = "/data/fastq"),
    .package = "ORFik"
  )
  expect_identical(dirExp("study1-homo_sapiens", "aligned"),
                   "/data/bam//study1-homo_sapiens/aligned")
  expect_identical(dirExp("study1-homo_sapiens", "trim"),
                   "/data/bam//study1-homo_sapiens/trim")
})

test_that("dirExp: an unrecognized rel_dir silently falls back to the raw fastq dir (documented current behavior)", {
  testthat::local_mocked_bindings(
    config = function(...) c(bam = "/data/bam", fastq = "/data/fastq"),
    .package = "ORFik"
  )
  expect_identical(dirExp("study1-homo_sapiens", "totally_unknown_rel_dir"),
                   "/data/fastq//study1-homo_sapiens")
})

test_that("file_statistics_internal summarizes count/size per (dir_type, format) pair from a real tempdir tree", {
  l <- tempfile(); dir.create(file.path(l, "aligned"), recursive = TRUE)
  writeLines(strrep("x", 1000), file.path(l, "aligned", "a.bam"))
  writeLines(strrep("x", 2000), file.path(l, "aligned", "b.bam"))

  result <- file_statistics_internal(dir_types = "aligned", formats = "bam", type = "processed", l = l)
  expect_identical(result$n_files[1], 2L)
  expect_true(result$total_size_GB[1] >= 0)
})

test_that("file_statistics_internal: no matching files gives a real 0, not NA/error", {
  l <- tempfile(); dir.create(l, recursive = TRUE)
  result <- file_statistics_internal(dir_types = "does_not_exist", formats = "bam", type = "processed", l = l)
  expect_identical(result$n_files[1], 0L)
  expect_identical(result$total_size_GB[1], 0)
})
