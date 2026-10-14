test_that("all_pending_failed: a stage loop stops waiting when every pending study has an error record", {
  d <- tempfile("err"); config <- list(error_dir = d)
  expect_false(all_pending_failed(config, c("PRJ1-homo_sapiens")))          # no error dir yet
  dir.create(d); saveRDS("x", file.path(d, "PRJ1.rds"))                     # records keyed by accession
  expect_true(all_pending_failed(config, "PRJ1-homo_sapiens"))
  expect_false(all_pending_failed(config, c("PRJ1-homo_sapiens", "PRJ2-homo_sapiens")))
  saveRDS("x", file.path(d, "PRJ2-homo_sapiens.rds"))                       # or by experiment name
  expect_true(all_pending_failed(config, c("PRJ1-homo_sapiens", "PRJ2-homo_sapiens")))
  expect_false(all_pending_failed(config, character()))
})
