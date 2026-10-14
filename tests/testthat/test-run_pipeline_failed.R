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

test_that("pending_exps is per experiment for multi-step stages (no NA from step-major indexing)", {
  done <- list("E1" = c(start = TRUE, fetch = FALSE), "E2" = c(start = TRUE, fetch = TRUE))
  local_mocked_bindings(all_substeps_done = function(config, steps, exp) unname(done[[exp]][steps]))
  expect_identical(pending_exps(list(), c("start", "fetch"), c("E1", "E2")), "E1")
  expect_identical(pending_exps(list(), c("start", "fetch"), "E1"), "E1")
  d <- tempfile("err"); dir.create(d); saveRDS("x", file.path(d, "E1.rds"))
  expect_true(all_pending_failed(list(error_dir = d), pending_exps(list(), c("start", "fetch"), "E1")))
})
