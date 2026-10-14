test_that("platform_guard skips long-read runs and records them", {
  d <- tempfile("proj"); dir.create(d)
  x <- data.table(Run = paste0("R", 1:4), BioProject = "P", Platform = c("ILLUMINA", "OXFORD_NANOPORE", "BGISEQ", "PACBIO_SMRT"))
  y <- suppressMessages(platform_guard(x, list(project = d)))
  expect_identical(y$Run, c("R1", "R3"))
  e <- fread(file.path(d, "excluded_runs_platform.csv"))
  expect_identical(e$Run, c("R2", "R4")); expect_true(all(e$reason == "requires_longread_processing"))
  expect_identical(platform_guard(x[, .(Run)], list(project = d)), x[, .(Run)])     # no Platform column: unchanged
})

test_that("platform_guard skips ABI SOLiD colour-space runs with their own reason", {
  d <- tempfile("proj"); dir.create(d)
  x <- data.table(Run = paste0("R", 1:3), Platform = c("ILLUMINA", "ABI_SOLID", "OXFORD_NANOPORE"))
  y <- suppressMessages(platform_guard(x, list(project = d)))
  expect_identical(y$Run, "R1")
  e <- fread(file.path(d, "excluded_runs_platform.csv"))
  expect_identical(e$reason, c("requires_colorspace_decoding", "requires_longread_processing"))
})
