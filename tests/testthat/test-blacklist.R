test_that("blacklist: old one-column format, study and run rows, filter and add", {
  d <- tempfile("bl"); dir.create(d); f <- file.path(d, "BLACKLIST.csv")
  fwrite(data.table(id = c("PRJ1", "PRJ1", "  ")), f)                          # old format + duplicate + blank
  config <- list(blacklist = f)
  b <- blacklist_read(config, sync = FALSE)
  expect_identical(names(b), c("id", "run_id", "reason"))
  expect_identical(nrow(b), 1L)
  expect_true(is.na(b$run_id))
  b <- suppressMessages(blacklist_add(config, id = c("PRJ2", "PRJ2", "PRJ1"), run_id = c("R5", "R6", NA), reason = "3'Ribo-seq"))
  expect_identical(nrow(b), 3L)                                                  # PRJ1 whole-study row not duplicated
  dt <- data.table(Run = paste0("R", 1:7), BioProject = c("PRJ1", "PRJ1", "PRJ2", "PRJ2", "PRJ2", "PRJ2", "PRJ3"))
  out <- blacklist_filter(dt, blacklist_read(config, sync = FALSE))
  expect_identical(out$Run, c("R3", "R4", "R7"))                                 # whole PRJ1 + runs R5, R6 removed
  expect_identical(attr(out, "removed"), 4L)
  expect_identical(nrow(fread(f)), 3L)
})
