test_that("export_sucessful_metadata refuses to shrink FINAL drastically; libtypes default follows preset", {
  d <- tempfile("v"); dir.create(d)
  out <- file.path(d, "FINAL_LIST.csv"); nxt <- file.path(d, "next.csv")
  fwrite(data.table(Run = paste0("R", 1:10), LIBRARYTYPE = "RNA", KEEP = TRUE), out)
  files <- data.table(Run = paste0("R", 1:10), LIBRARYTYPE = "RNA", KEEP = TRUE)
  expect_error(suppressMessages(export_sucessful_metadata(files, "RFP", out, nxt)), "Refusing")
  expect_identical(nrow(fread(out)), 10L)                       # untouched
  expect_false(file.exists(nxt))
  suppressMessages(export_sucessful_metadata(files, "RNA", out, nxt))
  expect_identical(nrow(fread(out)), 10L)
  expect_identical(formals(pipeline_validate_metadata)$libtypes, quote(libtype_long_to_short(config$preset)))
  expect_identical(libtype_long_to_short(c("Ribo-seq", "RNA-seq", "disome")), c("RFP", "RNA", "disome"))
})
