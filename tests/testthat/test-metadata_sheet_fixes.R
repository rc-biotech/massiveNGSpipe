test_that("sheet_apply_fixes edits only targeted cells, skips changed or missing runs, logs", {
  base <- tempfile("fix"); dir.create(file.path(base, "NGS_pipeline"), recursive = TRUE)
  f <- file.path(base, "NGS_pipeline", "next_round_manual.csv")
  fwrite(data.table(Run = c("R1", "R2", "R3"), CELL_LINE = c("k562", "primary_fibroblast", "HeLa"),
                    TISSUE = c("blood", "", "cervix"), REPLICATE = 1:3), f)
  fixes <- data.table(Run = c("R1", "R2", "R2", "R3", "R9"), column = c("CELL_LINE", "CELL_LINE", "TISSUE", "CELL_LINE", "CELL_LINE"),
                      old = c("k562", "primary_fibroblast", "", "hela", "x"), new = c("K562", "primary", "fibroblast", "HeLa", "y"))
  dry <- suppressMessages(sheet_apply_fixes(fixes, base_dir = base, use_google = FALSE))
  expect_identical(dry$status, c("would_apply", "would_apply", "would_apply", "changed_meanwhile", "run_missing"))
  expect_identical(fread(f)$CELL_LINE[1], "k562")                                   # dry run writes nothing
  res <- suppressMessages(sheet_apply_fixes(fixes, base_dir = base, use_google = FALSE, note = "test", dry_run = FALSE))
  s <- fread(f)
  expect_identical(s$CELL_LINE, c("K562", "primary", "HeLa"))
  expect_identical(s$TISSUE, c("blood", "fibroblast", "cervix"))
  expect_identical(s$REPLICATE, 1:3)                                                 # untouched column intact
  log <- fread(file.path(base, "NGS_pipeline", "metadata_fix_log.csv"))
  expect_identical(nrow(log), 3L)
  expect_true(all(log$note == "test"))
  # a second application finds the new values: the old values no longer match
  again <- suppressMessages(sheet_apply_fixes(fixes[1:2], base_dir = base, use_google = FALSE))
  expect_true(all(again$status == "changed_meanwhile"))
})

test_that("locale guard: damaged values detected, bare decimals refused on upload", {
  dt <- data.table(Run = paste0("R", 1:5), TIMEPOINT = c("0,5", "3333333333", "0.5", "CT04", "30min"), KEEP = TRUE)
  chk <- metadata_locale_check(dt)
  expect_identical(chk$problem, c("locale_damaged", "locale_damaged", "decimal_at_risk"))
  base <- tempfile("fix"); dir.create(file.path(base, "NGS_pipeline"), recursive = TRUE)
  fwrite(dt, file.path(base, "NGS_pipeline", "next_round_manual.csv"))
  expect_error(sheet_apply_fixes(data.table(Run = "R1", column = "TIMEPOINT", old = "0,5", new = "0.5"),
                                 base_dir = base, use_google = FALSE), "unit")
  ok <- suppressMessages(sheet_apply_fixes(data.table(Run = "R1", column = "TIMEPOINT", old = "0,5", new = "30min"),
                                           base_dir = base, use_google = FALSE))
  expect_identical(ok$status, "would_apply")
})
