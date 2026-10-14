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

test_that("sheet_remove_runs drops only the given rows; experiment_remove_runs rewrites the CSV and resets merge flags", {
  base <- tempfile("rm"); dir.create(file.path(base, "NGS_pipeline"), recursive = TRUE)
  f <- file.path(base, "NGS_pipeline", "next_round_manual.csv")
  fwrite(data.table(Run = c("R1", "R2", "R3"), BioProject = "P", KEEP = TRUE), f)
  res <- suppressMessages(sheet_remove_runs(c("R2", "R9"), base_dir = base, use_google = FALSE, dry_run = FALSE))
  expect_identical(res$status, c("removed", "run_missing"))
  expect_identical(fread(f)$Run, c("R1", "R3"))
  expect_identical(fread(file.path(base, "NGS_pipeline", "metadata_removed_rows.csv"))$Run, "R2")
  # experiment
  ed <- tempfile("exp"); dir.create(ed)
  writeLines(c('"name","P-x","","","","",""', '"gff","g.gtf.db","","","organism","X",""', '"fasta","g.fa","","","author","A",""',
               '"libtype","stage","rep","condition","fraction","filepath","Run"',
               '"RFP","leaf","1","WT","","/x/R1.bam","R1"', '"RFP","leaf","2","WT","","/x/R2.bam","R2"',
               '"RFP","leaf","3","WT","","/x/R3.bam","R3"'), file.path(ed, "P-x.csv"))
  steps <- c("pshifted", "merged_lib", "pcounts"); fl <- setNames(file.path(base, "flags", steps), steps)
  for (s in fl) { dir.create(s, recursive = TRUE); saveRDS(TRUE, file.path(s, "P-x.rds")) }
  config <- structure(list(flag = fl), class = c("list", "mNGSp_config"))
  testthat::local_mocked_bindings(is_config = function(x) TRUE)
  suppressMessages(experiment_remove_runs("P-x", "R2", config, exp_dir = ed, dry_run = FALSE))
  l <- readLines(file.path(ed, "P-x.csv"))
  expect_identical(length(l), 6L); expect_false(any(grepl("R2", l)))
  expect_true(file.exists(file.path(fl["pshifted"], "P-x.rds")))
  expect_false(file.exists(file.path(fl["merged_lib"], "P-x.rds")))
  expect_false(file.exists(file.path(fl["pcounts"], "P-x.rds")))
})
