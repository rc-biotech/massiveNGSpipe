test_that("repair_counttable_labels maps stale columns only by exact evidence", {
  skip_if_not_installed("SummarizedExperiment")
  proj <- tempfile("proj"); dir.create(proj); qdir <- tempfile("qc"); dir.create(qdir)
  runs <- c("SRR1", "SRR2", "SRR3"); lab <- c("RNA_WT_r1", "RNA_KO_r1", "RNA_KO_r2")
  fake_df <- data.frame(stage = "", rep = c(1, 1, 2), condition = c("WT", "KO", "KO"), fraction = "")
  local_mocked_bindings(read.experiment = function(...) fake_df, uniqueMappers = function(x) FALSE,
                        `uniqueMappers<-` = function(x, value) { attr(x, "um") <- value; x },
                        QCfolder = function(d) if (isTRUE(attr(d, "um"))) file.path(qdir, "none") else qdir,
                        bamVarName = function(df, skip.replicate = FALSE) if (skip.replicate) c("RNA_WT", "RNA_KO", "RNA_KO") else lab)
  local_mocked_bindings(runIDs = function(df) runs, .package = "ORFik")
  mk <- function(cols) SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:6, 2, 3, dimnames = list(NULL, cols))),
                                                                    colData = S4Vectors::DataFrame(SAMPLE = cols, row.names = cols))
  saveRDS(mk(c("RNA_WT_r1", "RNA_old_a", "RNA_old_b")), file.path(qdir, "countTable_cds.rds"))
  fwrite(data.table(time = "t", experiment = "E", Run = "SRR2", old = "RNA_old_a", new = "RNA_KO_r1"), file.path(proj, "relabel_log.csv"))
  cfg <- list(project = proj)
  r <- suppressMessages(repair_counttable_labels("E", cfg, exp_dir = qdir))
  expect_identical(r$status, rep("would_repair", 2)); expect_identical(r$method, c("relabel_log", "single_leftover"))
  r <- suppressMessages(repair_counttable_labels("E", cfg, exp_dir = qdir, dry_run = FALSE))
  se <- readRDS(file.path(qdir, "countTable_cds.rds"))
  expect_identical(colnames(se), lab); expect_identical(se$SAMPLE, c("RNA_WT", "RNA_KO", "RNA_KO"))
  expect_identical(SummarizedExperiment::assay(se)[, 2], 3:4)                      # counts move with their column
  expect_true(file.exists(file.path(proj, "counttable_repair_log.csv")))
  # two unexplained columns: refused, untouched
  saveRDS(mk(c("RNA_WT_r1", "x", "y")), file.path(qdir, "countTable_cds.rds")); unlink(file.path(proj, "relabel_log.csv"))
  r <- suppressMessages(repair_counttable_labels("E", cfg, exp_dir = qdir, dry_run = FALSE))
  expect_true(all(r$status == "unmapped_column")); expect_identical(colnames(readRDS(file.path(qdir, "countTable_cds.rds"))), c("RNA_WT_r1", "x", "y"))
  # explicit mapping resolves them
  r <- suppressMessages(repair_counttable_labels("E", cfg, mapping = data.table(label = c("x", "y"), Run = c("SRR3", "SRR2")), exp_dir = qdir, dry_run = FALSE))
  expect_identical(colnames(readRDS(file.path(qdir, "countTable_cds.rds"))), c("RNA_WT_r1", "RNA_KO_r2", "RNA_KO_r1"))
})
