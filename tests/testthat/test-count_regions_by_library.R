# count_regions_by_library() must give the same SummarizedExperiments as ORFik::countTable_regions() (real 2-sample study)
test_that("count_regions_by_library() is identical to ORFik::countTable_regions() for both mapper branches", {
  e <- "PRJEB46746-homo_sapiens"
  skip_if_not(file.exists(file.path(ORFik::config()["exp"], paste0(e, ".csv"))), "real test study not available")
  for (um in c(FALSE, TRUE)) {
    df <- ORFik::read.experiment(e, validate = FALSE, output.env = new.env()); ORFik::uniqueMappers(df) <- um
    out_new <- withr::local_tempdir(); out_old <- withr::local_tempdir()
    new <- count_regions_by_library(df, out.dir = out_new, BPPARAM = BiocParallel::MulticoreParam(2))
    old <- ORFik::countTable_regions(df, out.dir = out_old, lib.type = "cov", forceRemake = TRUE,
                                     BPPARAM = BiocParallel::SerialParam())
    for (r in names(new)) {
      expect_identical(new[[r]], old[[r]])
      expect_identical(qs2::qs_read(file.path(out_new, "QC_STATS", paste0("countTable_", r, ".qs"))),
                       qs2::qs_read(file.path(out_old, "QC_STATS", paste0("countTable_", r, ".qs"))))
    }
  }
})
