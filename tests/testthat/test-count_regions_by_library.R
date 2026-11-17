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

test_that("per-sample cache invalidates changed inputs, corruption and regions", {
  e <- "PRJEB46746-homo_sapiens"
  skip_if_not(file.exists(file.path(ORFik::config()["exp"], paste0(e, ".csv"))), "real test study not available")
  ns <- asNamespace("ORFik")
  original_path <- get("filepath", ns)
  original_region <- get("loadRegionCustom", ns)
  replace_binding <- function(name, value) {
    unlockBinding(name, ns); assign(name, value, ns); lockBinding(name, ns)
  }
  on.exit({replace_binding("filepath", original_path); replace_binding("loadRegionCustom", original_region)}, add = TRUE)
  for (um in c(FALSE, TRUE)) {
    d <- ORFik::read.experiment(e, validate = FALSE, output.env = new.env()); ORFik::uniqueMappers(d) <- um
    source <- original_path(d, "cov")
    skip_if_not(length(source) >= 2 && all(file.exists(source)), "coverage not available")
    root <- withr::local_tempdir(); paths <- file.path(root, basename(source))
    stopifnot(all(file.copy(source, paths)))
    replace_binding("filepath", function(df, type, ...) if (identical(type, "cov")) paths else original_path(df, type, ...))
    out <- file.path(root, "cache"); fresh <- file.path(root, "fresh")
    run <- function(cache, dir = out) count_regions_by_library(d, out.dir = dir, count_cache = cache, BPPARAM = BiocParallel::MulticoreParam(2))
    off <- run(FALSE, fresh); on <- run(TRUE)
    expect_identical(on, off)
    cache_files <- file.path(out, "QC_STATS", ".per_sample_counts", paste0(ORFik::runIDs(d), ".qs"))
    before <- lapply(cache_files, readBin, what = "raw", n = 1e8)
    Sys.setFileTime(cache_files, as.POSIXct("2000-01-01", tz = "UTC"))
    stamps <- file.info(cache_files)$mtime
    expect_identical(run(TRUE), off)
    expect_identical(file.info(cache_files)$mtime, stamps)
    file.copy(source[2], paths[1], overwrite = TRUE)
    updated <- run(TRUE); after <- lapply(cache_files, readBin, what = "raw", n = 1e8)
    expect_false(identical(before[[1]], after[[1]]))
    expect_identical(before[-1], after[-1])
    expect_gt(as.numeric(file.info(cache_files)$mtime[1]), as.numeric(stamps[1]))
    expect_identical(file.info(cache_files)$mtime[-1], stamps[-1])
    expect_identical(updated, run(FALSE, fresh))
    writeBin(charToRaw("bad"), cache_files[1]); expect_identical(run(TRUE), updated)
    expect_identical(lapply(cache_files, readBin, what = "raw", n = 1e8), after)
    replace_binding("loadRegionCustom", function(...) { x <- original_region(...); x[seq_len(max(0, length(x) - 1L))] })
    changed <- run(TRUE); expect_identical(changed, run(FALSE, fresh))
    expect_false(identical(lapply(cache_files, readBin, what = "raw", n = 1e8), after))
    replace_binding("loadRegionCustom", original_region); replace_binding("filepath", original_path)
  }
})
