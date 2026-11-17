#' Count tables for covRle libraries, parallel over libraries
#'
#' Same output as \code{ORFik::countTable_regions(df, lib.type = "cov", forceRemake = TRUE)} (identical
#' SummarizedExperiment objects, same files \code{<out.dir>/QC_STATS/countTable_<region>.qs}), but organised the fast
#' way (Håkon 2026-11-17; codex_tasks/count_step_20261117): the region ranges are loaded ONCE in the parent and shared with
#' forked workers; each worker imports one library and counts it against every region; the parent binds the per-library
#' vectors into the study tables. ORFik instead loads all libraries, then counts them serially per region (at most 4 useful
#' workers). Measured 48-sample study: 61 s -> 34 s (8 workers), 48 s -> 23 s (16 workers).
#' @param df an ORFik experiment at the wanted \code{uniqueMappers()} setting
#' @param out.dir,regions,library.names as in \code{ORFik::countTable_regions()}
#' @param BPPARAM BiocParallel backend (forking keeps the shared ranges copy-on-write)
#' @return named list of SummarizedExperiment, one per region (invisible)
#' @noRd
count_regions_by_library <- function(df, out.dir = ORFik::libFolder(df), regions = c("mrna", "leaders", "cds", "trailers"),
                                     library.names = ORFik::bamVarName(df), BPPARAM = BiocParallel::SerialParam()) {
  paths <- ORFik::filepath(df, "cov")
  stopifnot(length(paths) == length(library.names), all(file.exists(paths)))
  o <- asNamespace("ORFik")
  ranges <- lapply(regions, function(r) get("loadRegionCustom", o)(df, r, "tx", FALSE)); names(ranges) <- regions
  style <- GenomeInfoDb::seqinfo(df); um <- ORFik::uniqueMappers(df)
  per_lib <- BiocParallel::bplapply(seq_along(paths), function(i) {
    lib <- ORFik::fimport(paths[i], chrStyle = style, only_unique_mappers = um)
    lapply(ranges, function(x) ORFik::countOverlapsW(x, lib, weight = NULL))
  }, BPPARAM = BPPARAM)
  res <- lapply(seq_along(regions), function(j) {
    mat <- data.table::as.data.table(do.call(cbind, lapply(per_lib, `[[`, j)))
    colnames(mat) <- library.names
    get("SummarizedExperimentExp", o)(df, mat, ranges[[j]], library.names)
  }); names(res) <- regions
  qc <- file.path(out.dir, "QC_STATS"); dir.create(qc, showWarnings = FALSE, recursive = TRUE)
  for (j in seq_along(regions)) {
    dest <- file.path(qc, paste0("countTable_", regions[j], ".qs")); tmp <- tempfile(tmpdir = qc)
    qs2::qs_save(res[[j]], tmp); if (!file.rename(tmp, dest)) stop("Cannot write count table: ", dest)
  }
  invisible(res)
}
