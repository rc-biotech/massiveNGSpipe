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
#' @param count_cache reuse validated per-sample coverage counts (default FALSE)
#' @return named list of SummarizedExperiment, one per region (invisible)
#' @noRd
count_regions_by_library <- function(df, out.dir = ORFik::libFolder(df), regions = c("mrna", "leaders", "cds", "trailers"),
                                     library.names = ORFik::bamVarName(df), BPPARAM = BiocParallel::SerialParam(),
                                     count_cache = FALSE) {
  paths <- ORFik::filepath(df, "cov")
  stopifnot(length(paths) == length(library.names), all(file.exists(paths)))
  o <- asNamespace("ORFik")
  ranges <- lapply(regions, function(r) get("loadRegionCustom", o)(df, r, "tx", FALSE)); names(ranges) <- regions
  style <- GenomeInfoDb::seqinfo(df); um <- ORFik::uniqueMappers(df)
  if (count_cache) {
    ids <- as.character(ORFik::runIDs(df))
    if (length(ids) != length(paths) || anyNA(ids) || any(!nzchar(ids)) ||
        anyDuplicated(ids) || any(basename(ids) != ids) || any(ids %in% c(".", "..")))
      stop("Cache requires unique, safe Run identifiers")
    cache_dir <- file.path(out.dir, "QC_STATS", ".per_sample_counts")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    region_hash <- digest::digest(ranges, algo = "xxhash64")
  }
  per_lib <- BiocParallel::bplapply(seq_along(paths), function(i) {
    if (count_cache) {
      info <- file.info(paths[i])
      key <- list(schema = 1L, path = normalizePath(paths[i], mustWork = TRUE),
                  size = info$size, mtime = as.numeric(info$mtime),
                  content = digest::digest(file = paths[i], algo = "xxhash64"),
                  regions = region_hash, orfik = as.character(utils::packageVersion("ORFik")),
                  lib.type = "cov", chrStyle = style, only_unique_mappers = um)
      cache_file <- file.path(cache_dir, paste0(ids[i], ".qs"))
      entry <- tryCatch(qs2::qs_read(cache_file), error = function(e) NULL)
      valid <- is.list(entry) && identical(names(entry), c("key", "counts")) &&
        identical(entry$key, key) && is.list(entry$counts) &&
        identical(names(entry$counts), names(ranges)) &&
        all(vapply(seq_along(ranges), function(j) {
          x <- entry$counts[[j]]
          is.numeric(x) && is.null(dim(x)) && length(x) == length(ranges[[j]]) &&
            (is.null(names(x)) || identical(names(x), names(ranges[[j]]))) &&
            all(is.finite(x))
        }, logical(1)))
      if (valid) return(entry$counts)
    }
    lib <- ORFik::fimport(paths[i], chrStyle = style, only_unique_mappers = um)
    counts <- lapply(ranges, function(x) ORFik::countOverlapsW(x, lib, weight = NULL))
    if (count_cache) {
      tmp <- tempfile(tmpdir = cache_dir)
      on.exit(unlink(tmp), add = TRUE)
      qs2::qs_save(list(key = key, counts = counts), tmp)
      if (!file.rename(tmp, cache_file)) stop("Cannot write count cache: ", cache_file)
    }
    counts
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
