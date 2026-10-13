# Cached, incremental meta-motif tables for organism collections.
#
# Product (unchanged names/shape): resFolder/meta_collection_tables/meta_collection_table_<anchor>.fst,
# rows = 600 positions around the anchor (300 up, 300 down), one double column per Run = summed
# coverage over the anchor windows. Anchors: TIS, TTS and every amino-acid pair of letters_vec.
#
# Differences from the original meta_meta_motif_motif():
# * Anchor windows are selected once per organism and cached (motif_regions.qs); subsampling to
#   max_ranges uses a fixed seed (the original sampled unseeded, so tables were not reproducible).
# * Each library is loaded once and scored for all anchors with one coveragePerTiling() call,
#   summed per position via a matrix (identical values, ~4.5x faster than per-anchor scoring).
# * Manifest of library inputs: reruns only score new/changed libraries; removed libraries are
#   dropped. Tables are small, so they are rewritten whole (.tmp + rename): plain fst, no segments.
# Benchmarks: ~/shared_workspace/hakon_playground/megafst/notes/04_motif.md

MOTIF_FORMAT_VERSION <- 1L

#' Select and cache the anchor windows of one collection
#' @param df collection ORFik experiment
#' @param anchors character, "TIS", "TTS" and amino-acid pair strings
#' @param extension integer, nt up/downstream of the anchor
#' @param max_ranges integer, subsample anchors above this count
#' @param seed integer, subsampling seed
#' @return named list of GRangesList (all windows exactly 2 * extension wide)
#' @noRd
motif_select_regions <- function(df, anchors, extension = 300, max_ranges = 30000, seed = 1) {
  org_short <- gsub(" ", "_", tolower(organism(df)))
  txdb <- loadTxdb(df)
  subset <- integer()
  if (org_short == "homo_sapiens") {
    subset <- canonical_isoforms(df)
  } else try(subset <- filterTranscripts(txdb, NULL, 30, NULL, longestPerGene = TRUE), silent = TRUE)
  all_cds <- loadRegion(txdb, "cds")
  mrna <- loadRegion(txdb, "mrna")
  cds <- all_cds[names(all_cds) %in% subset]
  mrna <- mrna[names(mrna) %in% subset]
  stopifnot(length(cds) == length(mrna))
  mrna <- extendTrailersUntil(extendLeadersUntil(mrna, cds, extension), cds, extension)
  codon_seqs <- NULL
  out <- list()
  for (anchor_site in anchors) {
    anchor <- if (anchor_site == "TIS") {
      startSites(cds, TRUE, TRUE, TRUE)
    } else if (anchor_site == "TTS") {
      stopSites(cds, TRUE, TRUE, TRUE)
    } else {
      if (is.null(codon_seqs)) {
        seqs <- txSeqsFromFa(cds, df, is.sorted = TRUE, keep.names = TRUE)
        codon_seqs <- Biostrings::translate(seqs, if.fuzzy.codon = "solve")
      }
      hits <- IRangesList(Biostrings::vmatchPattern(anchor_site, codon_seqs))
      names(hits) <- seq_along(hits)
      hits <- hits[lengths(hits) > 0]
      width(hits) <- 1
      end(hits) <- IntegerList((end(hits) * 3) - 2L)
      start(hits) <- IntegerList((start(hits) * 3) - 2L)
      unlistGrl(pmapFromTranscriptF(hits, cds, removeEmpty = TRUE))
    }
    n_anchor <- length(anchor)
    windows <- windowPerGroup(anchor, tx = mrna, upstream = extension, downstream = extension - 1)
    windows <- windows[widthPerGroup(windows) == extension * 2]
    if (length(windows) > max_ranges) {
      set.seed(seed)
      windows <- windows[sample(length(windows), max_ranges, replace = FALSE)]
    }
    message(sprintf("--- %s: %d anchors, %d windows used", anchor_site, n_anchor, length(windows)))
    out[[anchor_site]] <- windows
  }
  out
}

#' Summed coverage per window position for every anchor, one library, one load
#' @param cov_path character, covRle file
#' @param windows_all GRangesList, all anchors' windows concatenated
#' @param anchor_index integer, anchor of each window in windows_all
#' @param anchors character, anchor names (levels of anchor_index)
#' @param width integer, window width
#' @return named list of numeric(width), one per anchor
#' @noRd
motif_score_library <- function(cov_path, windows_all, anchor_index, anchors, width) {
  cov <- fimport(cov_path)
  v <- unlist(coveragePerTiling(windows_all, cov, is.sorted = TRUE), use.names = FALSE)
  as_int <- is.integer(S4Vectors::runValue(v))   # integer coverage sums stay integer, as in the original
  m <- matrix(as.numeric(v), nrow = width)
  stopifnot(ncol(m) == length(anchor_index))
  setNames(lapply(seq_along(anchors), function(a) {
    s <- rowSums(m[, anchor_index == a, drop = FALSE])
    if (as_int && all(abs(s) <= .Machine$integer.max)) as.integer(s) else s
  }), anchors)
}

#' Build or incrementally update one collection's meta-motif tables
#' @param df collection ORFik experiment (all_samples-<Organism>)
#' @param letters_vec character, amino acids whose ordered pairs are anchors
#' @param extension,max_ranges,seed see motif_select_regions()
#' @param force logical, recompute regions and all libraries
#' @param batch_size integer, libraries scored per batch
#' @param BPPARAM BiocParallel param over libraries
#' @param out_dir character, output directory (default the collection's meta_collection_tables)
#' @inheritParams wait_for_memory
#' @return invisibly the plan mode ("none", "incremental", "full")
#' @export
meta_motif_update <- function(df, letters_vec = c("P", "R"), extension = 300, max_ranges = 30000,
                              seed = 1, force = FALSE, batch_size = 500, keep_free_gb = 32,
                              timeout_minutes = 120, covrle_memory_ratio = 11,
                              out_dir = file.path(dirname(QCfolder(df)), "meta_collection_tables"),
                              BPPARAM = BiocParallel::MulticoreParam(8, exportglobals = FALSE, log = FALSE)) {
  cov_paths <- filepath(df, "cov", base_folders = libFolder(df, mode = "all"), suffix_stem = c("_pshifted", ""))
  lib <- megafst_library_table(runIDs(df), cov_paths)
  org_dir <- out_dir
  dir.create(org_dir, showWarnings = FALSE, recursive = TRUE)
  lock <- filelock::lock(file.path(org_dir, ".motif.lock"), exclusive = TRUE, timeout = 0)
  if (is.null(lock)) stop("Another meta-motif build holds the lock for ", org_dir)
  on.exit(filelock::unlock(lock), add = TRUE)
  megafst_remove_stale_tmp(org_dir)

  pairs <- as.vector(outer(letters_vec, letters_vec, paste0))
  anchors <- c("TIS", "TTS", pairs)
  files <- setNames(file.path(org_dir, paste0("meta_collection_table_", anchors, ".fst")), anchors)
  regions_path <- file.path(org_dir, "motif_regions.qs")
  params_path <- file.path(org_dir, "motif_params.rds")
  manifest_path <- file.path(org_dir, "motif_manifest.fst")
  annot <- c(df@txdb, df@fafile)
  params <- list(format_version = MOTIF_FORMAT_VERSION, anchors = anchors, extension = extension,
                 max_ranges = max_ranges, seed = seed, annotation = annot,
                 annotation_mtime = as.numeric(file.info(annot)$mtime))
  params_old <- if (file.exists(params_path)) readRDS(params_path) else NULL
  regions_valid <- !force && !is.null(params_old) && identical(params_old, params) && file.exists(regions_path)
  if (regions_valid) {
    windows <- qs2::qs_read(regions_path)
  } else {
    message("-- Selecting motif regions (cached for later runs)")
    windows <- motif_select_regions(df, anchors, extension, max_ranges, seed)
    tmp <- paste0(regions_path, ".tmp"); qs2::qs_save(windows, tmp); stopifnot(file.rename(tmp, regions_path))
  }
  old <- if (regions_valid && file.exists(manifest_path) && all(file.exists(files)))
    fst::read_fst(manifest_path, as.data.table = TRUE) else NULL
  if (is.null(old)) {
    mode <- "full"; todo <- lib$Run; keep <- character()
  } else {
    m <- merge(lib, old, by = "Run", all.x = TRUE, suffixes = c("", "_old"), sort = FALSE)
    todo <- m[is.na(cov_path_old) | cov_path != cov_path_old | cov_mtime != cov_mtime_old |
                cov_size != cov_size_old]$Run
    keep <- setdiff(intersect(lib$Run, old$Run), todo)
    mode <- if (length(todo) || !setequal(old$Run, lib$Run) || !identical(old$Run, lib$Run)) "incremental" else "none"
  }
  message(sprintf("-- meta-motif plan: %s (%d to score, %d kept, %d removed)", mode, length(todo),
                  length(keep), if (is.null(old)) 0L else length(setdiff(old$Run, lib$Run))))
  if (mode == "none") return(invisible("none"))

  windows_all <- do.call(c, unname(windows))
  anchor_index <- rep(seq_along(anchors), lengths(windows[anchors]))
  width <- 2L * extension
  scores <- setNames(vector("list", length(anchors)), anchors)
  for (a in anchors) scores[[a]] <- list()
  todo_rows <- which(lib$Run %in% todo)
  batches <- split(todo_rows, ceiling(seq_along(todo_rows) / batch_size))
  workers <- BiocParallel::bpnworkers(BPPARAM)
  for (b in seq_along(batches)) {
    rows <- batches[[b]]
    need <- workers * max(lib$cov_size[rows]) * covrle_memory_ratio / 1024^3 +
      workers * length(windows_all) * width * 8 * 3 / 1024^3
    wait_for_memory(need, sprintf("meta-motif batch %d/%d", b, length(batches)), keep_free_gb, timeout_minutes)
    message(sprintf("--- scoring batch %d/%d: %d libraries", b, length(batches), length(rows)))
    res <- BiocParallel::bplapply(lib$cov_path[rows], motif_score_library, windows_all = windows_all,
                                  anchor_index = anchor_index, anchors = anchors, width = width,
                                  BPPARAM = BPPARAM)
    for (k in seq_along(rows)) for (a in anchors) scores[[a]][[lib$Run[rows[k]]]] <- res[[k]][[a]]
    rm(res); gc(FALSE)
  }
  for (a in anchors) {
    new_cols <- if (length(scores[[a]])) as.data.table(scores[[a]]) else NULL
    old_cols <- if (length(keep)) fst::read_fst(files[[a]], columns = keep, as.data.table = TRUE) else NULL
    tab <- if (is.null(old_cols)) new_cols else if (is.null(new_cols)) old_cols else cbind(old_cols, new_cols)
    stopifnot(nrow(tab) == width, setequal(names(tab), lib$Run))
    setcolorder(tab, lib$Run)
    tmp <- paste0(files[[a]], ".tmp"); fst::write_fst(tab, tmp); stopifnot(file.rename(tmp, files[[a]]))
  }
  tmp <- paste0(params_path, ".tmp"); saveRDS(params, tmp); stopifnot(file.rename(tmp, params_path))
  tmp <- paste0(manifest_path, ".tmp"); fst::write_fst(lib, tmp); stopifnot(file.rename(tmp, manifest_path))
  invisible(mode)
}
