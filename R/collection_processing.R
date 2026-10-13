#' Create megafst by coverage tiles
#'
#' Read + or - strand covRle over all libraries and process
save_coverage_tiles <- function(filepath, chromosomes, chromosome_lengths,
                                output_dir, tile_size, colnames, verbose = TRUE,
                                remake = FALSE, BPPARAM) {
  message("-- Loading file: ", basename(filepath))
  covrle_list_by_strand <- read_RDSQS(filepath, nthread = 5)
  message("-- Processing all chromosomes: ", basename(filepath))
  index_dt <- rbindlist(lapply(chromosomes, function(chr) {
    fst::threads_fst(5, reset_after_fork = FALSE)

    total <- chromosome_lengths[chr]
    n_tiles <- ceiling(total / tile_size)
    indices <- seq_len(n_tiles)
    if (verbose) message("Processing chromosome: ", chr, " (Tile pages: ", n_tiles, ")")
    part_files <- file.path(output_dir, paste0("coverage_", chr, "_part", indices))
    if (attr(covrle_list_by_strand, "strand") == "+") {
      file_paths <- paste0(part_files, "_forward",".fst")
    } else {
      file_paths <- paste0(part_files, "_reverse",".fst")
    }

    irl <- List(lapply(indices, function(i) {
      idx_from <- ((i - 1) * tile_size) + 1
      idx_to <- min(i * tile_size, total)
      IRanges(idx_from, idx_to)
    }))
    chr_subset <- List(lapply(covrle_list_by_strand, function(covrle) covrle[[chr]]))
    if (verbose & n_tiles > 1) cat("Tile Page:")
    lapply(indices, function(i) {
      if (verbose & n_tiles > 1) cat(" ", i, ", ", sep = "")
      if (!file.exists(file_paths[i]) | remake) {
        # TODO: resolutions go here ->
        write_fst(setnames(setDT(lapply(chr_subset, function(x) as.integer(x[irl[[i]]]))), colnames),
                  file_paths[i])
      }
      return(invisible(NULL))
    })
    if (verbose & n_tiles > 1) cat("\n")

    data.table(
      chr = chr,
      part = indices,
      start = start(unlist(irl)),
      end = end(unlist(irl)),
      file = file_paths
    )}))
}

#' Build or update the megafst (indexed coverage pages) and megacovrle of organism collections
#'
#' Incremental and resource-checked since 2026-10: libraries are streamed in batches from their
#' per-library covRle files, only new/changed libraries are appended/replaced on rerun (fst fork
#' column append/replace), full rebuilds are staged and published file by file, and every load
#' waits for available memory. Page layout, file names and values are unchanged; see
#' \code{megafst_update()} / \code{megacovrle_update()} and
#' ~/shared_workspace/hakon_playground/megafst for benchmarks.
#' @param all_exp,lib_type,organisms,min_samples as before
#' @param tile_size integer, page length (keep 2e6: larger pages slow cold reads, see notes)
#' @param force_remake_megacovrle logical, default FALSE (was TRUE): force a full megacovrle rebuild
#' @param force_remake_megafst logical, default FALSE: force a full page rebuild
#' @param make_megacovrle logical, default TRUE: also maintain the per-strand megacovrle files
#' @param fst_threads integer, fst threads per process
#' @param bp2 legacy chromosome interleave factor (page order)
#' @param BPPARAM BiocParallel param for page writing within a chromosome
#' @param batch_size integer, libraries loaded per batch
#' @param compress fst compression level
#' @param keep_free_gb,timeout_minutes RAM guard, see wait_for_memory()
#' @param ... passed to megafst_update() (e.g. compact_fraction, max_appends)
#' @return invisible(NULL)
fst_megaformat_make <- function(all_exp = list.experiments(validate = FALSE, pattern = "all_samples", libtypeExclusive = lib_type),
                                lib_type = "RFP",
                                organisms = unique(all_exp[!(organism %in% "all_species")][samples > min_samples][order(samples, decreasing = TRUE)]$name),
                                min_samples = 10,
                                tile_size = 2e6,
                                force_remake_megacovrle = FALSE,
                                force_remake_megafst = FALSE,
                                make_megacovrle = TRUE,
                                fst_threads = 5,
                                bp2 = 5, BPPARAM = BiocParallel::MulticoreParam(4, exportglobals = FALSE, log = FALSE),
                                batch_size = 500, compress = 50, keep_free_gb = 32, timeout_minutes = 120,
                                ...) {
  fst::threads_fst(fst_threads, reset_after_fork = FALSE)
  for (organism in organisms) {
    df <- read.experiment(organism, validate = FALSE)
    message("- ", organism(df))
    cov_paths <- try(filepath(df, "cov", suffix_stem = c("", "_pshifted"),
                              base_folders = libFolder(df, mode = "all")), silent = TRUE)
    if (is(cov_paths, "try-error") || !all(file.exists(cov_paths))) {
      message("-- Organism does not have cov RLEs for all libraries, skipping!")
      next
    }
    runs <- runIDs(df)
    if (make_megacovrle) {
      file_paths <- file.path(resFolder(df), "megalist_covrle",
                              paste0("all_covrles_", sub(" ", "_", organism(df)), c("_forward", "_reverse"), ".qs"))
      megacovrle_update(runs, cov_paths, file_paths, force = force_remake_megacovrle,
                        keep_free_gb = keep_free_gb, timeout_minutes = timeout_minutes,
                        collection_csv = file.path(config()["exp"], paste0(organism, ".csv")))
    }
    megafst_update(runs, cov_paths, seqinfo(df), file.path(resFolder(df), "collection_tables_indexed"),
                   tile_size = tile_size, batch_size = batch_size, compress = compress,
                   force = force_remake_megafst, keep_free_gb = keep_free_gb,
                   timeout_minutes = timeout_minutes, BPPARAM = BPPARAM, bp2 = bp2, ...)
  }
  message("Done")
  invisible(NULL)
}

#' Build or update the meta-motif tables of organism collections
#'
#' Cached and incremental since 2026-10, see \code{meta_motif_update()}: anchor windows are selected
#' once (seeded subsampling) and only new/changed libraries are scored on rerun.
#' @param all_exp,lib_type,organisms as before
#' @param force_remake_megacovrle logical, default FALSE (was TRUE): force recomputing regions and
#' every library (name kept for compatibility)
#' @param letters_vec character, amino acids whose ordered pairs are anchors
#' @param BPPARAM BiocParallel param over libraries
#' @param ... passed to meta_motif_update()
#' @return invisible(NULL)
meta_meta_motif_motif <- function(all_exp = list.experiments(validate = FALSE, pattern = "all_samples", libtypeExclusive = lib_type),
                                  lib_type = "RFP",
                                  organisms = unique(all_exp[!(organism %in% "all_species")][samples > 10][order(samples, decreasing = TRUE)]$name),
                                  force_remake_megacovrle = FALSE,
                                  letters_vec = c("P", "R"),
                                  BPPARAM = BiocParallel::MulticoreParam(8, exportglobals = FALSE, log = FALSE), ...) {
  for (org in organisms) {
    df <- read.experiment(org, validate = FALSE)
    message("- ", gsub(" ", "_", tolower(organism(df))))
    res <- try(meta_motif_update(df, letters_vec = letters_vec, force = force_remake_megacovrle,
                                 BPPARAM = BPPARAM, ...))
    if (is(res, "try-error")) warning("meta-motif failed for ", org, ": ", as.character(res))
  }
  invisible(NULL)
}
