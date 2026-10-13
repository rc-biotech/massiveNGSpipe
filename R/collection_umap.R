#' UMAP of collection samples by canonical gene counts
#'
#' One UMAP per organism collection (\code{all_samples-<organism>}), saved to
#' \code{<refFolder>/UMAP/UMAP_by_gene_counts[_<lib_type>].fst} for the RiboCrypt UMAP page.\cr
#' Same method as the earlier scater pipeline (logNormCounts -> runPCA -> runUMAP), implemented
#' directly so that only \code{uwot} is needed (no scater/SingleCellExperiment, no matrixStats pin):
#' \enumerate{
#'  \item library size factors (column sums / mean column sum), log2(count / size factor + 1)
#'  \item PCA (exact, centred, not scaled) on the \code{ntop} most variable genes
#'  \item \code{uwot::umap} on the PCA scores
#' }
#' Log counts are computed in row chunks, only the top variable genes are materialised.
#' The embedding is cached: when the count fingerprint (samples, genes, per-sample totals) and
#' parameters are unchanged, only the metadata columns are refreshed.
#' @param all_exp data.table of experiments (list.experiments)
#' @param lib_type character, library type of the collection (RFP, RNA, ...)
#' @param organisms collection experiment names
#' @param meta_rc data.table of curated metadata with a Run column
#' @param min_samples minimum number of samples to make a UMAP
#' @param ntop number of most variable genes used for the PCA
#' @param n_neighbors,min_dist uwot parameters (n_neighbors is capped at samples - 1)
#' @param seed integer, makes the embedding reproducible
#' @param threads threads for the uwot neighbour search
#' @param force logical, recompute even if the cached embedding is valid
#' @return invisible(NULL), files are written
UMAP_by_gene_counts <- function(all_exp = list.experiments(validate = FALSE, pattern = "all_samples", libtypeExclusive = lib_type),
                                lib_type = "RFP",
                                organisms = unique(all_exp[!(organism %in% "all_species")]$name),
                                meta_rc = fread("~/livemount/Bio_data/NGS_pipeline/FINAL_LIST.csv"),
                                min_samples = 10, ntop = 500, n_neighbors = 30, min_dist = 2,
                                seed = 1L, threads = 8, force = FALSE) {
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is needed for UMAP_by_gene_counts")
  for (org in unique(organisms)) {
    df_all <- read.experiment(org, validate = FALSE)
    message("- ", organism(df_all))
    if (nrow(df_all) < min_samples) {
      message("-- Too few samples to make UMAP, must be at least ", min_samples, "!")
      next
    }
    dir <- file.path(refFolder(df_all), "UMAP")
    type <- ifelse(lib_type == "RFP", "", paste0("_", lib_type))
    path <- file.path(dir, paste0("UMAP_by_gene_counts", type, ".fst"))
    cache_path <- sub("\\.fst$", "_inputs.rds", path)

    srrs <- ORFik::runIDs(df_all)
    counts <- umap_canonical_counts(df_all)
    params <- list(ntop = ntop, n_neighbors = n_neighbors, min_dist = min_dist, seed = seed,
                   method = "log2-libsize_exactPCA_uwot_v1")
    fingerprint <- list(samples = srrs, genes = rownames(counts),
                        totals = colSums(counts, na.rm = TRUE), params = params)
    old <- if (!force && file.exists(path) && file.exists(cache_path)) readRDS(cache_path)
    if (identical(old, fingerprint)) {
      message("-- Counts and parameters unchanged, reusing embedding (metadata refreshed)")
      dt_umap <- fst::read_fst(path, columns = c("UMAP 1", "UMAP 2"), as.data.table = TRUE)
    } else {
      dt_umap <- umap_from_counts(counts, ntop = ntop, n_neighbors = n_neighbors,
                                  min_dist = min_dist, seed = seed, threads = threads)
    }
    rm(counts)

    dt_umap <- cbind(dt_umap, umap_sample_metadata(srrs, meta_rc))
    dir.create(dir, FALSE, TRUE)
    tmp <- paste0(path, ".tmp")
    fst::write_fst(dt_umap, tmp)
    file.rename(tmp, path)
    saveRDS(fingerprint, cache_path)
    message("Done , saved UMAP fst to location: ", path)
  }
  return(invisible(NULL))
}

#' Canonical-isoform count matrix (genes x samples) of a collection
#' @noRd
umap_canonical_counts <- function(df_all) {
  count_table <- countTable(df_all, type = "summarized")
  keep <- rownames(count_table) %in% canonical_isoforms(df_all)
  counts <- SummarizedExperiment::assay(count_table)[keep, , drop = FALSE]
  if (ncol(counts) != nrow(df_all))
    stop("Count table has ", ncol(counts), " samples, experiment has ", nrow(df_all))
  counts
}

#' Log-normalise, PCA on top variable genes, UMAP
#'
#' Matches scater::logNormCounts (library size factors centred to mean 1, log2(x/sf + 1)),
#' scater::runPCA (ntop by variance, centred, unscaled, exact SVD) and scater::runUMAP on the
#' PCA scores. Log counts are formed per row chunk; only the ntop rows are kept.
#' @noRd
umap_from_counts <- function(counts, ntop = 500, n_neighbors = 30, min_dist = 2, seed = 1L,
                             threads = 8, chunk_rows = 2000) {
  n <- ncol(counts)
  lib_sizes <- colSums(counts, na.rm = TRUE)
  if (any(lib_sizes <= 0))
    stop("Samples with zero counts in canonical genes: ", paste(colnames(counts)[lib_sizes <= 0], collapse = ", "))
  sf <- lib_sizes / mean(lib_sizes)
  inv_sf <- 1 / sf
  log_rows <- function(rows) {
    x <- counts[rows, , drop = FALSE]
    storage.mode(x) <- "double"
    log2(sweep(x, 2, inv_sf, `*`) + 1)
  }
  chunks <- split(seq_len(nrow(counts)), ceiling(seq_len(nrow(counts)) / chunk_rows))
  vars <- unlist(lapply(chunks, function(rows) {
    x <- log_rows(rows)
    m <- rowMeans(x)
    rowSums((x - m)^2) / (n - 1)
  }), use.names = FALSE)
  top <- order(vars, decreasing = TRUE)[seq_len(min(ntop, length(vars)))]
  x <- t(log_rows(top))          # samples x genes
  x <- scale(x, center = TRUE, scale = FALSE)
  ncomp <- min(50, n - 1, ncol(x))
  s <- svd(x, nu = ncomp, nv = 0)
  pcs <- sweep(s$u, 2, s$d[seq_len(ncomp)], `*`)
  # pca = min(50, dim): what scater::runUMAP(dimred = "PCA") passed to uwot (re-centring rotation)
  emb <- uwot::umap(pcs, n_neighbors = min(n_neighbors, n - 1), min_dist = min_dist,
                    pca = min(50, dim(pcs)), n_threads = threads, n_sgd_threads = 0,
                    seed = seed, verbose = FALSE)
  dt <- as.data.table(emb)
  colnames(dt) <- c("UMAP 1", "UMAP 2")
  dt
}

#' Metadata columns of the UMAP table, in sample order
#' @noRd
umap_sample_metadata <- function(srrs, meta_rc) {
  m <- meta_rc[match(srrs, meta_rc$Run), ]
  if (all(is.na(m$Run))) stop("No samples of the collection found in the metadata")
  if (anyNA(m$Run)) warning(sum(is.na(m$Run)), " samples missing in the metadata")
  fix <- function(x) { x[is.na(x) | x == ""] <- "Unknown"; x }
  author <- m$AUTHOR; author[author %in% "Makar"] <- ""; author <- fix(author)
  cell_lines <- fix(m$CELL_LINE); tissues <- fix(m$TISSUE); inhibitors <- fix(m$INHIBITOR)
  data.table(sample = m$Run,
             cell_line = tolower(cell_lines),
             tissue = tolower(tissues),
             tissues_cell_lines = tolower(paste0(tissues, " | ", cell_lines)),
             inhibitors = tolower(inhibitors),
             BioProject = m$BioProject,
             author = author,
             studies = m$BioProject)
}
