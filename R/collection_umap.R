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
#' The embedding is cached: when the full count-matrix hash, sample and gene
#' identities, parameters and runtime versions are unchanged, only the metadata columns are refreshed.
#' @param all_exp data.table of experiments (list.experiments)
#' @param lib_type character, library type of the collection (RFP, RNA, ...)
#' @param organisms collection experiment names
#' @param meta_rc data.table of curated metadata with a Run column
#' @param min_samples minimum number of samples to make a UMAP
#' @param ntop number of most variable genes used for the PCA
#' @param n_neighbors,min_dist uwot parameters (n_neighbors is capped at samples - 1)
#' @param seed integer, makes the embedding reproducible
#' @param n_components integer vector, embedding dimensionalities (e.g. \code{3:6}); each is
#'   an independent fit from the same counts. Values other than 2 are saved
#'   separately (for example \code{UMAP_by_gene_counts_6d.fst}). UMAP axes from
#'   independent embeddings are not interchangeable.
#' @param threads threads for the uwot neighbour search
#' @param force logical, recompute even if the cached embedding is valid
#' @return invisible(NULL), files are written
UMAP_by_gene_counts <- function(all_exp = list.experiments(validate = FALSE, pattern = "all_samples", libtypeExclusive = lib_type),
                                lib_type = "RFP",
                                organisms = unique(all_exp[!(organism %in% "all_species")]$name),
                                meta_rc = fread("~/livemount/Bio_data/NGS_pipeline/FINAL_LIST.csv"),
                                min_samples = 10, ntop = 500, n_neighbors = 30, min_dist = 2,
                                seed = 1L, threads = 8, force = FALSE,
                                n_components = 2L) {
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is needed for UMAP_by_gene_counts")
  if (!length(n_components) || anyNA(n_components) || any(n_components < 2) ||
      any(n_components != as.integer(n_components)) || anyDuplicated(n_components))
    stop("n_components must be distinct integers of at least 2")
  n_components <- as.integer(n_components)
  for (org in unique(organisms)) {
    df_all <- read.experiment(org, validate = FALSE)
    message("- ", organism(df_all))
    if (nrow(df_all) < min_samples) {
      message("-- Too few samples to make UMAP, must be at least ", min_samples, "!")
      next
    }
    dir <- file.path(refFolder(df_all), "UMAP")
    type <- ifelse(lib_type == "RFP", "", paste0("_", lib_type))
    srrs <- ORFik::runIDs(df_all)
    counts <- umap_canonical_counts(df_all)
    counts_hash <- digest::digest(counts, algo="xxhash64")
    for (n_dim in n_components) {
      dim_suffix <- if (n_dim == 2L) "" else paste0("_", n_dim, "d")
      path <- file.path(dir, paste0("UMAP_by_gene_counts", type, dim_suffix, ".fst"))
      cache_path <- sub("\\.fst$", "_inputs.rds", path)
      params <- list(ntop = ntop, n_neighbors = n_neighbors, min_dist = min_dist, seed = seed,
                     method = "log2-libsize_exactPCA_uwot_v1")
      params$n_components <- n_dim
      fingerprint <- list(samples = srrs, genes = rownames(counts),
                          totals = colSums(counts, na.rm = TRUE), params = params)
      # Equal per-sample totals do not prove that gene counts are unchanged.
      fingerprint$counts_hash <- counts_hash
      fingerprint$runtime <- list(R=as.character(getRversion()), uwot=as.character(utils::packageVersion("uwot")),
                                  threads=threads, n_sgd_threads=0L)
      old <- if (!force && file.exists(path) && file.exists(cache_path)) readRDS(cache_path)
      if (identical(old, fingerprint)) {
        message("-- Counts and parameters unchanged, reusing embedding (metadata refreshed)")
        dt_umap <- fst::read_fst(path, columns = paste("UMAP",seq_len(n_dim)), as.data.table = TRUE)
      } else {
        dt_umap <- umap_from_counts(counts, ntop = ntop, n_neighbors = n_neighbors,
                                    min_dist = min_dist, seed = seed, threads = threads,
                                    n_components = n_dim)
      }

      dt_umap <- cbind(dt_umap, umap_sample_metadata(srrs, meta_rc))
      dir.create(dir, FALSE, TRUE)
      tmp <- tempfile(pattern=".umap-fst-",tmpdir=dir)
      fst::write_fst(dt_umap, tmp)
      if (!file.rename(tmp, path)) stop("Unable to publish UMAP file: ",path)
      tmp_inputs <- tempfile(pattern=".umap-inputs-",tmpdir=dir)
      saveRDS(fingerprint, tmp_inputs)
      if (!file.rename(tmp_inputs,cache_path)) stop("Unable to publish UMAP inputs: ",cache_path)
      message("Done , saved UMAP fst to location: ", path)
    }
    rm(counts)
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
                             threads = 8, chunk_rows = 2000, n_components = 2L) {
  n <- ncol(counts)
  if (length(n_components) != 1L || is.na(n_components) || n_components < 2 ||
      n_components != as.integer(n_components) || n_components >= n)
    stop("n_components must be an integer between 2 and the number of samples minus one")
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
  emb <- uwot::umap(pcs, n_neighbors = min(n_neighbors, n - 1), min_dist = min_dist, n_components = n_components,
                    pca = min(50, dim(pcs)), n_threads = threads, n_sgd_threads = 0,
                    seed = seed, verbose = FALSE)
  dt <- as.data.table(emb)
  colnames(dt) <- paste("UMAP",seq_len(n_components))
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

#' Describe UMAP axes by marker-gene programs
#'
#' Writes \code{<refFolder>/UMAP/UMAP_<d>d_axis_annotations.csv} (columns \code{axis, label, source_md5},
#' read by the RiboCrypt Observatory axis titles) and \code{UMAP_<d>d_axis_associations.csv}
#' (all module x axis correlations) for each embedding \code{UMAP_by_gene_counts[_<d>d].fst}.\cr
#' Module scores are computed from the same canonical-isoform counts as the embedding, so every
#' sample is scored: log2(count / size factor + 1), minus the per-sample mean of housekeeping genes,
#' then mean of the up genes minus mean of the down genes (\code{umap_marker_modules()}).
#' Each axis is described by the two modules with the largest absolute within-study correlation
#' (scores and coordinates centred per BioProject, studies with at least \code{min_study_samples}),
#' with the sign of the correlation; associations weaker than \code{weak} are marked "weak".
#' These are exploratory associations of one fitted embedding, not loadings; rerun after every refit
#' (the label file carries the md5 of the embedding it describes, RiboCrypt ignores stale labels).
#' @param df_all the collection experiment (all_samples-<organism>)
#' @param n_components integer vector, the embeddings to describe (2 = UMAP_by_gene_counts.fst)
#' @param modules data.table with columns module, up, down (lists of gene symbols)
#' @param min_study_samples studies with fewer samples are left out of the within-study correlation
#' @param weak absolute correlation below which an association is marked weak
#' @param min_gene_fraction minimum fraction of module genes present in the counts
#' @return invisible data.table of all associations
#' @export
UMAP_axis_annotations <- function(df_all, n_components = 2:6, modules = umap_marker_modules(),
                                  min_study_samples = 3, weak = 0.2, min_gene_fraction = 0.5) {
  dir <- file.path(refFolder(df_all), "UMAP")
  counts <- umap_canonical_counts(df_all)
  colnames(counts) <- ORFik::runIDs(df_all)   # same order, see UMAP_by_gene_counts
  sym <- symbols(df_all)
  genes <- sym$external_gene_name[match(rownames(counts), sym$ensembl_tx_name)]
  scores <- umap_module_scores(counts, genes, modules, min_gene_fraction)
  res <- list()
  for (n_dim in as.integer(n_components)) {
    path <- file.path(dir, paste0("UMAP_by_gene_counts", if (n_dim == 2L) "" else paste0("_", n_dim, "d"), ".fst"))
    if (!file.exists(path)) {
      message("-- No embedding ", path, ", skipped")
      next
    }
    source_md5 <- unname(tools::md5sum(path))
    u <- fst::read_fst(path, as.data.table = TRUE)
    if (!setequal(u$sample, colnames(counts)))
      stop("Embedding ", basename(path), " does not contain the collection samples; refit it first")
    assoc <- umap_axis_associations(u, scores, n_dim, min_study_samples)
    labels <- assoc[, .(label = umap_axis_label(module, within_study, weak)), by = axis][order(axis)]
    for (out in c("annotations", "associations")) {
      file <- file.path(dir, paste0("UMAP_", n_dim, "d_axis_", out, ".csv"))
      tmp <- tempfile(pattern = ".umap-axis-", tmpdir = dir)
      if (out == "annotations") {
        fwrite(data.table(axis = labels$axis, label = labels$label, source_md5 = source_md5), tmp)
      } else fwrite(cbind(dim = n_dim, assoc, source_md5 = source_md5), tmp)
      if (!file.rename(tmp, file)) stop("Unable to publish ", file)
    }
    message("-- ", n_dim, "D axes: ", paste(labels$axis, labels$label, sep = ": ", collapse = " | "))
    res[[length(res) + 1L]] <- cbind(dim = n_dim, assoc)
  }
  invisible(rbindlist(res))
}

#' Marker-gene modules used to describe UMAP axes (human symbols)
#'
#' Short-named programs from the dominant cell state panel (forks/dominant_cell_states):
#' proliferation, TOP translation, ISR, interferon, hypoxia, EMT, OXPHOS, the acute SASP core
#' and senescence commitment (CDKN2A up, LMNB1 down). Attribute \code{baseline} holds the
#' housekeeping genes used as per-sample reference.
#' @return data.table with columns module, up, down
#' @export
umap_marker_modules <- function() {
  m <- data.table(
    module = c("proliferation", "TOP translation", "ISR markers", "interferon", "hypoxia", "EMT",
               "OXPHOS", "acute SASP", "CDKN2A/LMNB1"),
    up = list(
      c("MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDK1", "PLK1", "RRM1",
        "RRM2", "FEN1", "PRIM2", "POLA1"),
      c("RPS6", "RPLP0", "RPL32", "RPL13A", "RPS3", "EIF3E", "EIF4B", "EEF2", "PABPC1", "NPM1",
        "NCL", "FBL", "MYC"),
      c("ATF4", "DDIT3", "ATF3", "ASNS", "PPP1R15A", "DDIT4", "CEBPB", "HERPUD1", "HSPA5", "XBP1"),
      c("STAT1", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IFI44", "CXCL10"),
      c("CA9", "VEGFA", "EGLN3", "SLC2A1", "BNIP3", "LDHA", "PDK1", "NDRG1"),
      c("VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "CDH2", "MMP2", "THY1", "POSTN"),
      c("NDUFA1", "NDUFS1", "SDHB", "UQCRC1", "COX5A", "ATP5F1A", "ATP5MC1", "PHB2"),
      c("SERPINE1", "IGFBP7", "TGFBI", "TIMP1", "CCND1", "IL6", "CDKN1A", "CXCL8"),
      "CDKN2A"),
    down = list(character(), character(), character(), character(), character(),
                c("CDH1", "EPCAM", "KRT8", "KRT18", "KRT19"), "SLC2A1", character(), "LMNB1"))
  attr(m, "baseline") <- c("ACTB", "GAPDH", "HPRT1", "TBP", "PPIA", "YWHAZ", "EEF1A1", "PSMB4",
                           "RPS18", "RPL27", "RPL30", "RPS27A", "TPT1", "RPL13A")
  m
}

#' Module scores per sample (samples x modules) from canonical counts, genes = symbol per count row
#' @noRd
umap_module_scores <- function(counts, genes, modules, min_gene_fraction = 0.5) {
  gene <- toupper(genes)
  baseline <- attr(modules, "baseline")
  wanted <- unique(c(unlist(modules$up), unlist(modules$down), baseline))
  rows <- which(gene %in% wanted)
  sf <- colSums(counts, na.rm = TRUE)
  sf <- sf / mean(sf)
  x <- counts[rows, , drop = FALSE]
  storage.mode(x) <- "double"
  x <- rowsum(x, gene[rows])           # one row per symbol
  x <- log2(sweep(x, 2, sf, `/`) + 1)
  base <- intersect(baseline, rownames(x))
  if (length(base) < 5) stop("Too few housekeeping genes found in the counts (", length(base), ")")
  x <- sweep(x, 2, colMeans(x[base, , drop = FALSE]), `-`)
  mean_of <- function(g) if (length(g)) colMeans(x[g, , drop = FALSE]) else 0
  scores <- list()
  for (i in seq_len(nrow(modules))) {
    up <- intersect(modules$up[[i]], rownames(x)); down <- intersect(modules$down[[i]], rownames(x))
    total <- length(modules$up[[i]]) + length(modules$down[[i]])
    if ((length(up) + length(down)) / total < min_gene_fraction || (!length(up) && !length(down))) {
      message("-- Module ", modules$module[i], " skipped (too few genes present)")
      next
    }
    scores[[modules$module[i]]] <- mean_of(up) - mean_of(down)
  }
  if (!length(scores)) stop("No marker module could be scored")
  data.table(sample = colnames(counts), as.data.table(scores))
}

#' Pearson and within-study correlation of each module with each axis
#' @noRd
umap_axis_associations <- function(u, scores, n_dim, min_study_samples = 3) {
  z <- merge(u[, c("sample", "BioProject", paste("UMAP", seq_len(n_dim))), with = FALSE], scores, by = "sample")
  if (nrow(z) != nrow(u)) stop(nrow(u) - nrow(z), " embedding samples have no module scores")
  modules <- setdiff(names(scores), "sample")
  study <- z$BioProject
  ok_study <- !is.na(study) & nzchar(study)
  n_study <- table(study[ok_study])
  use <- ok_study & as.vector(n_study[study]) >= min_study_samples & !is.na(as.vector(n_study[study]))
  centre <- function(v) v[use] - ave(v[use], study[use], FUN = mean)
  rbindlist(lapply(seq_len(n_dim), function(j) {
    a <- z[[paste("UMAP", j)]]
    rbindlist(lapply(modules, function(m) {
      v <- z[[m]]
      data.table(axis = j, module = m, n = length(v), n_within = sum(use),
                 pearson = stats::cor(v, a), within_study = stats::cor(centre(v), centre(a)))
    }))
  }))
}

#' Axis label: the two strongest modules by absolute within-study correlation
#' @noRd
umap_axis_label <- function(module, r, weak = 0.2) {
  top <- head(order(-abs(r)), 2)
  paste(sprintf("%s (%s%s)", module[top], ifelse(r[top] < 0, "-", "+"),
                ifelse(abs(r[top]) < weak, "; weak", "")), collapse = ", ")
}
