#' Build the per-sample QC table used by RiboCrypt (metadata_done_samples_extended_qc.csv)
#'
#' Online Ribo-seq only. Replaces the interactive QC scripts (shared_projects/QC/scripts/00_master_combine_all.R and
#' its inputs) with one cached function for the collection processing driver. For every sample in
#' \code{metadata} (normally \code{<project>/metadata_done_samples_extended.csv}) it adds:
#' \itemize{
#'  \item STAR alignment log columns (study \code{LOGS_SINGLE|LOGS/00_STAR_LOG_table.csv});
#'  \item frame usage over all read lengths and the top read length (study \code{QC_STATS/Ribo_frames_all.csv});
#'  \item read length min/max/mode and percent at mode (score-weighted, \code{ofst/read_length_distribution});
#'  \item BAM prefix metrics (layout, soft-clipped percent) and coverage metrics on unique mappers (biotype /
#'   exonic / intronic / intergenic shares, active protein-coding genes, sense percent, strandedness);
#'  \item TIS/TTS/PP peak features from the organism collection metacoverage tables
#'   (\code{all_samples-<Organism>/meta_collection_tables});
#'  \item terminal nucleotide proportions when \code{terminal_bias_dir} holds per-run \code{*_5prime_bias.tsv /
#'   *_3prime_bias.tsv} files.
#' }
#' The expensive per-run part (BAM + coverage) is cached in \code{cache_dir} and recomputed only when the run's BAM or
#' coverage file changed (size or modification time), so a rerun after reprocessing a few studies is fast. Samples
#' without inputs get NA. Columns keep the names of the previous table.
#' @param config pipeline config (online Ribo-seq project)
#' @param metadata data.table of done samples
#' @param out_file output csv (written atomically; NULL to skip writing)
#' @param cache_dir per-run cache directory
#' @param terminal_bias_dir optional directory with per-study folders of terminal bias tsv files
#' @param organisms optional subset of ScientificName values (default all)
#' @param workers parallel workers for the per-run part
#' @return the table, invisibly
#' @export
metadata_qc_table <- function(config = pipeline_config(),
                              metadata = fread(file.path(config$project, "metadata_done_samples_extended.csv")),
                              out_file = file.path(config$project, "metadata_done_samples_extended_qc.csv"),
                              cache_dir = file.path(config$project, "qc_cache"),
                              terminal_bias_dir = "~/livemount/shared_projects/QC/data/terminal_bias_results_joanna",
                              organisms = NULL, workers = 16) {
  meta <- qc_prepare_metadata(metadata)
  if (!is.null(organisms)) meta <- meta[ScientificName %in% organisms]
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  exp_dir <- ORFik::config()["exp"]
  feats <- list()
  for (org in names(sort(table(meta$ScientificName)))) {
    m_org <- meta[ScientificName == org]
    message("- QC table: ", org, " (", nrow(m_org), " samples, ", uniqueN(m_org$study), " studies)")
    studies <- unique(m_org$study)
    studies <- studies[file.exists(file.path(exp_dir, paste0(studies, ".csv")))]
    study_tabs <- lapply(studies, function(s) tryCatch(qc_study_features(s, m_org[study == s]$Run),
      error = function(e) { message("-- ", s, ": ", conditionMessage(e)); NULL }))
    st <- rbindlist(study_tabs, fill = TRUE)
    heavy <- qc_heavy_features(st, org, cache_dir, workers)
    motif <- tryCatch(qc_motif_features(org), error = function(e) { message("-- motif features: ", conditionMessage(e)); NULL })
    if (nrow(st)) {
      st <- st[, !c("bam", "cov", "experiment"), with = FALSE]
      if (!is.null(heavy) && nrow(heavy)) st <- merge(st, heavy, by = "Run", all.x = TRUE)
      if (!is.null(motif) && nrow(motif)) st <- merge(st, motif, by = "Run", all.x = TRUE)
    }
    feats[[org]] <- st
  }
  feat <- rbindlist(feats, fill = TRUE)
  term <- qc_terminal_bias(terminal_bias_dir, meta$Run)
  if (!is.null(term) && nrow(term)) feat <- merge(feat, term, by = "Run", all = TRUE)
  feat <- feat[!duplicated(Run)]
  out <- merge(meta, feat, by = "Run", all.x = TRUE, sort = FALSE)
  out <- qc_order_and_round(out, out_file)
  stopifnot(nrow(out) == nrow(meta), !anyDuplicated(out$Run))
  if (!is.null(out_file)) {
    tmp <- paste0(out_file, ".tmp"); fwrite(out, tmp); stopifnot(file.rename(tmp, out_file))
    message("- Wrote ", out_file, ": ", nrow(out), " samples, ", ncol(out), " columns")
  }
  invisible(out)
}

#' Bump when the per-run (heavy) metric definitions change, to invalidate cached rows
#' (v2, 2026-10-15: gene biotypes were not read, so rRNA/lncRNA shares were 0)
#' @noRd
QC_CACHE_VERSION <- "v2"
QC_BAM_CACHE_VERSION <- "bam_v1"   # BAM metrics cached separately (genome-wide window sample)

#' @noRd
qc_prepare_metadata <- function(metadata) {
  m <- copy(as.data.table(metadata))
  drop <- intersect(c("CT", "TC", "not_unique", "name", "KEEP", "UNIQUE", "CHECKED", "InsertSize", "InsertDev"), names(m))
  if (length(drop)) m[, (drop) := NULL]
  m[, mode := "online"]
  m[, study := paste(study_accession, reference_folder_name(ScientificName), sep = "-")]
  stopifnot(!anyDuplicated(m$Run))
  m
}

#' Cheap per-study features + per-run file paths for the heavy part
#' @noRd
qc_study_features <- function(study, runs) {
  df <- read.experiment(study, validate = FALSE, output.env = new.env())
  rid <- ORFik::runIDs(df)
  base <- data.table(Run = rid, experiment = study, short = bamVarName(df, skip.libtype = TRUE),
                     bam = df$filepath)
  d <- df; uniqueMappers(d) <- TRUE
  base[, cov := tryCatch(filepath(d, "cov", base_folders = libFolder(d, mode = "all")), error = function(e) NA_character_)]
  base <- base[Run %in% runs]
  aligned <- dirname(df$filepath[1])
  # STAR logs
  star_f <- file.path(aligned, c("LOGS_SINGLE", "LOGS", "LOGS_PAIRED"), "00_STAR_LOG_table.csv")
  star_f <- star_f[file.exists(star_f)]
  if (length(star_f)) {
    star <- rbindlist(lapply(star_f, fread), fill = TRUE)
    star[, Run := fifelse(grepl("[SED]RR[0-9]+", sample), sub("^.*?([SED]RR[0-9]+).*$", "\\1", sample, perl = TRUE), NA_character_)]
    star <- star[!is.na(Run) & !duplicated(Run)]
    base <- merge(base, star, by = "Run", all.x = TRUE)
  }
  # frame usage (all read lengths) and top read length
  fr_f <- file.path(aligned, "QC_STATS", "Ribo_frames_all.csv")
  if (file.exists(fr_f)) {
    fr <- fread(fr_f)
    fr[, short := sub("^[A-Za-z0-9]+_", "", fraction)]
    key <- if (any(fr$fraction %in% base$short)) "fraction" else "short"
    fr[, Run := base$Run[match(get(key), base$short)]]
    fr <- fr[!is.na(Run)]
    if (nrow(fr)) {
      fb <- dcast(fr[, .(p = sum(percent)), by = .(Run, frame)], Run ~ frame, value.var = "p")
      setnames(fb, setdiff(names(fb), "Run"), paste0("Frame_usage_", setdiff(names(fb), "Run")))
      fb <- merge(fb, fr[, .(top_readlength = length[which.max(percent)]), by = Run], by = "Run")
      base <- merge(base, fb, by = "Run", all.x = TRUE)
    }
  }
  # read length distribution (score weighted, pre-shift ofst)
  rl_dir <- file.path(aligned, "ofst", "read_length_distribution")
  rl <- rbindlist(lapply(base$Run, function(r) {
    f <- file.path(rl_dir, paste0(r, ".csv"))
    if (!file.exists(f)) return(NULL)
    x <- fread(f)[count > 0]
    if (!nrow(x)) return(NULL)
    data.table(Run = r, Read_Length_Min = min(x$read_length), Read_Length_Max = max(x$read_length),
               Read_Length_Mode = x$read_length[which.max(x$count)],
               Pct_Reads_At_Mode_Len = 100 * max(x$count) / sum(x$count))
  }))
  if (nrow(rl)) base <- merge(base, rl, by = "Run", all.x = TRUE)
  base[, short := NULL]
  base
}

#' File signature for cache invalidation
#' @noRd
qc_file_sig <- function(f) {
  i <- file.info(f)
  paste(ifelse(is.na(i$size), "NA", i$size), ifelse(is.na(i$mtime), "NA", format(as.numeric(i$mtime), digits = 15)))
}

#' Heavy per-run features (coverage biotype metrics cached per run; BAM metrics via qc_bam_features)
#' @noRd
qc_heavy_features <- function(st, org, cache_dir, workers) {
  if (is.null(st) || !nrow(st)) return(NULL)
  runs <- st[, .(Run, experiment, bam, cov)]
  runs[, sig := paste(QC_CACHE_VERSION, qc_file_sig(bam), qc_file_sig(cov))]
  runs[, cache := file.path(cache_dir, experiment, paste0(Run, ".rds"))]
  cached <- lapply(runs$cache, function(f) if (file.exists(f)) tryCatch(readRDS(f), error = function(e) NULL) else NULL)
  ok <- vapply(seq_along(cached), function(i) !is.null(cached[[i]]) && identical(cached[[i]]$sig, runs$sig[i]), logical(1))
  todo <- runs[!ok]
  message("-- heavy QC: ", sum(ok), " cached, ", nrow(todo), " to compute")
  if (nrow(todo)) {
    bins <- qc_species_bins(todo$experiment[1])
    res <- BiocParallel::bplapply(seq_len(nrow(todo)), function(i) {
      r <- todo[i]
      row <- qc_cov_metrics(r$cov, bins)
      dir.create(dirname(r$cache), showWarnings = FALSE, recursive = TRUE)
      saveRDS(list(sig = r$sig, row = row), r$cache)
      row
    }, BPPARAM = BiocParallel::MulticoreParam(workers = min(workers, nrow(todo)), stop.on.error = FALSE))
    res <- res[vapply(res, is.data.frame, logical(1))]
    rm(bins); gc()
  } else res <- list()
  out <- rbindlist(c(lapply(cached[ok], function(x) x$row), res), fill = TRUE)
  drop <- intersect(c("Sequencing_Layout", "Soft_Clipped_Reads_Pct"), names(out))
  if (length(drop)) out[, (drop) := NULL]
  bam <- qc_bam_features(runs, cache_dir, workers)
  merge(out, bam, by = "Run", all = TRUE)
}

#' BAM metrics per run, cached separately from the coverage metrics (signature: QC_BAM_CACHE_VERSION + BAM file)
#' @noRd
qc_bam_features <- function(runs, cache_dir, workers) {
  runs <- copy(runs)
  runs[, bsig := paste(QC_BAM_CACHE_VERSION, qc_file_sig(bam))]
  runs[, bcache := file.path(cache_dir, experiment, paste0(Run, "_bam.rds"))]
  cached <- lapply(runs$bcache, function(f) if (file.exists(f)) tryCatch(readRDS(f), error = function(e) NULL) else NULL)
  ok <- vapply(seq_along(cached), function(i) !is.null(cached[[i]]) && identical(cached[[i]]$sig, runs$bsig[i]), logical(1))
  todo <- runs[!ok]
  message("-- BAM QC: ", sum(ok), " cached, ", nrow(todo), " to compute")
  res <- if (nrow(todo)) BiocParallel::bplapply(seq_len(nrow(todo)), function(i) {
    r <- todo[i]
    row <- qc_bam_metrics(r$bam, r$Run)
    dir.create(dirname(r$bcache), showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(sig = r$bsig, row = row), r$bcache)
    row
  }, BPPARAM = BiocParallel::MulticoreParam(workers = min(workers, nrow(todo)), stop.on.error = FALSE)) else list()
  res <- res[vapply(res, is.data.frame, logical(1))]
  rbindlist(c(lapply(cached[ok], function(x) x$row), res), fill = TRUE)
}

#' BAM metrics over the whole BAM, streamed in chunks (weights from collapsed read names "_x<count>").
#' A coordinate-sorted prefix or genome-uniform windows are not representative of the reads (most reads sit on few
#' highly expressed loci), so every record is counted; results are cached per BAM signature.
#' @noRd
qc_bam_metrics <- function(bam, run, chunk = 2e6) {
  empty <- data.table(Run = run, Sequencing_Layout = NA_character_, Soft_Clipped_Reads_Pct = NA_real_)
  if (is.na(bam) || !file.exists(bam)) return(empty)
  tryCatch({
    bf <- Rsamtools::BamFile(bam, yieldSize = chunk); open(bf); on.exit(close(bf))
    p <- Rsamtools::ScanBamParam(what = c("qname", "flag", "cigar"))
    tot <- 0; clip <- 0; paired <- FALSE
    repeat {
      x <- Rsamtools::scanBam(bf, param = p)[[1]]
      if (!length(x$flag)) break
      w <- suppressWarnings(as.numeric(sub(".*_x", "", x$qname))); w[is.na(w)] <- 1
      tot <- tot + sum(w); clip <- clip + sum(w[grepl("S", x$cigar, fixed = TRUE)])
      paired <- paired || any(bitwAnd(x$flag, 1L) > 0)
    }
    if (tot == 0) return(empty)
    data.table(Run = run, Sequencing_Layout = ifelse(paired, "Paired-End", "Single-End"),
               Soft_Clipped_Reads_Pct = 100 * clip / tot)
  }, error = function(e) empty)
}

#' Species annotation bins for coverage metrics (built once per organism)
#' @noRd
qc_species_bins <- function(study) {
  df <- read.experiment(study, validate = FALSE, output.env = new.env())
  txdb <- loadTxdb(df)
  exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
  gtf <- tryCatch(ORFik:::getGtfPathFromTxdb(txdb, stop.error = FALSE), error = function(e) NULL)
  biotype <- NULL
  if (!is.null(gtf) && file.exists(gtf)) {
    g <- tryCatch(rtracklayer::import(gtf, feature.type = "gene"), error = function(e) NULL)
    bcol <- intersect(c("gene_biotype", "gene_type"), names(S4Vectors::mcols(g)))[1]   # Ensembl / GENCODE naming
    if (!is.null(g) && !is.na(bcol))
      biotype <- setNames(as.character(S4Vectors::mcols(g)[[bcol]]), as.character(g$gene_id))
  }
  if (is.null(biotype)) {   # no biotypes in the annotation: genes with a CDS are protein coding
    message("-- no gene biotypes in the annotation; protein coding = genes with a CDS, rRNA/lncRNA shares not available")
    cds_genes <- unique(txNamesToGeneNames(names(loadRegion(txdb, "cds")), txdb))
    biotype <- setNames(ifelse(names(exons_by_gene) %in% cds_genes, "protein_coding", "other"), names(exons_by_gene))
  }
  ex <- unlist(exons_by_gene)
  bt <- biotype[names(ex)]; bt[is.na(bt)] <- "Other_Exonic"
  macro <- rep("Other_Exonic", length(bt))
  macro[bt == "protein_coding"] <- "Protein_Coding"
  macro[bt %in% c("rRNA", "Mt_rRNA", "rRNA_pseudogene")] <- "rRNA"
  macro[bt == "lncRNA"] <- "lncRNA"
  macro[grepl("pseudogene", bt, ignore.case = TRUE) & !bt %in% "rRNA_pseudogene"] <- "Pseudogene"
  red <- unlist(GenomicRanges::reduce(S4Vectors::splitAsList(unname(ex), macro)))
  dj <- GenomicRanges::disjoin(red)
  hits <- GenomicRanges::findOverlaps(dj, red)
  hd <- data.table(q = S4Vectors::queryHits(hits), m = names(red)[S4Vectors::subjectHits(hits)])[, .(n = uniqueN(m), m = m[1]), by = q]
  dj$category <- "Unknown"
  dj$category[hd[n == 1]$q] <- hd[n == 1]$m
  dj$category[hd[n > 1]$q] <- "Overlapping_Exonic"
  all_exons <- GenomicRanges::reduce(dj)
  introns <- tryCatch(GenomicRanges::setdiff(GenomicRanges::reduce(unlist(loadRegion(txdb, "intron"))), all_exons),
                      error = function(e) GenomicRanges::GRanges())
  flat <- function(region) tryCatch(GenomicRanges::reduce(unlist(loadRegion(txdb, region))), error = function(e) GenomicRanges::GRanges())
  anti <- all_exons; BiocGenerics::strand(anti) <- BiocGenerics::invertStrand(BiocGenerics::strand(all_exons))
  pc_ids <- names(biotype)[biotype == "protein_coding"]
  list(bins = S4Vectors::splitAsList(dj, dj$category), introns = introns, all_exons = all_exons, anti = anti,
       leaders = flat("leaders"), cds = flat("cds"), trailers = flat("trailers"),
       pc_genes = exons_by_gene[names(exons_by_gene) %in% pc_ids])
}

#' Coverage metrics from a covRle file (same definitions as the former QC script)
#' @noRd
qc_cov_metrics <- function(cov_file, b) {
  empty <- data.table(Run = NA_character_, Total_Cov_Signal = NA_real_)
  if (is.na(cov_file) || !file.exists(cov_file)) return(empty)
  tryCatch({
    cov <- qs2::qs_read(cov_file, nthreads = 1)
    cnt <- function(f) if (length(f) == 0) integer(0) else countOverlaps(f, cov)   # empty regions (e.g. no UTRs) error otherwise
    total <- sum(sum(cov))
    cat_c <- cnt(b$bins)
    g <- function(n) if (n %in% names(cat_c)) cat_c[[n]] else 0
    exonic <- sum(cat_c); intron <- sum(cnt(b$introns))
    intergenic <- max(0, total - exonic - intron)
    pc_gene <- cnt(b$pc_genes)
    anti <- sum(cnt(b$anti))
    sense <- 100 * exonic / max(1, exonic + anti)
    lead <- sum(cnt(b$leaders)); cds <- sum(cnt(b$cds)); trail <- sum(cnt(b$trailers)); mrna <- max(1, lead + cds + trail)
    pct <- function(x) 100 * x / max(1, total)
    data.table(Run = NA_character_, Total_Cov_Signal = total,
               ProteinCoding_Pct = pct(g("Protein_Coding")), rRNA_Pct = pct(g("rRNA")), lncRNA_Pct = pct(g("lncRNA")),
               Pseudogene_Pct = pct(g("Pseudogene")), OtherExonic_Pct = pct(g("Other_Exonic")),
               Overlapping_Exon_Pct = pct(g("Overlapping_Exonic")), Exonic_Pct = pct(exonic), Intron_Pct = pct(intron),
               Intergenic_Pct = pct(intergenic),
               Leader_Pct_mRNA = 100 * lead / mrna, CDS_Pct_mRNA = 100 * cds / mrna, Trailer_Pct_mRNA = 100 * trail / mrna,
               Active_PC_Genes_min10 = sum(pc_gene >= 10), Active_PC_Genes_min100 = sum(pc_gene >= 100),
               Median_Expressed_PC_Level = { v <- pc_gene[pc_gene > 0]; if (length(v)) stats::median(v) else 0 },
               Sense_Mapping_Pct = sense,
               Library_Strandedness = if (sense > 85) "Forward (Sense)" else if (sense < 15) "Reverse (Antisense)" else "Unstranded")
  }, error = function(e) empty)
}

#' TIS/TTS/PP peak features from the organism collection metacoverage tables (columns = Run ids)
#' @noRd
qc_motif_features <- function(org) {
  coll <- paste0("all_samples-", gsub(" ", "_", org))
  if (!file.exists(file.path(ORFik::config()["exp"], paste0(coll, ".csv")))) return(NULL)
  dir <- file.path(ORFik::config()["bam"], coll, "meta_collection_tables")
  f <- file.path(dir, paste0("meta_collection_table_", c("TIS", "TTS", "PP"), ".fst")); names(f) <- c("TIS", "TTS", "PP")
  if (!all(file.exists(f))) return(NULL)
  center <- 300; tight <- seq(-2, 2) + center; tight7 <- seq(-2, 4) + center; tight7a <- seq(-4, 2) + center
  long <- seq(1, 40); up <- long + center - 70; down <- long + center + 15; ratio <- length(long) / length(tight)
  TTS <- fst::read_fst(f["TTS"], as.data.table = TRUE); TIS <- fst::read_fst(f["TIS"], as.data.table = TRUE)
  PP <- fst::read_fst(f["PP"], as.data.table = TRUE)
  num_cols <- function(dt) names(dt)[vapply(dt, is.numeric, logical(1))]
  runs <- Reduce(intersect, list(num_cols(TIS), num_cols(TTS), num_cols(PP)))
  TIS <- TIS[, ..runs]; TTS <- TTS[, ..runs]; PP <- PP[, ..runs]
  cs <- function(dt, rows) colSums(dt[rows])
  drop_off <- cs(TTS, tight - 3) / pmax(1, cs(TTS, down) / ratio)
  drop_in <- cs(TIS, tight) / pmax(1, cs(TIS, up) / ratio)
  distr <- function(dt, rows, nm) { m <- as.matrix(dt[rows]); s <- sweep(m, 2L, colSums(m), "/"); s[!is.finite(s)] <- 0
    s <- as.data.table(t(s)); setnames(s, nm); s }
  data.table(Run = runs, drop_off = drop_off, drop_in = drop_in, drop_in_off_ratio = drop_in / drop_off,
             tis_peak_strength = cs(TIS, tight) / pmax(1, cs(TIS, down) / ratio),
             tts_peak_strength = cs(TTS, tight - 3) / pmax(1, cs(TTS, up) / ratio),
             pp_peak_strength = cs(PP, tight + 3) / pmax(1, cs(PP, down) / ratio),
             distr(TIS, tight7, paste0("TIS_distr_", tight7 - center)),
             distr(TTS, tight7a - 3, paste0("TTS_distr_A_site_", tight7a - center + 1)),
             tis_peak_strength_vs_15ntpluss = cs(TIS, tight) / pmax(1, cs(TIS, tight + 15)),
             cds_strength_vs_15ntpluss = cs(TIS, down + 4) / pmax(1, cs(TIS, tight + 15) / (1 / ratio)))
}

#' Terminal nucleotide proportions (first base at the read 5' / 3' end) from per-run bias tsv files
#' @noRd
qc_terminal_bias <- function(dir, runs) {
  if (is.null(dir) || !dir.exists(dir)) return(NULL)
  files <- list.files(dir, "_[53]prime_bias\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(NULL)
  info <- data.table(f = files, flank = ifelse(grepl("_5prime_bias", files), "5'", "3'"),
                     Run = sub("_[53]prime_bias\\.tsv$", "", basename(files)))
  info <- info[Run %in% runs]
  if (!nrow(info)) return(NULL)
  obs <- rbindlist(lapply(seq_len(nrow(info)), function(i) {
    dt <- fread(info$f[i]); col <- if (info$flank[i] == "5'") "read_5prime_3nt" else "read_3prime_3nt"
    dt[, base := substr(get(col), 1L, 1L)][, .(o = sum(count)), by = base][, .(Run = info$Run[i], key = paste0(base, info$flank[i]), p = o / sum(o))]
  }))
  dcast(obs[!duplicated(paste(Run, key))], Run ~ key, value.var = "p")
}

#' Column order of the previous table, rounded numerics
#' @noRd
qc_order_and_round <- function(out, out_file) {
  if (!is.null(out_file) && file.exists(out_file)) {
    prev <- names(fread(out_file, nrows = 0))
    setcolorder(out, c(intersect(prev, names(out)), setdiff(names(out), prev)))
  }
  num <- names(out)[vapply(out, function(x) is.numeric(x) && !is.integer(x), logical(1))]
  if (length(num)) out[, (num) := lapply(.SD, round, 2), .SDcols = num]
  out
}
