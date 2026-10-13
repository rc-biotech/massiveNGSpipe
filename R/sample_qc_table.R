#' One QC row per sample of a processed study
#'
#' Collects the per-sample statistics that are otherwise spread over a study's files into one Run-keyed table,
#' so cross-study outlier detection (wrong adapter/barcode/UMI, failed P-shifting, non-Ribo-seq libraries) is
#' cheap: fastp reports (raw/trimmed reads, adapter-trimmed %, trimmed length), the adapter/barcode table,
#' STAR logs (unique / multi / too-short %; contaminant % from the rRNA/tRNA/ncRNA depletion step), aligned (pre-shift ofst) and P-shifted read-length distributions
#' (score-weighted: median, mode, % inside \code{accepted_lengths}), the frame table (frame-0 % over all
#' lengths and over accepted lengths) and the P-site offsets used. Missing sources give NA (older runs).
#' @param exp experiment name
#' @param accepted_lengths read lengths the pipeline P-shifts (config$accepted_lengths_rpf)
#' @param exp_dir directory of experiment CSVs
#' @return data.table, one row per Run
#' @export
sample_qc_table <- function(exp, accepted_lengths = c(20, 21, 25:33), exp_dir = ORFik::config()["exp"]) {
  df <- read.experiment(file.path(exp_dir, paste0(exp, ".csv")), validate = FALSE, output.env = new.env())
  runs <- ORFik::runIDs(df)
  al <- libFolder(df)
  study_dir <- dirname(al)
  out <- data.table(experiment = exp, Run = runs, label = bamVarName(df))
  num <- function(x) suppressWarnings(as.numeric(x))
  # fastp
  fp <- rbindlist(lapply(runs, function(r) {
    f <- file.path(study_dir, "trim", paste0("report_", r, ".json"))
    if (!file.exists(f)) return(NULL)
    j <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
    if (is.null(j)) return(NULL)
    s <- j$summary
    data.table(Run = r, raw_reads = s$before_filtering$total_reads, trimmed_reads = s$after_filtering$total_reads,
               raw_mean_length = s$before_filtering$read1_mean_length, trimmed_mean_length = s$after_filtering$read1_mean_length,
               adapter_trimmed_pct = if (!is.null(j$adapter_cutting$adapter_trimmed_reads))
                 100 * j$adapter_cutting$adapter_trimmed_reads / max(1, s$before_filtering$total_reads) else NA_real_)
  }), fill = TRUE)
  if (nrow(fp)) out <- merge(out, fp, by = "Run", all.x = TRUE, sort = FALSE)
  # adapter / barcode table
  ab <- file.path(study_dir, "trim", "adapter_barcode_table.csv")
  if (file.exists(ab)) {
    x <- fread(ab, colClasses = "character")
    x <- x[, .(Run = id, adapter, barcode_detected, barcode5p_size = num(barcode5p_size), barcode3p_size = num(barcode3p_size),
               no_adapter_removed_pct = num(`reads_no_adapter_removed fastp(%)`))]
    out <- merge(out, x, by = "Run", all.x = TRUE, sort = FALSE)
  }
  # STAR logs (collapsed reads): genome (aligned/LOGS, LOGS_SINGLE, LOGS_PAIRED by pipeline era) and contaminants
  star_log <- function(dirs, r) {
    for (d in dirs) {
      f <- list.files(d, paste0(r, "_Log\\.final\\.out$"), full.names = TRUE)[1]
      if (!is.na(f)) return(readLines(f))
    }
    NULL
  }
  g <- function(l, key) num(sub("%", "", trimws(sub(".*\\|", "", grep(key, l, value = TRUE, fixed = TRUE)[1]))))
  genome_dirs <- list.dirs(al, recursive = FALSE)[grepl("^LOGS", basename(list.dirs(al, recursive = FALSE)))]
  contam_dirs <- list.dirs(file.path(study_dir, "contaminants_depletion"), recursive = FALSE)
  contam_dirs <- contam_dirs[grepl("^LOGS", basename(contam_dirs))]
  st <- rbindlist(lapply(runs, function(r) {
    l <- star_log(genome_dirs, r); cl <- star_log(contam_dirs, r)
    if (is.null(l) && is.null(cl)) return(NULL)
    data.table(Run = r,
               star_input = if (!is.null(l)) g(l, "Number of input reads") else NA_real_,
               star_unique_pct = if (!is.null(l)) g(l, "Uniquely mapped reads %") else NA_real_,
               star_multi_pct = if (!is.null(l)) g(l, "% of reads mapped to multiple loci") else NA_real_,
               star_too_short_pct = if (!is.null(l)) g(l, "% of reads unmapped: too short") else NA_real_,
               contaminant_input = if (!is.null(cl)) g(cl, "Number of input reads") else NA_real_,
               contaminant_pct = if (!is.null(cl)) g(cl, "Uniquely mapped reads %") + g(cl, "% of reads mapped to multiple loci") +
                 g(cl, "% of reads mapped to too many loci") else NA_real_)
  }), fill = TRUE)
  if (nrow(st)) out <- merge(out, st, by = "Run", all.x = TRUE, sort = FALSE)
  # read length distributions (score-weighted), aligned = ofst, shifted = pshifted
  rl <- function(dir, prefix) rbindlist(lapply(runs, function(r) {
    f <- file.path(al, dir, "read_length_distribution", paste0(r, ".csv"))
    if (!file.exists(f)) return(NULL)
    x <- fread(f)
    if (!nrow(x) || sum(x$count) == 0) return(NULL)
    x <- x[order(read_length)]; cs <- cumsum(x$count) / sum(x$count)
    res <- data.table(Run = r, x$read_length[which(cs >= 0.5)[1]], x$read_length[which.max(x$count)],
                      100 * sum(x$count[x$read_length %in% accepted_lengths]) / sum(x$count))
    setnames(res, c("Run", paste0(prefix, c("_median_length", "_mode_length", "_pct_in_accepted"))))
    res
  }), fill = TRUE)
  for (p in list(c("ofst", "aligned"), c("pshifted", "pshifted"))) {
    x <- rl(p[1], p[2])
    if (nrow(x)) out <- merge(out, x, by = "Run", all.x = TRUE, sort = FALSE)
  }
  # frame table (keyed by label without libtype)
  fr <- file.path(QCfolder(df), "Ribo_frames_all.csv")
  if (file.exists(fr)) {
    x <- fread(fr)
    short <- bamVarName(df, skip.libtype = TRUE)
    x[, Run := runs[match(fraction, short)]]
    x <- x[!is.na(Run)]
    f0 <- x[, .(frame0_pct = 100 * sum(score[frame == 0]) / sum(score),
                frame0_pct_accepted = 100 * sum(score[frame == 0 & length %in% accepted_lengths]) /
                  max(1, sum(score[length %in% accepted_lengths]))), by = Run]
    out <- merge(out, f0, by = "Run", all.x = TRUE, sort = FALSE)
  }
  # P-site offsets used
  sh <- file.path(al, "pshifted", "shifting_table.rds")
  if (file.exists(sh)) {
    s <- tryCatch(readRDS(sh), error = function(e) NULL)
    if (is.list(s) && length(s) == length(runs)) {
      out <- merge(out, data.table(Run = runs,
        shift_lengths = vapply(s, function(t) nrow(as.data.frame(t)), numeric(1)),
        shift_offsets = vapply(s, function(t) { t <- as.data.frame(t); paste0(t[[1]], ":", t[[2]], collapse = ",") }, character(1))),
        by = "Run", all.x = TRUE, sort = FALSE)
    }
  }
  out[match(runs, Run)]
}

#' Sample QC tables for many experiments, robust to individual failures
#' @param exps experiment names
#' @param ... passed to \code{sample_qc_table()}
#' @param BPPARAM BiocParallel parameters
#' @return data.table of all samples, plus attribute "failed" (experiment, error)
#' @export
sample_qc_tables <- function(exps, ..., BPPARAM = BiocParallel::SerialParam()) {
  res <- BiocParallel::bplapply(exps, function(e) tryCatch(sample_qc_table(e, ...), error = function(err) conditionMessage(err)),
                                BPPARAM = BPPARAM)
  failed <- vapply(res, is.character, logical(1))
  out <- rbindlist(res[!failed], fill = TRUE)
  setattr(out, "failed", data.table(experiment = exps[failed], error = unlist(res[failed])))
  out
}
