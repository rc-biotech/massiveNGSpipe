# P-shift failure-cause diagnostics.
#
# Nothing in the pipeline today classifies WHY a P-shift is bad -- only a
# binary good/warning periodicity verdict exists (periodicity_check_flag(),
# R/shifting_helpers.R), and nothing downstream ever reads it back. This
# file adds a single, progressively-enriched per-experiment record that
# each pipeline stage writes its own piece of as soon as the relevant
# signal is available (raw read count after trim, alignment rate after
# align, periodicity after pshift), plus a classifier that turns those
# signals into one likely-cause label for triage.

#' bam_dir (i.e. \code{conf["bam"]}) for an already-built ORFik experiment
#'
#' The inverse direction of \code{\link{qc_diagnostics_path}}: some
#' callers (e.g. \code{bad_pshifting_report()}) only have the ORFik
#' experiment object \code{df}, not the original \code{conf["bam"]}
#' string. \code{resFolder(df) == conf["bam"]/aligned} for
#' massiveNGSpipe-created experiments (verified live against a real
#' processed experiment), so its parent directory is \code{conf["bam"]}.
#' @param df an ORFik experiment
#' @return character path
#' @noRd
bam_dir_from_df <- function(df) {
  dirname(resFolder(df))
}

#' True (score-weighted) read-length distribution for one ofst file
#'
#' Every existing interactive tool that computes a length distribution
#' (shared_scripts/.../shift_study_app.R, mNGSp_manual_reshift.R) does a
#' plain \code{table(readWidths(x))} over an ofst's rows -- but ofst rows
#' are collapsed/deduplicated sequences, each carrying its true read
#' count in \code{mcols()$score} (same collapse-weighting issue already
#' fixed for alignment metrics, see \code{get_expanded_alignment_metrics()},
#' R/bam_utils.R). Verified live on real production data: the row-count
#' view and the score-weighted view can peak at completely different
#' read lengths. This weights by \code{score} to get the TRUE
#' distribution; row-count-only ofst files (no \code{score} column, e.g.
#' pre-dating multiplicity tracking) fall back to weight 1 per row.
#' @param ofst_path character, path to a \code{.ofst} file
#' @return data.table(read_length, count, percent), sorted by
#' read_length. \code{count} is the total (score-weighted) read count.
#' A 0-row table (not an error) for an empty/0-read ofst file.
#' @noRd
read_length_distribution <- function(ofst_path) {
  x <- ORFik::fimport(ofst_path)
  if (length(x) == 0) {
    return(data.table::data.table(read_length = integer(), count = numeric(), percent = numeric()))
  }
  score <- S4Vectors::mcols(x)$score
  if (is.null(score)) score <- rep(1, length(x))
  dt <- data.table::data.table(read_length = ORFik::readWidths(x), score = score)
  agg <- dt[, .(count = sum(score)), by = read_length][order(read_length)]
  agg[, percent := round(100 * count / sum(count), 2)]
  agg[]
}

#' Save one sample's read-length distribution CSV
#' @param ofst_path character, the .ofst (or _pshifted.ofst) file to summarize
#' @param out_dir character, destination folder (e.g.
#' \code{<ofst_dir>/read_length_distribution})
#' @param run_id character, sample identifier for the per-sample filename
#' @return invisible(data.table), this sample's own distribution (so a
#' caller building its own aggregate doesn't need to re-read the CSV back)
#' @noRd
save_read_length_distribution <- function(ofst_path, out_dir, run_id) {
  dist <- read_length_distribution(ofst_path)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(dist, file.path(out_dir, paste0(run_id, ".csv")))
  invisible(dist)
}

#' Combine every per-sample CSV in a read_length_distribution folder into
#' one study-level \code{00_aggregated.csv}
#'
#' \code{00_}-prefixed, matching the existing aggregate-file convention
#' already used for \code{00_STAR_LOG_table.csv}/\code{00_STAR_LOG_plot.pdf}
#' (\code{aligned/LOGS*/}). Re-reads whatever per-sample CSVs currently
#' exist in \code{out_dir} rather than taking them as an argument, so a
#' resumed run's aggregate is always built from every sample done so far,
#' not just the ones newly converted in this call.
#' @param out_dir character, the read_length_distribution folder
#' @return invisible(NULL)
#' @noRd
aggregate_read_length_distribution <- function(out_dir) {
  files <- list.files(out_dir, pattern = "\\.csv$", full.names = TRUE)
  files <- files[basename(files) != "00_aggregated.csv"]
  if (length(files) == 0) return(invisible(NULL))
  all_dt <- data.table::rbindlist(lapply(files, data.table::fread), fill = TRUE)
  if (nrow(all_dt) == 0) return(invisible(NULL))
  agg <- all_dt[, .(count = sum(count)), by = read_length][order(read_length)]
  agg[, percent := round(100 * count / sum(count), 2)]
  data.table::fwrite(agg, file.path(out_dir, "00_aggregated.csv"))
  invisible(NULL)
}

#' Canonical QC-diagnostics file path for one experiment
#'
#' Computable identically before or after the ORFik experiment object
#' exists -- verified directly against a real processed experiment:
#' \code{QCfolder(df) == file.path(conf["bam"], "aligned", "QC_STATS/")}.
#' This lets every pipeline stage from \code{pipeline_trim()} onward write
#' its own diagnostics using only \code{conf["bam"]}, without needing to
#' construct (or wait for) the ORFik experiment object itself.
#' @param bam_dir character, \code{conf["bam"]} for this experiment
#' @return character path (not guaranteed to exist yet)
#' @noRd
qc_diagnostics_path <- function(bam_dir) {
  file.path(bam_dir, "aligned", "QC_STATS", "qc_diagnostics.rds")
}

#' Read-modify-write update of one experiment's QC diagnostics record
#'
#' Never clobbers fields written by an earlier pipeline stage -- each
#' named argument in \code{...} overwrites only that one field.
#' @param path character, see \code{\link{qc_diagnostics_path}}
#' @param ... named fields to set/overwrite
#' @return invisible(list), the full updated record
#' @noRd
update_qc_diagnostics <- function(path, ...) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  existing <- if (file.exists(path)) readRDS(path) else list()
  new_fields <- list(...)
  existing[names(new_fields)] <- new_fields
  existing$last_updated <- Sys.time()
  saveRDS(existing, path)
  invisible(existing)
}

#' Read one experiment's QC diagnostics record
#' @param path character, see \code{\link{qc_diagnostics_path}}
#' @return list, \code{list()} if nothing has been recorded yet
#' @noRd
read_qc_diagnostics <- function(path) {
  if (!file.exists(path)) return(list())
  readRDS(path)
}

#' Flag an experiment as too-few-raw-reads, from already-computed trim stats
#'
#' Reuses \code{ORFik::trimming.table()}'s own aggregation (the same one
#' \code{save_report()} already relies on) rather than re-parsing
#' individual fastp JSON files.
#' @param bam_dir character, \code{conf["bam"]} for this experiment
#' @param trimmed_dir character, this experiment's trim output directory
#' (holds the per-run fastp JSON reports)
#' @param min_raw_reads numeric, threshold below which a run counts as
#' too-few-reads (see \code{config$min_raw_reads_pshift})
#' @return invisible(list), the updated diagnostics record
#' @noRd
check_too_few_reads <- function(bam_dir, trimmed_dir, min_raw_reads) {
  path <- qc_diagnostics_path(bam_dir)
  stats <- tryCatch(ORFik::trimming.table(trimmed_dir), error = function(e) NULL)
  if (is.null(stats) || nrow(stats) == 0 || is.null(stats$raw_reads)) {
    return(update_qc_diagnostics(path, too_few_reads = NA, raw_reads = NA_real_,
                                 trim_reads = NA_real_))
  }
  raw_reads <- sum(stats$raw_reads, na.rm = TRUE)
  trim_reads <- if (!is.null(stats$trim_reads)) sum(stats$trim_reads, na.rm = TRUE) else NA_real_
  update_qc_diagnostics(path, too_few_reads = raw_reads < min_raw_reads,
                        raw_reads = raw_reads, trim_reads = trim_reads)
}

#' Flag an experiment as likely-wrong-organism, from already-written STAR stats
#'
#' Mirrors \code{status_per_study()}'s own
#' \code{full_process_SINGLE.csv}-then-\code{full_process.csv} fallback
#' lookup (factored out as \code{\link{full_process_csv_path}} so the two
#' call sites can't drift apart).
#' @param bam_dir character, \code{conf["bam"]} for this experiment
#' @param min_alignment_rate numeric, percent, threshold below which
#' counts as likely-wrong-organism (see \code{config$min_alignment_rate_pshift})
#' @return invisible(list), the updated diagnostics record
#' @noRd
check_alignment_rate <- function(bam_dir, min_alignment_rate) {
  path <- qc_diagnostics_path(bam_dir)
  csv_path <- full_process_csv_path(bam_dir)
  if (is.null(csv_path)) {
    # File genuinely doesn't exist yet -- distinct from "checked and passed".
    return(update_qc_diagnostics(path, wrong_organism = NA, mapping_rate_pct = NA_real_))
  }
  stats <- tryCatch(data.table::fread(csv_path), error = function(e) NULL)
  rate_col <- intersect(c("Uniquely mapped reads %", "total mapped reads %"), colnames(stats))
  if (is.null(stats) || nrow(stats) == 0 || length(rate_col) == 0) {
    return(update_qc_diagnostics(path, wrong_organism = NA, mapping_rate_pct = NA_real_))
  }
  rate <- mean(suppressWarnings(as.numeric(gsub("%", "", stats[[rate_col[1]]]))), na.rm = TRUE)
  update_qc_diagnostics(path, wrong_organism = !is.na(rate) && rate < min_alignment_rate,
                        mapping_rate_pct = rate)
}

#' Path to an experiment's full_process*.csv, or NULL if neither exists yet
#'
#' Same fallback order \code{status_per_study()} already uses
#' (\code{R/pipeline_progress.R}): \code{full_process_SINGLE.csv} first,
#' else \code{full_process.csv}.
#' @param bam_dir character, \code{conf["bam"]} for this experiment
#' @return character path, or NULL
#' @noRd
full_process_csv_path <- function(bam_dir) {
  single <- file.path(bam_dir, "full_process_SINGLE.csv")
  plain <- file.path(bam_dir, "full_process.csv")
  if (file.exists(single)) return(single)
  if (file.exists(plain)) return(plain)
  NULL
}

#' Flag an experiment's adapter/barcode trimming as likely-suspect
#'
#' From \code{adapter_barcode_table.csv} (written by \code{pipeline_trim()}
#' for RFP libraries only) -- a high fraction of reads with no adapter
#' removed at all is the closest existing proxy signal for a
#' wrong/mismatched adapter or barcode configuration.
#' @param bam_dir character, \code{conf["bam"]} for this experiment
#' @param trimmed_dir character, this experiment's trim output directory
#' (holds \code{adapter_barcode_table.csv})
#' @param max_no_adapter_removed_pct numeric, percent, threshold above
#' which counts as likely-bad-adapter-barcode
#' @return invisible(list), the updated diagnostics record
#' @noRd
check_adapter_barcode_quality <- function(bam_dir, trimmed_dir, max_no_adapter_removed_pct) {
  path <- qc_diagnostics_path(bam_dir)
  csv_path <- file.path(trimmed_dir, "adapter_barcode_table.csv")
  if (!file.exists(csv_path)) {
    return(update_qc_diagnostics(path, bad_adapter_barcode = NA, no_adapter_removed_pct = NA_real_))
  }
  stats <- tryCatch(data.table::fread(csv_path), error = function(e) NULL)
  col <- "reads_no_adapter_removed fastp(%)"
  if (is.null(stats) || nrow(stats) == 0 || !(col %in% colnames(stats))) {
    return(update_qc_diagnostics(path, bad_adapter_barcode = NA, no_adapter_removed_pct = NA_real_))
  }
  pct <- mean(suppressWarnings(as.numeric(stats[[col]])), na.rm = TRUE)
  update_qc_diagnostics(path, bad_adapter_barcode = !is.na(pct) && pct > max_no_adapter_removed_pct,
                        no_adapter_removed_pct = pct)
}

#' Pick one likely primary cause from an experiment's recorded diagnostics
#'
#' Fixed precedence, not a score: \code{too_few_reads > wrong_organism >
#' adapter_barcode_umi > unknown}. Too-few-reads is checked first because
#' a tiny denominator makes the alignment-rate signal meaningless on its
#' own; wrong-organism second because very low mapping despite adequate
#' reads is a fairly unambiguous signal; adapter/barcode last because
#' it's the least direct signal and only meaningful once the two harder
#' failures are ruled out. All underlying signals stay in the record
#' regardless of which one wins -- this only picks a label for triage.
#' @param bam_dir character, \code{conf["bam"]} for this experiment
#' @param config the mNGSp config object (unused today, accepted for a
#' future where the precedence itself becomes configurable)
#' @return character, the chosen \code{primary_cause}
#' @noRd
classify_bad_shift_cause <- function(bam_dir, config = NULL) {
  path <- qc_diagnostics_path(bam_dir)
  diag <- read_qc_diagnostics(path)
  cause <- if (isTRUE(diag$too_few_reads)) {
    "too_few_reads"
  } else if (isTRUE(diag$wrong_organism)) {
    "wrong_organism"
  } else if (isTRUE(diag$bad_adapter_barcode)) {
    "adapter_barcode_umi"
  } else {
    "unknown"
  }
  update_qc_diagnostics(path, primary_cause = cause)
  cause
}

#' Per-length frame preference with gene and site concentration (diagnostic)
#'
#' A strong frame preference at some read length can come from one locus (PRJNA430238, 2026-10-21: the apparent
#' 25-27 nt phase was a single Phf21a 5' end, inflated by counting each read once per isoform). For each length this
#' reports the CDS-body reads, how many genes contribute, the share of the top gene and of the top single position,
#' and the best-frame fraction with and without the top gene. Pass one CDS per gene (e.g. longest) so reads are not
#' counted once per isoform. Diagnostic only: never use it to drop data automatically.
#' @param reads GRanges (ofst import; read multiplicity in mcols score, else 1 per read)
#' @param cds GRangesList, one CDS per gene
#' @param lengths integer, read lengths to report
#' @param end "five" (5' ends) or "three" (3' ends)
#' @param flank integer, positions excluded at the CDS start and end
#' @return data.table(len, reads, genes, top_gene, top_gene_frac, top_site_frac, best_frame, frame_frac,
#'  frame_frac_without_top, concentrated)
#' @export
frame_concentration_table <- function(reads, cds, lengths = 20:40, end = c("five", "three"), flank = 30L) {
  end <- match.arg(end)
  empty <- data.table(len = integer(), reads = numeric(), genes = integer(), top_gene = character(), top_gene_frac = numeric(),
                      top_site_frac = numeric(), best_frame = integer(), frame_frac = numeric(),
                      frame_frac_without_top = numeric(), concentrated = logical())
  if (!length(reads) || !length(cds)) return(empty)
  w <- S4Vectors::mcols(reads)$score; if (is.null(w)) w <- rep(1, length(reads))
  len <- ORFik::readWidths(reads); k <- len %in% lengths
  reads <- reads[k]; w <- w[k]; len <- len[k]
  if (!length(reads)) return(empty)
  p1 <- GenomicRanges::resize(GenomicRanges::granges(reads), 1L, fix = if (end == "five") "start" else "end")
  m <- GenomicFeatures::mapToTranscripts(p1, cds)
  if (!length(m)) return(empty)
  xi <- S4Vectors::mcols(m)$xHits; tx <- as.character(GenomicRanges::seqnames(m))
  cw <- ORFik::widthPerGroup(cds, TRUE)[tx]
  d <- data.table(len = len[xi], w = w[xi], tx = tx, p = GenomicRanges::start(m), cw = cw)[p > flank & p < cw - flank]
  if (!nrow(d)) return(empty)
  d[, frame := (p - 1L) %% 3L]
  bestfrac <- function(fr, ww) { s <- tapply(ww, fr, sum); if (!length(s)) NA_real_ else max(s) / sum(s) }
  out <- d[, {
    g <- .SD[, .(w = sum(w)), by = tx][order(-w)]
    site <- .SD[, .(w = sum(w)), by = .(tx, p)]
    s <- tapply(w, frame, sum)
    rest <- tx != g$tx[1]
    .(reads = sum(w), genes = nrow(g), top_gene = g$tx[1], top_gene_frac = g$w[1] / sum(w),
      top_site_frac = max(site$w) / sum(w), best_frame = as.integer(names(s)[which.max(s)]),
      frame_frac = max(s) / sum(s), frame_frac_without_top = if (any(rest)) bestfrac(frame[rest], w[rest]) else NA_real_)
  }, by = len][order(len)]
  out[, concentrated := top_gene_frac > 0.3 | top_site_frac > 0.2]
  out[]
}
