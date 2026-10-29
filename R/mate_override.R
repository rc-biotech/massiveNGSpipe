#' Opt-in mate selection for paired-end Ribo-seq whose footprint is in read 2
#'
#' Some paired-end libraries carry the footprint only in read 2 (read 1 is adapter / index + poly-G), e.g. PRJNA1301971
#' (2026-10-29: R2 = GAT + 5 nt UMI | reverse complement of the footprint | 2 nt | AGATCGGAAGAGCGT). The default paired
#' path uses one read file per run and cannot recover these. For runs listed in \code{<trim dir>/mate_override_manual.csv}
#' (columns \code{Run}, \code{mate} 1/2, \code{reverse} TRUE/FALSE, \code{adapter}, \code{trim_front}, \code{trim_tail},
#' \code{source}, optional \code{note}) the trim stage trims and collapses only the chosen mate with trimcollapse
#' (\code{trim_front} nt removed from the 5' end, \code{trim_tail} nt removed after the adapter cut), names the output
#' by Run, reverse-complements it when \code{reverse} is TRUE, and records what was done in
#' \code{<trim dir>/.mate_override/<Run>.json}. Alignment and experiment creation then treat these runs as single-end;
#' the metadata (LibraryLayout) is not changed. Default is off: without the table nothing changes.
#' @name mate_override
NULL

#' @noRd
mate_override_path <- function(trimmed_dir) file.path(trimmed_dir, "mate_override_manual.csv")

#' Mate override table of a study (strict parsing)
#' @param trimmed_dir the study's trim directory
#' @return data.table (empty when no table exists)
#' @noRd
mate_overrides <- function(trimmed_dir) {
  f <- mate_override_path(trimmed_dir)
  if (!file.exists(f)) return(data.table(Run = character()))
  x <- fread(f, colClasses = "character")
  need <- c("Run", "mate", "reverse", "adapter", "trim_front", "trim_tail", "source")
  if (!all(need %in% names(x))) stop(f, " must have columns ", paste(need, collapse = ", "))
  if (anyDuplicated(x$Run)) stop(f, ": duplicated Run entries: ", paste(unique(x$Run[duplicated(x$Run)]), collapse = ", "))
  if (!all(x$mate %in% c("1", "2"))) stop(f, ": mate must be 1 or 2")
  rev <- toupper(trimws(x$reverse))
  if (!all(rev %in% c("TRUE", "FALSE"))) stop(f, ": reverse must be TRUE or FALSE")
  if (!all(grepl("^([ACGT]+|disable)$", x$adapter))) stop(f, ": adapter must be a DNA sequence or 'disable'")
  if (!all(grepl("^[0-9]+$", x$trim_front) & grepl("^[0-9]+$", x$trim_tail))) stop(f, ": trim_front/trim_tail must be integers >= 0")
  if (any(!nzchar(x$source))) stop(f, ": source (evidence) is required for every row")
  x[, `:=`(mate = as.integer(mate), reverse = rev == "TRUE", trim_front = as.integer(trim_front), trim_tail = as.integer(trim_tail))]
  x
}

#' Runs of a study that use a mate override
#' @noRd
mate_override_runs <- function(trimmed_dir) mate_overrides(trimmed_dir)$Run

#' Treat overridden runs as single-end in an in-memory run table (metadata files are not touched)
#' @param runs data.table with Run, LibraryLayout
#' @param trimmed_dir the study's trim directory
#' @return runs, with LibraryLayout "SINGLE" for overridden runs
#' @noRd
apply_mate_override_layout <- function(runs, trimmed_dir) {
  ov <- mate_override_runs(trimmed_dir)
  if (!length(ov)) return(runs)
  runs <- copy(runs)
  runs[Run %in% ov, LibraryLayout := "SINGLE"]
  runs
}

#' Trim + collapse one run from its overridden mate
#' @param run Run id
#' @param filenames the run's raw read files (length 2 for paired)
#' @param target_dir the study folder (trimcollapse writes to <target_dir>/trim)
#' @param trimmed_dir the study's trim directory
#' @return one-row data.table for the barcode table / trim marker
#' @noRd
trim_with_mate_override <- function(run, filenames, target_dir, trimmed_dir) {
  o <- mate_overrides(trimmed_dir)[Run == run]
  stopifnot(nrow(o) == 1)
  file <- filenames[o$mate]
  if (is.na(file) || !file.exists(file)) stop("mate_override: read ", o$mate, " file missing for ", run)
  out <- trimcollapse_run(file, target_dir, adapter = o$adapter, trim.front = o$trim_front, tail_after_adapter = o$trim_tail)
  json <- file.path(trimmed_dir, paste0("report_", trimcollapse_id(file), ".json"))
  final <- file.path(trimmed_dir, "SINGLE", paste0("collapsed_trimmed_", run, ".fasta.gz"))
  final_json <- file.path(trimmed_dir, paste0("report_", run, ".json"))
  if (out != final) { if (!file.rename(out, final)) stop("mate_override: cannot rename ", out) }
  if (json != final_json) { if (!file.rename(json, final_json)) stop("mate_override: cannot rename ", json) }
  n <- NA_integer_
  if (o$reverse) {
    seqs <- Biostrings::readDNAStringSet(final)
    rc <- Biostrings::reverseComplement(seqs)
    tmp <- file.path(dirname(final), paste0(".", basename(final), ".mate_tmp", Sys.getpid(), ".gz"))
    Biostrings::writeXStringSet(rc, tmp, compress = TRUE)
    if (!file.rename(tmp, final)) { unlink(tmp); stop("mate_override: cannot replace ", final) }
    n <- length(rc)
  }
  write_json_atomic(list(Run = run, input = basename(file), mate = o$mate, reverse = o$reverse, adapter = o$adapter,
                         trim_front = o$trim_front, trim_tail = o$trim_tail, source = o$source,
                         output = basename(final), output_md5 = file_md5(final), n_seqs = n,
                         time = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                    file.path(trimmed_dir, ".mate_override", paste0(run, ".json")))
  fast_trim_state$collapsed <- final
  data.table(id = run, adapter = o$adapter, barcode_detected = o$trim_front + o$trim_tail > 0,
             barcode5p_size = o$trim_front, barcode3p_size = o$trim_tail, mate_override = o$mate, reverse = o$reverse)
}
