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
mate_override_runs <- function(trimmed_dir) {
  auto <- mate_override_auto_path(trimmed_dir)
  unique(c(mate_overrides(trimmed_dir)$Run, if (file.exists(auto)) fread(auto, colClasses = "character")$Run))
}

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

#' @noRd
mate_override_auto_path <- function(trimmed_dir) file.path(trimmed_dir, "mate_override_auto.csv")

#' Which mate of a paired-end Ribo-seq run holds the footprint (first reads only, fast)
#'
#' Read 1 is unusable when most reads start with one fixed sequence (adapter / index read-through, e.g. PRJNA1301971:
#' 96% ACGTCTGAACTCCAGTCAC + index, then poly-G). Read 2 is then used when its reads carry an insert followed by an
#' Illumina adapter seed. Anything else keeps read 1 (the default paired path).
#' @param file1,file2 read 1 / read 2 fastq(.gz)
#' @param n reads sampled from the start of each file
#' @return list(mate, reverse, r1_top_prefix, r2_top_prefix, r2_insert_adapter, reason)
#' @noRd
paired_footprint_mate <- function(file1, file2, n = 1e5, prefix_nt = 20, junk_min = 0.5, insert_min = 0.3) {
  seqs <- function(f) { con <- if (grepl("\\.gz$", f)) gzfile(f, "r") else file(f, "r"); on.exit(close(con))
    l <- readLines(con, n = 4 * n); l[seq(2, length(l), by = 4)] }
  top_prefix <- function(s) { p <- substr(s[nchar(s) >= prefix_nt], 1, prefix_nt); if (!length(p)) return(0); max(table(p)) / length(p) }
  insert_adapter <- function(s) { pos <- regexpr("AGATCGGAAGAGC", s, fixed = TRUE); mean(pos >= 18) }
  s1 <- seqs(file1); s2 <- seqs(file2)
  r1 <- top_prefix(s1); r2 <- top_prefix(s2); i2 <- insert_adapter(s2)
  use2 <- r1 >= junk_min && r2 < 0.2 && i2 >= insert_min
  list(mate = if (use2) 2L else 1L, reverse = use2, r1_top_prefix = round(r1, 3), r2_top_prefix = round(r2, 3),
       r2_insert_adapter = round(i2, 3),
       reason = sprintf("auto: read 1 top %d-nt prefix %.0f%% of reads, read 2 insert+adapter %.0f%%", prefix_nt, 100 * r1, 100 * i2))
}

#' Trim + collapse a paired run from read 2 (automatic footprint-mate detection)
#'
#' Read 2 goes through the normal single-end barcode/UMI detection and trimcollapse path (via a Run-named symlink),
#' the collapsed reads are reverse-complemented back to the RNA sense and the decision is recorded in
#' \code{mate_override_auto.csv} (used like the manual table by alignment and experiment creation).
#' @noRd
trim_from_read2_auto <- function(study_sample, file2, target_dir, trimmed_dir, mode, decision) {
  run <- study_sample$Run
  tmp <- file.path(trimmed_dir, ".mate2_input", run); dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  link <- file.path(tmp, paste0(run, if (grepl("\\.gz$", file2)) ".fastq.gz" else ".fastq"))
  if (!file.symlink(normalizePath(file2), link)) stop("mate override: cannot link read 2 for ", run)
  ss <- copy(study_sample)[, LibraryLayout := "SINGLE"]
  fast_trim_state$collapsed <- NULL
  barcode_dt <- run_barcode_detection_and_trim(ss, tmp, target_dir, trimmed_dir, mode, link, NULL)
  final <- fast_trim_state$collapsed
  if (is.null(final) || !file.exists(final))
    stop("mate override (auto) for ", run, ": read 2 did not go through the fast trim path (adapter '", barcode_dt$adapter,
         "'); add a manual row to ", mate_override_path(trimmed_dir))
  seqs <- Biostrings::readDNAStringSet(final)
  tmpf <- file.path(dirname(final), paste0(".", basename(final), ".mate_tmp", Sys.getpid(), ".gz"))
  Biostrings::writeXStringSet(Biostrings::reverseComplement(seqs), tmpf, compress = TRUE)
  if (!file.rename(tmpf, final)) { unlink(tmpf); stop("mate override: cannot replace ", final) }
  row <- data.table(Run = run, mate = 2L, reverse = TRUE, adapter = barcode_dt$adapter,
                    trim_front = barcode_dt$barcode5p_size, trim_tail = barcode_dt$barcode3p_size, source = decision$reason)
  f <- mate_override_auto_path(trimmed_dir)
  old <- if (file.exists(f)) fread(f, colClasses = "character")[Run != run] else NULL
  fwrite(rbind(old, row[, lapply(.SD, as.character)], fill = TRUE), f)
  write_json_atomic(c(decision, list(Run = run, output = basename(final), output_md5 = file_md5(final), n_seqs = length(seqs),
                                     time = format(Sys.time(), tz = "UTC", usetz = TRUE))),
                    file.path(trimmed_dir, ".mate_override", paste0(run, ".json")))
  fast_trim_state$collapsed <- final
  barcode_dt[, `:=`(mate_override = 2L, reverse = TRUE)]
  barcode_dt
}
