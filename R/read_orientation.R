#' Opt-in read orientation for libraries sequenced from the RNA 3' side
#'
#' Some Ribo-seq libraries are read antisense to the RNA (read 1 starts at the 3' adapter side, e.g. reverse-oriented
#' small-RNA or oligo-dT designs; 2026-10-14 barcode audit: PRJNA453992, PRJNA835330). Their alignments land on the
#' opposite strand, which breaks P-shifting and counting. For runs explicitly listed in
#' \code{<trim dir>/reads_reverse_manual.csv} (columns \code{Run}, \code{reverse} TRUE/FALSE, \code{source}, \code{note}),
#' the collapsed FASTA is reverse-complemented right after collapsing and before the sample's collapsed flag is set.
#' Headers (and so the collapsed read-multiplicity weights) and lengths are unchanged. A per-run transaction sidecar
#' (\code{.orientation/<Run>.json}) makes it idempotent and crash-safe (see \code{orient_collapsed_fasta}). Default is off; paired-end runs are refused (orientation of read 1 must be confirmed separately).
#' @name read_orientation
NULL

#' @noRd
reads_reverse_manual_path <- function(trimmed_dir) file.path(trimmed_dir, "reads_reverse_manual.csv")

#' Per-run orientation transaction sidecar (no shared file between parallel collapse workers)
#' @noRd
orientation_state_path <- function(trimmed_dir, run) file.path(trimmed_dir, ".orientation", paste0(run, ".json"))

#' Runs of a study marked for reverse-complementing (strict parsing)
#' @param trimmed_dir the study's trim directory
#' @return character vector of Run ids with reverse == TRUE
#' @noRd
reads_reverse_runs <- function(trimmed_dir) {
  f <- reads_reverse_manual_path(trimmed_dir)
  if (!file.exists(f)) return(character())
  x <- fread(f, colClasses = "character")
  if (!all(c("Run", "reverse") %in% names(x))) stop(f, " must have columns Run, reverse")
  v <- toupper(trimws(x$reverse))
  if (any(!v %in% c("TRUE", "FALSE"))) stop(f, ": reverse must be TRUE or FALSE (got ", paste(unique(x$reverse[!v %in% c("TRUE", "FALSE")]), collapse = ", "), ")")
  if (anyDuplicated(x$Run)) stop(f, ": duplicated Run entries: ", paste(unique(x$Run[duplicated(x$Run)]), collapse = ", "))
  x$Run[v == "TRUE"]
}

#' @noRd
file_md5 <- function(path) unname(tools::md5sum(path))

#' @noRd
write_json_atomic <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  jsonlite::write_json(obj, tmp, auto_unbox = TRUE)
  if (!file.rename(tmp, path)) { unlink(tmp); stop("cannot write ", path) }
}

#' Reverse-complement one collapsed FASTA if its run is marked (idempotent, crash-safe)
#'
#' Transaction per run in \code{<trim>/.orientation/<Run>.json}: the input md5 and the planned output md5 are recorded
#' (state "prepared") before the atomic replacement, and state "done" after it. On a retry the current file's md5
#' decides: equal to the recorded output -> already flipped (state repaired to done); equal to the recorded input or
#' a new file (re-collapsed) -> flip. Each run has its own sidecar, so parallel collapse workers never share state.
#' @param collapsed_file path to collapsed_<run>.fasta(.gz)
#' @param run Run id
#' @param trimmed_dir the study's trim directory (holds the manual table and sidecars)
#' @param paired logical, TRUE refuses (read-1 orientation must be confirmed separately)
#' @return "flipped", "already_flipped" or "not_marked", invisibly
#' @noRd
orient_collapsed_fasta <- function(collapsed_file, run, trimmed_dir, paired = FALSE) {
  if (!run %in% reads_reverse_runs(trimmed_dir)) return(invisible("not_marked"))
  if (paired) stop("reads_reverse_manual: paired-end run ", run, " refused; confirm read-1 orientation separately")
  stopifnot(file.exists(collapsed_file))
  sp <- orientation_state_path(trimmed_dir, run)
  st <- if (file.exists(sp)) jsonlite::read_json(sp) else NULL
  cur <- file_md5(collapsed_file)
  if (!is.null(st) && identical(st$output_md5, cur)) {
    if (!identical(st$state, "done")) { st$state <- "done"; write_json_atomic(st, sp) }   # crash after replace: repair
    return(invisible("already_flipped"))
  }
  seqs <- Biostrings::readDNAStringSet(collapsed_file)
  rc <- Biostrings::reverseComplement(seqs)
  stopifnot(identical(names(rc), names(seqs)), identical(width(rc), width(seqs)))
  gz <- grepl("\\.gz$", collapsed_file)
  tmp <- file.path(dirname(collapsed_file), paste0(".", basename(collapsed_file), ".orient_tmp", Sys.getpid(), if (gz) ".gz" else ""))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  Biostrings::writeXStringSet(rc, tmp, compress = gz)
  st <- list(Run = run, file = basename(collapsed_file), input_md5 = cur, output_md5 = file_md5(tmp), n_seqs = length(rc),
             state = "prepared", time = format(Sys.time(), tz = "UTC", usetz = TRUE))
  write_json_atomic(st, sp)                                         # planned output recorded BEFORE the replacement
  if (!file.rename(tmp, collapsed_file)) stop("orientation: cannot replace ", collapsed_file, " (original left unchanged)")
  st$state <- "done"; write_json_atomic(st, sp)
  invisible("flipped")
}

#' Study summary of orientation transactions (read outside parallel workers)
#' @param trimmed_dir the study's trim directory
#' @return data.table, one row per run sidecar
#' @noRd
orientation_summary <- function(trimmed_dir) {
  f <- list.files(file.path(trimmed_dir, ".orientation"), "\\.json$", full.names = TRUE)
  rbindlist(lapply(f, function(x) as.data.table(jsonlite::read_json(x))), fill = TRUE)
}

#' Mark runs for reverse-complementing (or undo) and reset their processing from collapse
#'
#' Writes/updates \code{reads_reverse_manual.csv} (upsert by Run) and an audit row in
#' \code{<project>/read_orientation_log.csv}, then clears the runs' per-sample markers from \code{collapsed} onward
#' and the experiment flags from \code{collapsed} onward, so the next standard run re-collapses (needs the trimmed
#' reads; for online data with deleted intermediates also clear \code{trim} / \code{fetch}, e.g. via
#' \code{apply_trim_fix_to_sample()}), re-orients, re-aligns and re-shifts with fresh offsets.
#' @param exp experiment name
#' @param runs Run ids
#' @param reverse logical (TRUE = reverse-complement, FALSE = native orientation)
#' @param config pipeline config
#' @param source,note evidence (e.g. "raw reads: insert + RC small-RNA 5' adapter") and free text
#' @return invisible(NULL)
#' @export
set_read_orientation <- function(exp, runs, reverse, config, source, note = "") {
  df <- read.experiment(exp, validate = FALSE, output.env = new.env())
  stopifnot(all(runs %in% ORFik::runIDs(df)), is.logical(reverse), length(reverse) == 1, nzchar(source))
  trimmed_dir <- file.path(dirname(libFolder(df)), "trim")
  f <- reads_reverse_manual_path(trimmed_dir)
  cur <- if (file.exists(f)) fread(f, colClasses = "character") else data.table(Run = character(), reverse = character(),
                                                                                 source = character(), note = character())
  cur <- cur[!Run %in% runs]
  new <- data.table(Run = runs, reverse = as.character(reverse), source = source, note = note)
  fwrite(rbind(cur, new, fill = TRUE), f)
  log_f <- file.path(config$project, "read_orientation_log.csv")
  fwrite(cbind(data.table(time = format(Sys.time(), tz = "UTC", usetz = TRUE), experiment = exp), new),
         log_f, append = file.exists(log_f))
  steps <- names(config$flag); steps <- steps[which(steps == "collapsed"):length(steps)]
  for (s in intersect(steps, c("collapsed", "aligned", "ofst", "ofst_unique", "covrle", "bigwig")))
    for (r in runs) unlink(file.path(sample_flag_dir(config, s, exp), paste0(r, ".rds")))
  remove_flag_all_exp(config, steps = steps, exps = exp)
  # stale state for the affected runs only (siblings untouched): orientation sidecars, and their P-shift offsets in the
  # shifting tables (all and unique mappers), so reuse_shifts_if_existing cannot reuse opposite-strand offsets
  unlink(orientation_state_path(trimmed_dir, runs))
  for (sh in list.files(libFolder(df), "^shifting_table\\.rds$", recursive = TRUE, full.names = TRUE)) {
    tab <- readRDS(sh)
    drop <- vapply(names(tab), function(n) any(vapply(runs, function(r) grepl(r, basename(n), fixed = TRUE), logical(1))), logical(1))
    if (any(drop)) { file.copy(sh, paste0(sh, ".bak_orientation_", format(Sys.time(), "%Y%m%d_%H%M%S"))); saveRDS(tab[!drop], sh) }
  }
  message(exp, ": ", length(runs), " runs set to reverse = ", reverse, "; flags reset from collapsed; their shift offsets dropped")
  invisible(NULL)
}

#' Strand QC: fraction of CDS-overlapping reads on the sense strand, per sample
#'
#' Ribo-seq footprints should be sense to the CDS. A sample with mostly antisense CDS reads was probably sequenced from
#' the RNA 3' side (see \code{set_read_orientation()}). This only FLAGS; it never changes orientation. Uses the pre-shift
#' ofst (score-weighted), up to \code{max_rows} alignments per sample; samples with few CDS reads stay "unknown".
#' @param exp experiment name
#' @param max_rows maximum alignment rows read per sample
#' @param min_cds_reads minimum weighted CDS-overlapping reads for a verdict
#' @param n_samples check at most this many samples (first ones), e.g. 2 for a server-wide sweep
#' @return data.table(Run, cds_reads, sense_frac, verdict: sense / antisense / mixed / unknown)
#' @export
strand_qc <- function(exp, max_rows = 2e6, min_cds_reads = 1000, n_samples = Inf) {
  df <- read.experiment(exp, validate = FALSE, output.env = new.env())
  if (nrow(df) > n_samples) df <- df[seq_len(n_samples), ]
  cds <- unlist(loadRegion(df, "cds"), use.names = FALSE)
  rbindlist(lapply(seq_len(nrow(df)), function(i) {
    f <- tryCatch(filepath(df[i, ], "ofst"), error = function(e) NA_character_)
    if (is.na(f) || !file.exists(f)) return(data.table(Run = ORFik::runIDs(df)[i], cds_reads = NA_real_, sense_frac = NA_real_, verdict = "unknown"))
    x <- ORFik::fimport(f)
    if (length(x) > max_rows) x <- x[sort(sample.int(length(x), max_rows))]
    gr <- GenomicRanges::granges(x); w <- if (!is.null(S4Vectors::mcols(x)$score)) S4Vectors::mcols(x)$score else rep(1, length(x))
    sense <- IRanges::overlapsAny(gr, cds, ignore.strand = FALSE)
    any_ <- IRanges::overlapsAny(gr, cds, ignore.strand = TRUE)
    anti <- any_ & !sense
    tot <- sum(as.numeric(w[sense | anti])); sf <- if (tot > 0) sum(as.numeric(w[sense])) / tot else NA_real_
    v <- if (tot < min_cds_reads) "unknown" else if (sf >= 0.8) "sense" else if (sf <= 0.2) "antisense" else "mixed"
    data.table(Run = ORFik::runIDs(df)[i], cds_reads = tot, sense_frac = round(sf, 3), verdict = v)
  }))
}
