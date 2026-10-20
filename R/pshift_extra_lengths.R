#' Extend P-site shift tables with extra periodic read lengths (additive)
#'
#' The fixed accepted read lengths (\code{config$accepted_lengths_rpf}, default 20, 21, 25:33) leave periodic reads
#' unshifted (median 12.6-32.6% of aligned reads in human samples; see the cleanup project's C4.1 evidence). This adds
#' offsets for extra lengths that pass ORFik's periodicity test, per sample, and never changes existing rows: every
#' (length, offset) already in a sample's table is kept as is (including manually reviewed shifts), so existing
#' P-sites are unchanged and a sample can only gain reads.
#' @param df ORFik experiment (Ribo-seq)
#' @param shifting_table named list of per-sample shift tables (fraction = read length, offsets_start), as
#'  \code{shifts_load()} returns; names are the samples' ofst paths
#' @param extra_lengths read lengths to test (default 22:24 and 34:40)
#' @param strict.fft passed to \code{ORFik::detectRibosomeShifts}
#' @param BPPARAM BiocParallel parameters (one sample per worker)
#' @return list of extended tables, same names; attribute "added" = data.table(sample, length, offset)
#' @export
extend_shift_table <- function(df, shifting_table, extra_lengths = c(22:24, 34:40), strict.fft = TRUE,
                               BPPARAM = BiocParallel::SerialParam()) {
  files <- names(shifting_table)
  stopifnot(length(files) == nrow(df), !anyNA(files))
  ext <- BiocParallel::bplapply(seq_along(files), function(i) {
    tab <- as.data.table(shifting_table[[i]])
    test <- setdiff(extra_lengths, tab$fraction)
    if (!length(test)) return(list(tab = tab, added = NULL))
    reads <- ORFik::fimport(files[i])
    found <- tryCatch(as.data.table(ORFik::detectRibosomeShifts(reads, df, accepted.lengths = test, must.be.periodic = TRUE,
                                                                strict.fft = strict.fft, verbose = FALSE)),
                      error = function(e) data.table())
    if (!nrow(found)) return(list(tab = tab, added = NULL))
    found <- found[fraction %in% test, .(fraction, offsets_start)]
    list(tab = rbind(tab, found, fill = TRUE)[order(fraction)], added = data.table(sample = files[i], length = found$fraction,
                                                                                   offset = found$offsets_start))
  }, BPPARAM = BPPARAM)
  out <- lapply(ext, function(x) as.data.frame(x$tab))
  names(out) <- files
  setattr(out, "added", rbindlist(lapply(ext, `[[`, "added")))
  out
}

#' Re-P-shift a study with its existing shift table extended by extra periodic lengths
#'
#' Loads the study's current shift table (\code{shifting_table.rds}), extends it with
#' \code{extend_shift_table()} and re-shifts all samples with the extended table (ORFik
#' \code{shiftFootprintsByExperiment(shift.list = ...)}). With \code{out_dir} = NULL the official
#' \code{pshifted/} folder of the experiment's current mapper mode is written (production; reset downstream flags
#' from pshift afterwards via the normal pipeline); otherwise a scratch folder (trailing slash added).
#' @param exp experiment name, or an ORFik experiment (keeps its \code{uniqueMappers} setting)
#' @param extra_lengths read lengths to add when periodic
#' @param out_dir NULL (official pshifted folder) or a scratch directory
#' @param BPPARAM BiocParallel parameters
#' @return the added rows (sample, length, offset), invisibly
#' @export
pshift_add_extra_lengths <- function(exp, extra_lengths = c(22:24, 34:40), out_dir = NULL,
                                     BPPARAM = BiocParallel::SerialParam()) {
  df <- if (is(exp, "experiment")) exp else read.experiment(exp, validate = FALSE, output.env = new.env())
  shifting_table <- shifts_load_safe(df, TRUE)
  if (is.null(shifting_table) || length(shifting_table) != nrow(df) ||
      any(vapply(shifting_table, function(x) is.null(x) || !NROW(x), logical(1))))
    stop(name(df), ": no complete existing shift table; run the normal pshift step first")
  ext <- extend_shift_table(df, shifting_table, extra_lengths, BPPARAM = BPPARAM)
  added <- attr(ext, "added")
  message(name(df), ": added ", NROW(added), " (length, offset) rows over ", nrow(df), " samples")
  if (!NROW(added)) return(invisible(added))
  out <- if (is.null(out_dir)) paste0(sub("/+$", "", libFolder(df)), "/pshifted/") else paste0(sub("/+$", "", out_dir), "/")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ORFik::shiftFootprintsByExperiment(df, out.dir = out, output_format = "ofst", shift.list = ext,
                                     accepted.lengths = sort(unique(c(unlist(lapply(ext, `[[`, "fraction"))))),
                                     BPPARAM = BPPARAM)
  if (is.null(out_dir)) {
    pshift_record_applied_offsets(df)
    save_pshifted_length_distributions(df)
  }
  invisible(added)
}

#' Opt-in pipeline hook: extend a freshly shifted study with extra periodic lengths
#'
#' Called by \code{pipeline_pshift()} after a successful shift when \code{config$pshift_extra_lengths} is set
#' (default NULL = unchanged behaviour). Additive: existing offsets are kept, samples only gain reads.
#' A failure returns the error so the pshifted flag is not set and the step is retried.
#' @noRd
pshift_extra_lengths_hook <- function(df, res, extra_lengths, BPPARAM) {
  if (is.null(extra_lengths) || !length(extra_lengths) || inherits(res, "error")) return(res)
  tryCatch({ pshift_add_extra_lengths(df, extra_lengths, BPPARAM = BPPARAM); res },
           error = function(e) { message(name(df), ": extra-length extension failed: ", conditionMessage(e)); e })
}
