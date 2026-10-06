# End-to-end, single-call workflow for the "redetect -> fix -> rerun ->
# verify" cycle this session built up by hand, one bash/R snippet at a
# time, across several real barcode fixes (PRJNA926112, PRJNA770650,
# PRJEB50305). Collapsing that into a few package functions means one
# Rscript launch + one background wait per candidate instead of several
# separate round trips, and a small status file to read back instead of
# tailing a multi-hundred-line pipeline log each time.
#
# Companion to apply_trim_fix_to_sample() (R/manual_trim_fix.R), which
# this builds on rather than duplicates.

#' Is anything actively touching this experiment's raw/processed files
#' right now?
#'
#' A cheap, read-only safety check: greps the process table for a
#' fastp/STAR invocation whose command line mentions this experiment's
#' name. Meant to run BEFORE \code{\link{apply_trim_fix_and_rerun}} (and
#' is called inside it automatically) -- fixing and rerunning a sample
#' while something else is genuinely mid-flight on the SAME study risks
#' exactly the kind of file race this session hit more than once by
#' accident.
#' @param exp character, experiment name (e.g. "PRJNA926112-homo_sapiens")
#' @return logical
#' @export
sample_currently_processing <- function(exp) {
  procs <- tryCatch(system2("pgrep", c("-af", "fastp|STAR"), stdout = TRUE, stderr = FALSE),
                    warning = function(w) character(), error = function(e) character())
  any(grepl(exp, procs, fixed = TRUE))
}

#' Most recent modification time of an experiment's real files, with
#' known AI/pipeline-diagnostic noise excluded
#'
#' "Real" here excludes files this package's own diagnostic/backfill
#' passes touch independently of any actual reprocessing (read-length
#' distributions, reshift-applied markers, the shared summary-statistics
#' CSVs, and any bare .rds flag/marker file) -- without excluding these,
#' a routine diagnostic backfill makes a study look "recently touched"
#' when nothing about its real data changed. See the file-mtime safety
#' checks done by hand earlier this session for the exact same
#' exclusion list, now centralized here.
#' @param exp character, experiment name
#' @param config the mNGSp config object
#' @return POSIXct, or NA if no files are found
#' @export
last_real_touch <- function(exp, config = pipeline_config()) {
  study_accession <- sub("-.*", "", exp)
  raw_dir <- file.path(config$config["fastq"], exp)
  bam_dir <- file.path(config$config["bam"], exp)
  dirs <- c(raw_dir, bam_dir)[dir.exists(c(raw_dir, bam_dir))]
  if (length(dirs) == 0) return(as.POSIXct(NA))

  files <- unlist(lapply(dirs, function(d) list.files(d, recursive = TRUE, full.names = TRUE)))
  noise <- c("read_length_distribution", "auto_reshift_applied\\.rds$", "shifting_table\\.rds$",
            "summary_statistics", "\\.rds$")
  for (pat in noise) files <- files[!grepl(pat, files)]
  if (length(files) == 0) return(as.POSIXct(NA))
  max(file.info(files)$mtime, na.rm = TRUE)
}

#' Rank candidate samples for a barcode/trim fix, excluding
#' already-handled and currently-busy studies
#'
#' Automates the per-candidate triage this session did by hand: reads
#' \code{find_barcode_adapter_outliers.R}'s own per-sample output,
#' restricts to barcode-detection outliers, keeps one (the flagged)
#' sample per study, computes what fraction of its siblings DID detect
#' a barcode (the strongest confidence signal observed in practice --
#' a near-unanimous sibling majority with one holdout is a much
#' stronger signal than a borderline length difference alone), excludes
#' any (exp, run) already present in \code{manual_trim_fix_log.csv}
#' (apply_trim_fix_to_sample()'s own audit log -- so a sample already
#' fixed, or already investigated and intentionally declined, is never
#' re-suggested), and excludes anything currently mid-processing
#' (\code{\link{sample_currently_processing}}) or touched more recently
#' than \code{min_days_untouched} (\code{\link{last_real_touch}}).
#'
#' Read-only -- does not redetect, fix, or touch anything. Pass the
#' returned candidates to \code{\link{redetect_barcode_for_sample}} to
#' get real evidence before deciding, and always sanity-check a
#' candidate's FULL sibling table (not just this one summary row)
#' before applying a fix -- a majority pattern can still be
#' misleading when barcode sizes vary by sub-batch within one study
#' (confirmed live, PRJNA770650-homo_sapiens: a uniform-looking 51/52
#' majority actually split into at least 3 different per-batch barcode
#' sizes).
#' @param config the mNGSp config object
#' @param outliers_per_sample_path character, path to
#' \code{find_barcode_adapter_outliers.R}'s per-sample output CSV
#' @param min_majority_frac numeric, default 0.5. Minimum fraction of
#' siblings that must have a detected barcode for a study to qualify
#' @param min_days_untouched numeric, default 0 (no staleness
#' requirement). Set > 0 to additionally require the study's real files
#' (see \code{\link{last_real_touch}}) to predate now by at least this
#' many days -- useful when avoiding interference with other concurrent
#' work matters, not just avoiding already-fixed samples
#' @param n integer, default 10. Max rows to return, ranked best first
#' @return data.table(exp, run, study_accession, organism,
#' majority_frac, n_study_samples, n_study_barcode_true,
#' trim_mean_length, study_median_len, days_untouched), sorted by
#' majority_frac descending
#' @export
barcode_fix_candidates <- function(config, outliers_per_sample_path, min_majority_frac = 0.5,
                                   min_days_untouched = 0, n = 10) {
  dt <- data.table::fread(outliers_per_sample_path)
  dt <- dt[flag_reason %in% c("barcode_outlier", "barcode_and_length_outlier")]
  dt[, majority_frac := n_study_barcode_true / n_study_samples]
  dt <- unique(dt[order(-majority_frac)], by = "study_accession")
  dt <- dt[majority_frac > min_majority_frac]
  dt[, exp := paste0(study_accession, "-", gsub(" ", "_", tolower(ScientificName)))]

  already_handled_path <- file.path(config$project, "manual_trim_fix_log.csv")
  if (file.exists(already_handled_path)) {
    handled <- data.table::fread(already_handled_path)[, .(exp, run)]
    dt <- dt[!handled, on = c(exp = "exp", raw_library = "run")]
  }

  dt <- dt[order(-majority_frac)]
  if (nrow(dt) == 0) return(dt)
  dt <- dt[seq_len(min(nrow(dt), max(n * 4, 20)))] # bound how many get an mtime/process check below

  dt[, currently_active := vapply(exp, sample_currently_processing, logical(1))]
  dt <- dt[currently_active == FALSE]
  dt[, last_touch := vapply(exp, function(e) as.character(last_real_touch(e, config)), character(1))]
  dt[, days_untouched := round(as.numeric(difftime(Sys.time(), as.POSIXct(last_touch), units = "days")), 1)]
  dt <- dt[is.na(days_untouched) | days_untouched >= min_days_untouched]

  dt <- dt[order(-majority_frac)][seq_len(min(nrow(dt), n))]
  dt[, .(exp, run = raw_library, study_accession, organism = ScientificName, majority_frac,
        n_study_samples, n_study_barcode_true, trim_mean_length, study_median_len, days_untouched)]
}

#' Force-redetect one sample's adapter/barcode parameters, bypassing
#' \code{barcode_detector_single()}'s own early-exit heuristic
#'
#' Thin wrapper around the exact call this session used by hand for
#' every real fix: \code{check_at_mean_size = 0} guarantees real
#' change-point detection always runs, regardless of how short the
#' sample's existing trimmed length looks (the root cause behind every
#' confirmed fix so far -- see \code{R/fastq_helpers.R}'s own early-exit
#' documentation). Downloads the sample's raw fastq if needed
#' (\code{redownload_raw_if_needed = TRUE}) -- budget a few minutes,
#' this is a real network + fastp call, not a metadata lookup.
#'
#' This is evidence-gathering only -- it does not decide whether the
#' result is actually correct for this sample (see
#' \code{\link{barcode_fix_candidates}}'s own warning about per-batch
#' variation) and does not write anything. Read the result's
#' \code{consensus_string_5p}/\code{consensus_string_3p} against the
#' sample's siblings' own \code{adapter_barcode_table.csv} before
#' trusting it.
#' @param exp character, experiment name
#' @param run character, the run id to redetect
#' @param config the mNGSp config object
#' @return data.table, one row (same shape as
#' \code{adapter_barcode_table.csv})
#' @export
redetect_barcode_for_sample <- function(exp, run, config = pipeline_config()) {
  final_list <- data.table::fread(config$complete_metadata)
  study_org <- final_list[Run == run]
  if (nrow(study_org) == 0) stop("Run not found in ", config$complete_metadata, ": ", run)
  df <- read.experiment(exp, validate = FALSE)
  bam_root <- bam_dir_from_df(df)
  fastq_dir <- file.path(config$config["fastq"], exp)
  trimmed_dir <- file.path(bam_root, "trim")
  barcode_detector_single(study_org[1, .(Run, LibraryLayout)], fastq_dir, bam_root, trimmed_dir,
                          redownload_raw_if_needed = TRUE, check_at_mean_size = 0)
}

#' Apply a trim/barcode fix to one sample and immediately rerun the
#' pipeline for just that study, end to end, in one call
#'
#' Combines, in order: a pre-flight check
#' (\code{\link{sample_currently_processing}}, aborts rather than
#' racing with something else already touching this study),
#' \code{\link{apply_trim_fix_to_sample}} (writes the manual override,
#' backfills sibling markers, clears flags), building a config/pipelines
#' pair scoped to ONLY this one study's accession (never the full
#' catalog), \code{run_pipeline()}, and a sibling-file-integrity check
#' comparing each OTHER sample's own files' modification times before
#' and after -- any change there is a bug (this exact bug has happened
#' twice this session: \code{pipeline_align()}'s old
#' \code{delete_collapsed_files} branch, and \code{pipeline_cleanup()}'s
#' old self-move-deletes-file bug) and should stop you from trusting the
#' run, not just a cosmetic warning.
#'
#' Writes a small result list to
#' \code{<config$project>/ai_fix_status/<exp>.rds} as well as returning
#' it, so a later check only needs one cheap \code{readRDS()} instead of
#' re-deriving flags or tailing a pipeline log.
#' @inheritParams apply_trim_fix_to_sample
#' @param wait numeric, passed to \code{run_pipeline()}'s own poll
#' interval, default 10
#' @return invisible list(exp, run, success, elapsed_mins, final_flags,
#' touched_sibling_files, note) -- \code{success} is TRUE only when
#' every flag is done AND \code{touched_sibling_files} is empty
#' @export
apply_trim_fix_and_rerun <- function(exp, run, config, barcode5p_size = NULL,
                                     barcode3p_size = NULL, adapter = NULL, note = "", wait = 10) {
  if (sample_currently_processing(exp))
    stop("Refusing to proceed: ", exp, " is currently being processed by something else. ",
        "Re-check with sample_currently_processing() once it's done.")

  # Named distinctly from the "study_accession" COLUMN used below --
  # data.table resolves an i-expression's bare column-name references
  # against its own columns first, so `dt[study_accession == study_accession]`
  # with a same-named local variable would silently compare the column
  # to itself (always TRUE for every row), not to this value.
  target_accession <- sub("-.*", "", exp)
  df_current <- read.experiment(exp, validate = FALSE)
  bam_dir <- bam_dir_from_df(df_current)
  all_runs <- ORFik::runIDs(df_current)
  siblings <- setdiff(all_runs, run)
  before <- sibling_file_snapshot(bam_dir, siblings)

  apply_trim_fix_to_sample(exp, run, config, barcode5p_size = barcode5p_size,
                           barcode3p_size = barcode3p_size, adapter = adapter, note = note)

  scratch_meta <- data.table::fread(config$complete_metadata)[study_accession == target_accession]
  scratch_path <- tempfile(fileext = ".csv")
  data.table::fwrite(scratch_meta, scratch_path)

  t0 <- Sys.time()
  pipelines <- pipeline_init_all(config, complete_metadata = scratch_path,
                                 only_complete_genomes = TRUE, gene_symbols = FALSE,
                                 show_status_per_exp = FALSE, verbose_load_annotation = FALSE)
  run_pipeline(pipelines, config, wait = wait)
  elapsed_mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)

  after <- sibling_file_snapshot(bam_dir, siblings)
  touched_sibling_files <- changed_files(before, after)

  final_flags <- stats::setNames(vapply(names(config$flag), function(s) step_is_done(config, s, exp), logical(1)),
                                 names(config$flag))
  result <- list(exp = exp, run = run, success = all(final_flags) && length(touched_sibling_files) == 0,
                 elapsed_mins = elapsed_mins, final_flags = final_flags,
                 touched_sibling_files = touched_sibling_files, note = note)

  status_dir <- file.path(config$project, "ai_fix_status")
  dir.create(status_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(result, file.path(status_dir, paste0(exp, ".rds")))
  invisible(result)
}

#' Named mtime vector for every sibling run's own files under an
#' experiment's bam directory
#'
#' Excludes files under directories this package's own pshift/pcounts
#' steps legitimately rewrite for EVERY sample on EVERY run, by design
#' -- they have no per-sample resume (ORFik's own internal per-experiment
#' dispatch for pshift; countTable_regions() needs every sample together
#' for pcounts), documented in \code{convert_per_sample()}'s own roxygen
#' and confirmed live: fixing one sample and rerunning correctly
#' regenerates `pshifted/*_pshifted.ofst`, `read_length_distribution/`
#' under pshifted, and `QC_STATS/` (frame tables, count tables, plots)
#' for every sibling too. That is expected, not a bug -- counting it as
#' a "touched sibling file" would flag every single real fix as a
#' false-positive failure (confirmed live, PRJNA1071171-homo_sapiens,
#' 2026-10-06). What should still be scrutinized here is anything with
#' real per-sample resume (BAM, collapsed fasta, trimmed fastq, ofst,
#' covrle, bigwig) -- those must never change for an already-done
#' sibling.
#' @param bam_dir character, the experiment's bam directory
#' @param siblings character vector of run ids to match against
#' filenames
#' @return named numeric vector (names are file paths, values are
#' as.numeric(mtime)); empty if no matching files exist
#' @noRd
sibling_file_snapshot <- function(bam_dir, siblings) {
  if (length(siblings) == 0 || !dir.exists(bam_dir)) return(stats::setNames(numeric(), character()))
  files <- list.files(bam_dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl(paste(siblings, collapse = "|"), basename(files))]
  files <- files[!grepl("/pshifted/|/QC_STATS/", files)]
  if (length(files) == 0) return(stats::setNames(numeric(), character()))
  stats::setNames(as.numeric(file.info(files)$mtime), files)
}

#' Which files changed (by mtime, or by appearing/disappearing) between
#' two \code{\link{sibling_file_snapshot}} results
#' @param before,after named numeric vectors from
#' \code{\link{sibling_file_snapshot}}
#' @return character vector of changed file paths, empty if none
#' @noRd
changed_files <- function(before, after) {
  all_paths <- union(names(before), names(after))
  changed <- vapply(all_paths, function(p) {
    # Exact numeric comparison, not all.equal(): mtimes are epoch
    # seconds (~1.7e9), so all.equal()'s default RELATIVE tolerance
    # (~1.5e-8) swamps any real, few-second change -- a file genuinely
    # touched a second apart would still compare as "equal". A path
    # missing from either side (appeared/disappeared) also counts as
    # changed; indexing by name on a plain atomic vector errors (not
    # NULL) for a missing name, so check membership first.
    p_in_before <- p %in% names(before)
    p_in_after <- p %in% names(after)
    if (!p_in_before || !p_in_after) return(TRUE)
    before[[p]] != after[[p]]
  }, logical(1))
  all_paths[changed]
}
