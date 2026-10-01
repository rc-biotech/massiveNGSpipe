# Cross-study aggregation of the per-experiment qc_diagnostics.rds records
# (see R/pshift_diagnostics.R). Deliberately a separate file/function from
# status_per_study()/progress_report() (R/pipeline_progress.R) --
# completion-tracking and QC-verdict-tracking are different concerns.

#' `x` unless it's NULL, in which case `na`
#' @noRd
or_na <- function(x, na = NA) if (is.null(x)) na else x

#' Report P-shift QC verdicts + likely failure cause across all studies
#'
#' Reads back every experiment's \code{qc_diagnostics.rds} (written
#' progressively by \code{check_too_few_reads()}/\code{check_alignment_rate()}/
#' \code{check_adapter_barcode_quality()}/\code{periodicity_check_flag()}/
#' \code{classify_bad_shift_cause()} as each pipeline stage completes) into
#' one table for triage. A study that hasn't reached the trim step yet
#' simply has no row here -- this reports what's been recorded, not every
#' study \code{pipelines} contains.
#' @param pipelines the full pipelines list (see \code{pipeline_init_all()})
#' @param config the mNGSp config object
#' @param only_flagged logical, default FALSE. If TRUE, drop rows with no
#' recorded problem at all: \code{periodicity_status} not in
#' \code{c("warning", "no_data")} AND \code{primary_cause} in
#' \code{c(NA, "unknown")}.
#' @return a data.table, one row per experiment with a recorded
#' diagnostics file; zero rows (not an error) if nothing has been
#' recorded yet
#' @export
qc_verdict_report <- function(pipelines, config, only_flagged = FALSE) {
  rows <- list()
  for (pipeline in pipelines) {
    for (organism in names(pipeline$organisms)) {
      conf <- pipeline$organisms[[organism]]$conf
      experiment <- unname(conf["exp"])
      diag <- read_qc_diagnostics(qc_diagnostics_path(conf["bam"]))
      if (length(diag) == 0) next
      rows[[experiment]] <- data.table::data.table(
        experiment = experiment,
        organism = organism,
        too_few_reads = or_na(diag$too_few_reads),
        raw_reads = or_na(diag$raw_reads, NA_real_),
        wrong_organism = or_na(diag$wrong_organism),
        mapping_rate_pct = or_na(diag$mapping_rate_pct, NA_real_),
        bad_adapter_barcode = or_na(diag$bad_adapter_barcode),
        no_adapter_removed_pct = or_na(diag$no_adapter_removed_pct, NA_real_),
        periodicity_status = or_na(diag$periodicity_status, NA_character_),
        primary_cause = or_na(diag$primary_cause, NA_character_),
        last_updated = or_na(diag$last_updated, as.POSIXct(NA))
      )
    }
  }
  result <- data.table::rbindlist(rows, fill = TRUE)
  if (only_flagged && nrow(result) > 0) {
    has_periodicity_problem <- result$periodicity_status %in% c("warning", "no_data")
    has_classified_cause <- !is.na(result$primary_cause) & result$primary_cause != "unknown"
    result <- result[has_periodicity_problem | has_classified_cause]
  }
  result[]
}

#' One study's full P-shift diagnostic picture, read entirely from
#' precomputed files
#'
#' A pure read: never loads an ofst file, never recomputes anything --
#' everything here is already written progressively by the pipeline
#' itself (\code{qc_diagnostics.rds} by \code{R/pshift_diagnostics.R}'s
#' check functions; \code{Ribo_frames_all.csv} by \code{shift_qc()},
#' \code{R/shifting_helpers.R}; the two \code{00_aggregated.csv} length
#' distributions by \code{convert_bam_to_ofst_with_length_dist()}/
#' \code{save_pshifted_length_distributions()},
#' \code{R/pipeline_preset_steps_sub.R}). Cheap enough to call on demand
#' to decide whether and how a study needs reshifting: comparing
#' \code{length_dist_preshift} against \code{accepted_lengths} shows
#' whether footprint sizes genuinely fall outside the configured range
#' (an \code{accepted_lengths_rpf} fix, not a reshift); comparing
#' \code{length_dist_preshift} against \code{length_dist_postshift}
#' shows how much got dropped by that filter; \code{frames} shows
#' whether periodicity is actually poor (an offset fix) once the data
#' itself looks reasonable.
#' @param exp character, experiment name
#' @param config the mNGSp config object
#' @return a list: \code{exp}, \code{qc} (\code{qc_diagnostics.rds}
#' contents, \code{list()} if not yet written), \code{frames}
#' (data.table from \code{Ribo_frames_all.csv}, or \code{NULL}),
#' \code{length_dist_preshift}/\code{length_dist_postshift} (data.table
#' from the respective \code{00_aggregated.csv}, or \code{NULL} if not
#' yet written), \code{accepted_lengths} (\code{config$accepted_lengths_rpf})
#' @export
pshift_triage <- function(exp, config = pipeline_config()) {
  df <- read.experiment(exp, validate = FALSE)
  bam_dir <- bam_dir_from_df(df)

  read_dist_or_null <- function(dir_of_filetype) {
    path <- tryCatch(file.path(dirname(filepath(df, dir_of_filetype)[1]),
                               "read_length_distribution", "00_aggregated.csv"),
                     error = function(e) NA_character_)
    if (is.na(path) || !file.exists(path)) return(NULL)
    data.table::fread(path)
  }

  frames_path <- file.path(QCfolder(df), "Ribo_frames_all.csv")

  list(
    exp = exp,
    qc = read_qc_diagnostics(qc_diagnostics_path(bam_dir)),
    frames = if (file.exists(frames_path)) data.table::fread(frames_path) else NULL,
    length_dist_preshift = read_dist_or_null("ofst"),
    length_dist_postshift = read_dist_or_null("pshifted"),
    accepted_lengths = config$accepted_lengths_rpf
  )
}
