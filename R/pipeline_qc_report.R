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
