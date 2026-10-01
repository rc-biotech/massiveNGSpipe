# Scriptable (non-Shiny) manual P-site reshift. Formalizes the pattern
# already used by hand in shared_scripts/.../mNGSp_manual_reshift.R and
# R/shiny_app_reshifting.R's own save handler, so a reshift done this way
# is indistinguishable from one done interactively -- same shift-table
# shape, same manually_checked_shifts.rds marker -- and the Shiny app,
# progress_report(), and the rest of the pipeline all keep working
# unchanged afterward. The only functions in this package that write/
# mutate real shift data outside the main pipeline flow, hence their own
# file for discoverability/review.

#' Apply an explicit per-length-fraction P-site offset table to one
#' experiment, outside the interactive Shiny editor
#'
#' Does the real shift immediately (unlike the Shiny app's own save
#' handler, which only edits \code{shifting_table.rds} and defers the
#' actual reshift to the next regular pipeline run) -- appropriate for a
#' one-off, explicit, already-decided fix where immediate feedback
#' matters. \code{shifts_save()} still persists the same
#' \code{shifting_table.rds} \code{shifts_load_safe()} reads on a future
#' resume, and the Shiny app would also read, so neither is left stale.
#'
#' Clears \code{"pshifted"} and every downstream flag before re-marking
#' \code{"pshifted"} done, so the regular pipeline naturally re-runs
#' \code{valid_pshift}/\code{merged_lib}/etc. on its own next pass --
#' same intent as the Shiny app's own
#' \code{remove_flag_all_exp(steps = steps_to_remove, ...)}, computed
#' dynamically here from \code{names(config$flag)} rather than a
#' hardcoded step-count, so it stays correct regardless of preset.
#' @param exp character, experiment name
#' @param offsets data.table(fraction, offsets_start) -- the exact shape
#' \code{ORFik::shiftFootprintsByExperiment()} already expects (same
#' shape \code{R/shiny_app_reshifting.R}'s own save handler builds)
#' @param config the mNGSp config object
#' @param note character, free-text reason, appended to the audit log
#' (see \code{\link{log_manual_reshift}})
#' @return invisible(NULL)
#' @export
manual_reshift <- function(exp, offsets, config = pipeline_config(), note = "") {
  stopifnot(is.data.frame(offsets) && all(c("fraction", "offsets_start") %in% names(offsets)))
  df <- read.experiment(exp, validate = FALSE)
  ofst_files <- filepath(df, "ofst")
  shifts <- stats::setNames(rep(list(offsets), length(ofst_files)), ofst_files)

  ORFik::shiftFootprintsByExperiment(df, output_format = "ofst", shift.list = shifts)
  ORFik::shifts_save(shifts, file.path(libFolder(df), "pshifted"))
  ORFik::shiftPlots(df, output = "auto", plot.ext = ".png")
  shift_qc(df, max_no_adapter_removed_pct = config$max_no_adapter_removed_pct)
  save_pshifted_length_distributions(df)
  dir.create(QCfolder(df), showWarnings = FALSE, recursive = TRUE)
  saveRDS(TRUE, file.path(QCfolder(df), "manually_checked_shifts.rds"))
  log_manual_reshift(exp, offsets, note, config)

  downstream_steps <- names(config$flag)
  downstream_steps <- downstream_steps[which(downstream_steps == "pshifted"):length(downstream_steps)]
  remove_flag_all_exp(config, steps = downstream_steps, exps = exp)
  set_flag(config, "pshifted", exp)
  invisible(NULL)
}

#' Append one row to this project's manual-reshift audit log
#'
#' A small, append-only trail of every \code{\link{manual_reshift}} call
#' -- matters more once something other than an operator's own hands is
#' calling it.
#' @param exp character, experiment name
#' @param offsets data.table(fraction, offsets_start), as applied
#' @param note character, free-text reason
#' @param config the mNGSp config object (for \code{config$project})
#' @return invisible(NULL)
#' @noRd
log_manual_reshift <- function(exp, offsets, note, config) {
  log_path <- file.path(config$project, "manual_reshift_log.csv")
  row <- data.table::data.table(
    timestamp = as.character(Sys.time()),
    exp = exp,
    offsets = paste(offsets$fraction, offsets$offsets_start, sep = ":", collapse = ";"),
    note = note
  )
  dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(row, log_path, append = file.exists(log_path))
  invisible(NULL)
}
