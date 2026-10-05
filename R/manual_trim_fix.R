# Scriptable per-sample trim/barcode fix, outside the interactive flow.
# Companion to R/manual_reshift.R (same audit-log/flag-clearing spirit,
# for the trim/adapter/barcode stage instead of P-shift). Formalizes
# writing barcodes_manual.csv/adapters_manual.csv for one sample plus
# the per-sample-marker + experiment-flag clearing a real retrim needs
# -- see this session's redetect_barcode_params.R (aa_fix_scripts/) for
# the companion read-only "what are the correct parameters" step this
# is meant to apply once you've confirmed them.

#' Apply a corrected adapter/barcode parameter set to one sample
#'
#' Writes (or updates) this one sample's row in \code{barcodes_manual.csv}
#' and/or \code{adapters_manual.csv} -- upserting: any other sample's
#' existing row in those files is preserved, see
#' \code{\link{upsert_manual_override_row}} -- then clears this
#' sample's per-sample markers for every per-sample-tracked step
#' (\code{trim}, \code{collapsed}, \code{aligned}, \code{ofst},
#' \code{ofst_unique}, \code{covrle}, \code{bigwig} -- see
#' \code{R/pipeline_sample_flags.R}), and clears the EXPERIMENT-level
#' flags from \code{trim} onward (same pattern
#' \code{\link{manual_reshift}} already uses for P-shift) so the next
#' regular pipeline run naturally re-enters every stage.
#'
#' Per-sample-tracked stages then skip every already-done sibling and
#' reprocess only this one sample (verified: \code{barcode_detector_single()}
#' reads the manual-override file before its own early-exit check and
#' returns the manually-specified sizes regardless of which internal
#' branch fires, so this does not depend on that function's own
#' behavior being bug-free -- see the collapsed-fasta/early-exit
#' investigation this session). Stages with no per-sample tracking
#' (\code{pshift}, \code{pcounts}, \code{merged_lib}, etc) have no
#' finer granularity to use and necessarily redo the whole study, since
#' they build one object across every sample at once -- this is
#' expected, not a bug in this function.
#' @param exp character, experiment name
#' @param run character, the one sample/run id to fix
#' @param config the mNGSp config object
#' @param barcode5p_size,barcode3p_size numeric, corrected barcode
#' sizes. Both or neither -- leaves \code{barcodes_manual.csv} untouched
#' for this sample if both are NULL (default)
#' @param adapter character, corrected adapter sequence. NULL (default)
#' leaves \code{adapters_manual.csv} untouched for this sample
#' @param note character, free-text reason, appended to the audit log
#' (see \code{\link{log_trim_fix}})
#' @return invisible(NULL)
#' @export
apply_trim_fix_to_sample <- function(exp, run, config, barcode5p_size = NULL,
                                     barcode3p_size = NULL, adapter = NULL, note = "") {
  stopifnot(is.null(barcode5p_size) == is.null(barcode3p_size))
  stopifnot(!is.null(barcode5p_size) || !is.null(adapter))

  df <- read.experiment(exp, validate = FALSE)
  trimmed_dir <- file.path(bam_dir_from_df(df), "trim")
  dir.create(trimmed_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(barcode5p_size)) {
    upsert_manual_override_row(
      file.path(trimmed_dir, "barcodes_manual.csv"),
      data.table::data.table(Run = run, barcode5p_size = barcode5p_size,
                             barcode3p_size = barcode3p_size))
  }
  if (!is.null(adapter)) {
    upsert_manual_override_row(
      file.path(trimmed_dir, "adapters_manual.csv"),
      data.table::data.table(Run = run, adapter = adapter))
  }

  per_sample_steps <- intersect(
    c("trim", "collapsed", "aligned", "ofst", "ofst_unique", "covrle", "bigwig"),
    names(config$flag))
  for (step_id in per_sample_steps) {
    marker <- file.path(sample_flag_dir(config, step_id, exp), paste0(run, ".rds"))
    if (file.exists(marker)) file.remove(marker)
  }

  downstream_steps <- names(config$flag)
  downstream_steps <- downstream_steps[which(downstream_steps == "trim"):length(downstream_steps)]
  remove_flag_all_exp(config, steps = downstream_steps, exps = exp)

  log_trim_fix(exp, run, barcode5p_size, barcode3p_size, adapter, note, config)
  invisible(NULL)
}

#' Upsert one row into a manual-override CSV by Run
#'
#' \code{barcodes_manual_assign_table()}/\code{adapters_manual_assign_table()}
#' (\code{R/fastq_helpers.R}) always overwrite the whole file with
#' exactly the runs given -- fine for their own use (assigning the same
#' value to every run in a study), wrong here (fixing one sample must
#' not silently drop every other sample's existing override row).
#' @param path character, the manual-override CSV path
#' @param new_row data.table, one row, with a \code{Run} column
#' @return invisible(NULL)
#' @noRd
upsert_manual_override_row <- function(path, new_row) {
  combined <- if (file.exists(path)) {
    existing <- data.table::fread(path)
    data.table::rbindlist(list(existing[Run != new_row$Run], new_row), fill = TRUE)
  } else new_row
  data.table::fwrite(combined, path)
  invisible(NULL)
}

#' Append one row to this project's trim-fix audit log
#' @inheritParams apply_trim_fix_to_sample
#' @return invisible(NULL)
#' @noRd
log_trim_fix <- function(exp, run, barcode5p_size, barcode3p_size, adapter, note, config) {
  log_path <- file.path(config$project, "manual_trim_fix_log.csv")
  row <- data.table::data.table(
    timestamp = as.character(Sys.time()), exp = exp, run = run,
    barcode5p_size = if (is.null(barcode5p_size)) NA_real_ else barcode5p_size,
    barcode3p_size = if (is.null(barcode3p_size)) NA_real_ else barcode3p_size,
    adapter = if (is.null(adapter)) NA_character_ else adapter,
    note = note
  )
  dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(row, log_path, append = file.exists(log_path))
  invisible(NULL)
}
