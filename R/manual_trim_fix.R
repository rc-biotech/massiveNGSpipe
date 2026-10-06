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
#'
#' Before clearing anything, this also backfills a sibling marker for
#' every OTHER sample in the experiment that doesn't already have one,
#' for any per-sample step whose experiment-level flag is currently
#' TRUE. This matters for any study processed before the per-sample-
#' resume feature existed: such a study has an experiment-level "done"
#' flag but no per-sample markers at all, for any sample -- without the
#' backfill, \code{samples_done()} would read as empty for every
#' sibling too, and \code{pipeline_trim()}/\code{pipeline_collapse()}
#' would try to reprocess the whole study again, which for trim means
#' re-resolving already-deleted raw files for the untouched siblings
#' and erroring immediately (confirmed live, 2026-10-05,
#' PRJNA926112-homo_sapiens).
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

  # Backfill sibling markers for any step that predates this session's
  # per-sample-resume feature: a study fully processed before per-sample
  # flags existed has an experiment-level "done" flag but ZERO per-sample
  # markers for any of its samples -- the fixed sample included. Without
  # this, clearing just the fixed sample's (nonexistent) marker below
  # leaves samples_done() empty for every sibling too, so pipeline_trim()/
  # pipeline_collapse() treat the WHOLE study as not-yet-done and try to
  # reprocess every sample -- for trim specifically, re-resolving long-
  # since-deleted raw files for already-finished siblings, which errors
  # immediately and then (via run_pipeline()'s per-experiment retry-skip
  # design) looks like the whole pipeline has hung. Confirmed live on
  # PRJNA926112-homo_sapiens (processed 2026-07, before this feature
  # existed). Only backfill while the step's own experiment-level flag is
  # still TRUE -- checked before remove_flag_all_exp() clears it below --
  # since that is the only evidence a sample was actually done.
  # "trim" markers store a one-row data.table (not just TRUE) --
  # pipeline_trim() reconstructs adapter_barcode_table.csv via
  # rbindlist(sample_flag_values(config, "trim", experiment)), which
  # errors ("Item 1 of input is not a data.frame...") on a plain TRUE
  # backfill. Reuse that sibling's own existing row from the
  # already-on-disk table (still correct/historical for the sibling,
  # only the fixed sample's row there is wrong) instead of a dummy value.
  existing_trim_table <- NULL
  trim_table_path <- file.path(trimmed_dir, "adapter_barcode_table.csv")
  if (file.exists(trim_table_path)) existing_trim_table <- data.table::fread(trim_table_path)

  # Unique-mapper variants of ofst/covrle/bigwig are tracked under their
  # own per-sample-flag namespace (convert_per_sample(), called twice --
  # once per mapper mode -- inside pipeline_create_ofst()/
  # pipeline_create_covrle()/pipeline_create_bigwig(), R/pipeline_preset_steps_sub.R)
  # but have NO experiment-level flag of their own: both passes share
  # the SAME parent flag ("ofst"/"covrle"/"bigwig"), set once after both
  # finish. Gate their backfill on that shared parent flag instead of
  # (nonexistent) step_is_done(config, "ofst_unique", ...).
  unique_variant_parents <- c(ofst_unique = "ofst", covrle_unique = "covrle", bigwig_unique = "bigwig")
  all_per_sample_steps <- c(per_sample_steps, names(unique_variant_parents))

  all_runs <- ORFik::runIDs(df)
  sibling_runs <- setdiff(all_runs, run)
  for (step_id in all_per_sample_steps) {
    gate_step <- if (step_id %in% names(unique_variant_parents)) unique_variant_parents[[step_id]] else step_id
    if (!step_is_done(config, gate_step, exp)) next
    to_backfill <- setdiff(sibling_runs, samples_done(config, step_id, exp))
    for (sibling in to_backfill) {
      value <- TRUE
      if (step_id == "trim") {
        # Always a data.table-compatible value for "trim" (never plain
        # TRUE), regardless of whether this sibling has a row in the
        # existing table: an experiment can predate the per-sample-
        # resume feature badly enough that adapter_barcode_table.csv
        # itself is missing most siblings' rows (confirmed live,
        # PRJNA770650-homo_sapiens: 9 of 52 rows present) -- the
        # rbindlist() crash this guards against would otherwise still
        # happen for every sibling without a row, not just when the
        # table is entirely absent.
        value <- if (!is.null(existing_trim_table) && sibling %in% existing_trim_table$id)
          existing_trim_table[id == sibling] else data.table::data.table(id = sibling)
      }
      set_sample_flag(config, step_id, exp, sibling, value = value)
    }
  }

  for (step_id in all_per_sample_steps) {
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
