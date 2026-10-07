# Per-sample progress markers.
#
# Parallel to, but distinct from, the experiment-level flags in
# R/pipeline_flags.R. Existing flags only track completion at (study
# accession x organism) = "experiment" granularity; these give finer
# sample/run-level detail for the stages that process samples in a plain
# in-process loop (fetch, trim, align, ofst, covrle, bigwig -- the last
# three via convert_per_sample(), R/pipeline_preset_steps_sub.R), so the
# checklist in R/pipeline_checklist.R can show "M/K samples done" for
# whichever experiment a stage is currently working on, AND (see
# samples_done()/sample_flag_values() below) so a resumed run can skip
# samples already completed in an earlier, interrupted attempt instead of
# redoing the whole experiment from scratch.
#
# Not every stage can offer this: stages that hand off to ORFik's own
# internal BiocParallel dispatch (pshift, valid_pshift) or that need every
# sample together in one call (pcounts' countTable_regions(), which builds
# one SummarizedExperiment across all samples at once; contam, which
# processes a whole folder in one call) have no per-sample loop to hook
# into, and simply have no markers -- the checklist falls back to
# study-level counts only for those, and resume for those stages already
# works at the existing experiment-level flag granularity (no per-sample
# resume needed there). All functions here are internal (unexported),
# same as the existing flag primitives in pipeline_flags.R.
#
# The all-mappers and split-unique-mappers passes of ofst/covrle/bigwig
# are tracked under separate step ids ("ofst" vs "ofst_unique", etc.) --
# a single experiment can have one pass fully done and the other only
# partway through, and each needs its own independent resume point.

#' Directory holding one experiment's per-sample markers for one step
#' @param config the mNGSp config object
#' @param step_id character, e.g. "fetch", "trim", or "aligned"
#' @param experiment character, the "<accession>-<assembly_name>" string
#' @return character, a directory path (not guaranteed to exist yet)
#' @noRd
sample_flag_dir <- function(config, step_id, experiment)
  file.path(config$project, "sample_flags", step_id, experiment)

#' Mark one sample done for one experiment/step
#' @inheritParams sample_flag_dir
#' @param run character, the sample/run id (e.g. SRR accession)
#' @param value what to store in the marker. Usually just \code{TRUE}, but
#' \code{pipeline_trim()} stores the actual per-run \code{barcode_dt} row
#' here instead, so a resumed run can reconstruct the full barcode table
#' across both previously-done and newly-done samples (see
#' \code{\link{sample_flag_values}}).
#' @return invisible(NULL)
#' @noRd
set_sample_flag <- function(config, step_id, experiment, run, value = TRUE) {
  d <- sample_flag_dir(config, step_id, experiment)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  saveRDS(value, file.path(d, paste0(run, ".rds")))
  invisible(NULL)
}

#' Explicit, opt-in "force a clean restart" for one experiment/step
#'
#' Clears sample markers so the next attempt redoes every sample,
#' regardless of what an earlier attempt completed. NOT called implicitly
#' on every stage entry (see \code{\link{samples_done}} for why): a
#' resumed run should skip already-done samples, not silently wipe the
#' evidence of what was already done. Call this deliberately (e.g.
#' alongside \code{remove_flag_all_exp_from()}, pipeline_flags.R) when a
#' real from-scratch reprocess of a step is wanted.
#' @inheritParams sample_flag_dir
#' @return invisible(NULL)
#' @noRd
reset_sample_flags <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (dir.exists(d)) unlink(d, recursive = TRUE)
  invisible(NULL)
}

#' Number of samples marked done for this experiment/step
#' @inheritParams sample_flag_dir
#' @return integer, 0L if the marker directory doesn't exist yet
#' @noRd
n_samples_done <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (!dir.exists(d)) return(0L)
  length(list.files(d, pattern = "\\.rds$"))
}

#' Run ids already marked done for this experiment/step
#'
#' The actual resume-support lookup: which samples can be skipped on a
#' resumed run. One \code{list.files()} call, cheap.
#' @inheritParams sample_flag_dir
#' @return character vector of run ids, \code{character()} if the marker
#' directory doesn't exist yet
#' @noRd
samples_done <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (!dir.exists(d)) return(character())
  sub("\\.rds$", "", list.files(d, pattern = "\\.rds$"))
}

#' Read back every marker's stored value for this experiment/step
#'
#' E.g. to reconstruct a full per-sample table (\code{barcode_dt} in
#' \code{pipeline_trim()}) across a run that resumed partway through --
#' some markers were written in an earlier, interrupted attempt, others
#' just now, and this reads all of them uniformly.
#' @inheritParams sample_flag_dir
#' @return list of stored marker values, named by run id (the marker
#' filename without its \code{.rds} extension), \code{list()} if the
#' marker directory doesn't exist yet
#' @noRd
sample_flag_values <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (!dir.exists(d)) return(list())
  files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  stats::setNames(lapply(files, readRDS), sub("\\.rds$", "", basename(files)))
}
