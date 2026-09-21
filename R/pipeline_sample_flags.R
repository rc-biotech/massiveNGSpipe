# Per-sample progress markers.
#
# Parallel to, but distinct from, the experiment-level flags in
# R/pipeline_flags.R. Existing flags only track completion at (study
# accession x organism) = "experiment" granularity; these give finer
# sample/run-level detail for the stages that process samples in a plain
# in-process loop (currently: trim, align), so the checklist in
# R/pipeline_checklist.R can show "M/K samples done" for whichever
# experiment a stage is currently working on.
#
# Not every stage can offer this: stages that hand off to ORFik's own
# internal BiocParallel dispatch (pshift, valid_pshift, pcounts) or that
# process a whole folder in one call (contam) have no per-sample loop to
# hook into, and simply have no markers -- the checklist falls back to
# study-level counts only for those. All functions here are internal
# (unexported), same as the existing flag primitives in pipeline_flags.R.

sample_flag_dir <- function(config, step_id, experiment)
  file.path(config$project, "sample_flags", step_id, experiment)

# config: the mNGSp config object
# step_id: character, e.g. "trim" or "aligned"
# experiment: character, the "<accession>-<assembly_name>" string
# run: character, the sample/run id (e.g. SRR accession)
set_sample_flag <- function(config, step_id, experiment, run) {
  d <- sample_flag_dir(config, step_id, experiment)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  saveRDS(TRUE, file.path(d, paste0(run, ".rds")))
  invisible(NULL)
}

# Clear sample markers for one experiment before (re)starting a stage.
# A restart today always redoes a whole not-yet-done experiment from
# scratch (no per-sample resume exists), so any markers left over from an
# earlier, interrupted attempt at this experiment must not be shown as
# "already done" for a new attempt.
reset_sample_flags <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (dir.exists(d)) unlink(d, recursive = TRUE)
  invisible(NULL)
}

# Returns the number of samples marked done for this experiment/step.
n_samples_done <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (!dir.exists(d)) return(0L)
  length(list.files(d, pattern = "\\.rds$"))
}
