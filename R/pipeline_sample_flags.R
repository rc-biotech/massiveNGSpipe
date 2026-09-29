# Per-sample progress markers.
#
# Parallel to, but distinct from, the experiment-level flags in
# R/pipeline_flags.R. Existing flags only track completion at (study
# accession x organism) = "experiment" granularity; these give finer
# sample/run-level detail for the stages that process samples in a plain
# in-process loop (currently: fetch, trim, align), so the checklist in
# R/pipeline_checklist.R can show "M/K samples done" for whichever
# experiment a stage is currently working on, AND (see samples_done()/
# sample_flag_values() below) so a resumed run can skip samples already
# completed in an earlier, interrupted attempt instead of redoing the
# whole experiment from scratch.
#
# Not every stage can offer this: stages that hand off to ORFik's own
# internal BiocParallel dispatch (pshift, valid_pshift, pcounts) or that
# process a whole folder in one call (contam) have no per-sample loop to
# hook into, and simply have no markers -- the checklist falls back to
# study-level counts only for those, and resume for those stages already
# works at the existing experiment-level flag granularity (no per-sample
# resume needed there). All functions here are internal (unexported),
# same as the existing flag primitives in pipeline_flags.R.

sample_flag_dir <- function(config, step_id, experiment)
  file.path(config$project, "sample_flags", step_id, experiment)

# config: the mNGSp config object
# step_id: character, e.g. "trim" or "aligned"
# experiment: character, the "<accession>-<assembly_name>" string
# run: character, the sample/run id (e.g. SRR accession)
# value: what to store in the marker. Usually just TRUE, but
# pipeline_trim() stores the actual per-run barcode_dt row here instead,
# so that a resumed run can reconstruct the full barcode table across
# both previously-done and newly-done samples (see sample_flag_values()).
set_sample_flag <- function(config, step_id, experiment, run, value = TRUE) {
  d <- sample_flag_dir(config, step_id, experiment)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  saveRDS(value, file.path(d, paste0(run, ".rds")))
  invisible(NULL)
}

# Explicit, opt-in "force a clean restart" path -- clears sample markers
# for one experiment so the next attempt redoes every sample, regardless
# of what an earlier attempt completed. NOT called implicitly anymore on
# every stage entry (see samples_done() below for why): a resumed run
# should skip already-done samples, not silently wipe the evidence of
# what was already done. Call this deliberately (e.g. alongside
# remove_flag_all_exp_from(), pipeline_flags.R:248) when a real from-
# scratch reprocess of a step is wanted.
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

# Character vector of run ids already marked done for this experiment/step
# -- the actual resume-support lookup. One list.files() call, cheap.
samples_done <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (!dir.exists(d)) return(character())
  sub("\\.rds$", "", list.files(d, pattern = "\\.rds$"))
}

# Read back every marker's stored value for this experiment/step, e.g. to
# reconstruct a full per-sample table (barcode_dt in pipeline_trim())
# across a run that resumed partway through -- some markers were written
# in an earlier, interrupted attempt, others just now, and this reads all
# of them uniformly.
sample_flag_values <- function(config, step_id, experiment) {
  d <- sample_flag_dir(config, step_id, experiment)
  if (!dir.exists(d)) return(list())
  lapply(list.files(d, pattern = "\\.rds$", full.names = TRUE), readRDS)
}
