# Shared fixture builders for massiveNGSpipe unit tests.
#
# These build the minimal-but-structurally-real config/pipelines shapes
# that most massiveNGSpipe orchestration functions expect, backed by real
# tempdir() directories (not mocks) so flag/file-based functions behave
# exactly as they do in production, just against throwaway paths.

#' Minimal real mNGSp config for tests, with real flag directories created
#' under a tempdir() project.
fake_config <- function(project = tempfile("mNGSp_test_"), preset = "RNA-seq",
                        mode = "local", contam = FALSE, session_dir = NULL,
                        extra = list()) {
  flags <- pipeline_flags(project, mode = mode, preset = preset,
                          contam = contam, create_dirs = TRUE)
  flag_steps <- if (preset == "empty") list() else flag_grouping(flags)
  config <- list(
    project = project,
    preset = preset,
    mode = mode,
    flag = flags,
    flag_steps = flag_steps,
    config = c(ref = file.path(project, "ref"), bam = file.path(project, "bam"),
              fastq = file.path(project, "fastq"), exp = file.path(project, "exp")),
    session_dir = session_dir,
    error_dir = NULL,
    discord_webhook = NULL,
    stop_downloading_new_data_at_drive_usage = 92,
    compress_raw_data = TRUE
  )
  # utils::modifyList (not c()) so that names present in `extra` REPLACE
  # the default above rather than being appended as a duplicate --
  # config$x on a list with two same-named elements silently returns
  # only the first one, which would make `extra` overrides no-ops.
  utils::modifyList(config, extra)
}

#' Minimal fake `pipelines` list matching the shape used across the
#' package: named list of accession -> list(accession, organisms, study).
fake_pipelines <- function(accession = "PRJNA000001", organism = "Homo sapiens",
                           runs = data.table::data.table(
                             Run = c("SRR001", "SRR002"),
                             LibraryLayout = c("SINGLE", "SINGLE"),
                             LIBRARYTYPE = c("RFP", "RFP"),
                             ScientificName = organism
                           ),
                           exp_name = paste0(accession, "-", gsub(" ", "_", tolower(organism)))) {
  stats::setNames(list(list(
    accession = accession,
    organisms = stats::setNames(list(list(conf = c(exp = exp_name))), organism),
    study = runs
  )), accession)
}

#' Mark every step of `steps` done for `experiment` in `config` (helper
#' for tests that need a "this experiment is already fully done" state).
fake_mark_all_done <- function(config, steps, experiment) {
  for (s in steps) set_flag(config, s, experiment)
}
