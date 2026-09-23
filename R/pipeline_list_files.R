# Generalizes list_files_of_type() (one pipeline, plain-directory steps
# only) to the whole `pipelines` list, and to steps whose real output isn't
# a fixed directory lookup:
# - "collapsed" files live inside the "trim" directory's SINGLE/PAIRED
#   subfolders (see pipeline_collapse()), not their own directory.
# - Everything from "exp" onward has an ORFik experiment object once that
#   flag is set, and ORFik::filepath(df, type) is the already-correct,
#   already-used-throughout-massiveNGSpipe way to get per-run paths for it
#   (ofst/pshifted/cov/bigwig) -- reused here instead of re-deriving
#   directory conventions a second time.
# - merged_lib/pcounts/valid_pshift are study-level (or QC-marker) outputs,
#   not one file per run; their real locations were found by reading
#   pipeline_merge.R/pipeline_preset_steps_sub.R/shifting_helpers.R
#   directly, not guessed.

# Maps a massiveNGSpipe step name to how its files should be located.
# kind: "dir" (list_dirs_of_pipeline() key, pre-experiment steps),
#       "collapsed" (special-cased "trim" subfolders),
#       "exp_marker" (no real per-run file, just whether the experiment
#       record itself exists),
#       "exp" (ORFik::filepath(df, key) once "exp" is done),
#       "exp_dir" (a fixed subfolder under libFolder(df) once "exp" is done).
step_file_location <- function(step) {
  switch(step,
        start = , fetch = list(kind = "dir", key = "fastq"),
        trim = list(kind = "dir", key = "trim"),
        collapsed = list(kind = "collapsed"),
        contam = list(kind = "dir", key = "contam"),
        aligned = , cleanbam = list(kind = "dir", key = "aligned", pattern = "\\.bam$"),
        exp = list(kind = "exp_marker"),
        ofst = , cigar_collapse = list(kind = "exp", key = "ofst"),
        pshifted = list(kind = "exp", key = "pshifted"),
        covrle = list(kind = "exp", key = "cov"),
        bigwig = list(kind = "exp", key = "bigwig"),
        valid_pshift = list(kind = "exp_dir", subdir = "QC_STATS",
                            pattern = "^(good|warning)\\.rds$"),
        merged_lib = list(kind = "exp_dir", subdir = "pshifted_merged"),
        pcounts = list(kind = "exp_dir", subdir = "QC_STATS",
                      pattern = "^countTable_"),
        stop("Unknown step for list_files_pipelines(): '", step,
             "'. Valid steps: ", paste(names(step_file_location_all()), collapse = ", "))
  )
}

# All step names step_file_location() understands (for error messages / discovery).
step_file_location_all <- function() {
  steps <- c("start", "fetch", "trim", "collapsed", "contam", "aligned", "cleanbam",
            "exp", "ofst", "cigar_collapse", "pshifted", "covrle", "bigwig",
            "valid_pshift", "merged_lib", "pcounts")
  setNames(lapply(steps, step_file_location), steps)
}

#' List existing files for one massiveNGSpipe pipeline step, across all pipelines
#'
#' A single, standardized entry point for "which studies/experiments have
#' files for step X", instead of checking flags alone or re-deriving each
#' step's file location by hand. Read-only: never downloads, creates, or
#' deletes anything.
#'
#' @param pipelines the pipelines list (see [pipeline_init_all])
#' @param config the mNGSp config object (see [pipeline_config])
#' @param step character, one of "start", "fetch", "trim", "collapsed",
#' "contam", "aligned", "cleanbam", "exp", "ofst", "cigar_collapse",
#' "pshifted", "covrle", "bigwig", "valid_pshift", "merged_lib", "pcounts".
#' @param full_names logical, default TRUE. If FALSE, returns basenames only.
#' @return a named list, one element per experiment (`"<accession>-<organism>"`),
#' each a character vector of existing file paths for that step
#' (`character(0)` if none exist yet, e.g. the step hasn't run for that
#' experiment). For `step = "exp"`, the "file" is a stand-in: the experiment
#' name itself if its ORFik experiment record exists, since an ORFik
#' experiment isn't stored as one single file on disk.
#' @export
#' @examples
#' \dontrun{
#' config <- pipeline_config()
#' pipelines <- pipeline_init_all(config)
#' list_files_pipelines(pipelines, config, "collapsed")
#' # Which studies have any collapsed files at all:
#' has_collapsed <- vapply(list_files_pipelines(pipelines, config, "collapsed"),
#'                         function(x) length(x) > 0, logical(1))
#' }
list_files_pipelines <- function(pipelines, config, step, full_names = TRUE) {
  loc <- step_file_location(step)
  out <- list()
  for (pipeline in pipelines) {
    dirs_per_org <- list_dirs_of_pipeline(pipeline)
    for (organism in names(pipeline$organisms)) {
      conf <- pipeline$organisms[[organism]]$conf
      experiment <- unname(conf["exp"])
      dirs <- dirs_per_org[[organism]]

      files <- switch(loc$kind,
        "dir" = {
          d <- dirs[[loc$key]]
          if (!dir.exists(d)) character(0) else
            list.files(d, pattern = loc$pattern, full.names = full_names)
        },
        "collapsed" = {
          d <- file.path(dirs[["trim"]], c("SINGLE", "PAIRED"))
          d <- d[dir.exists(d)]
          if (length(d) == 0) character(0) else list.files(d, full.names = full_names)
        },
        "exp_marker" = {
          if (step_is_done(config, "exp", experiment)) experiment else character(0)
        },
        "exp" = {
          if (!step_is_done(config, "exp", experiment)) character(0) else {
            df <- try(ORFik::read.experiment(experiment, validate = FALSE,
                                             output.env = new.env()), silent = TRUE)
            if (is(df, "try-error")) character(0) else {
              # unlist(): some types (e.g. "bigwig") return a list of
              # forward/reverse-strand path pairs per run, not one path
              # per run like "ofst"/"pshifted"/"cov" do.
              p <- unlist(ORFik::filepath(df, loc$key), use.names = FALSE)
              p <- p[file.exists(p)]
              if (!full_names) basename(p) else p
            }
          }
        },
        "exp_dir" = {
          if (!step_is_done(config, "exp", experiment)) character(0) else {
            df <- try(ORFik::read.experiment(experiment, validate = FALSE,
                                             output.env = new.env()), silent = TRUE)
            if (is(df, "try-error")) character(0) else {
              d <- file.path(ORFik::libFolder(df), loc$subdir)
              if (!dir.exists(d)) character(0) else
                list.files(d, pattern = loc$pattern, full.names = full_names)
            }
          }
        }
      )
      out[[experiment]] <- files
    }
  }
  return(out)
}
