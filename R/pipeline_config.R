#' Initial config setup
#'
#' Set up all paths, functions to be run and flag directories.
#' Also adds parallel processing settings and google integration.
#' @inheritParams libtype_flags
#' @param project_dir where will specific pipeline outputs be put. Default:
#'  file.path(dirname(config)[1], "NGS_pipeline").
#'  If you need seperat pipelines, please use different locations.
#' @param config path, default \code{ORFik::config()}, where will
#' fastq, bam, references and ORFik experiments go
#' @param complete_metadata path, default: file.path(project_dir, "FINAL_LIST.csv")
#' Where should completed valid metadata be stored as csv?
#' @param backup_metadata path, default: file.path(project_dir, "BACKUP_LIST.csv").
#' The complete list of all unique runs checked in this pipelines,
#'  even ones deleted earlier. Useful for check of what has been done before.
#' @param temp_metadata  path, default: file.path(project_dir, "next_round_manual.csv").
#' The intermediate file used after curation, but before it is validated. This is
#' the where you update the current new metadata annotations for final approval into the
#' complete_metadata. This syncs automatically to the google_url sheet if included.
#' @param blacklist path, Which studies to ignore for this config,
#' default: file.path(project_dir, "BLACKLIST.csv")\cr
#' A csv with 1 column id, which gives BioProject IDs
#' @param google_url url or sheet object for google sheet to use. Set to NULL to
#' not use google sheet.
#' @param flags named character vector, with cut points, where can
#' the pipeline continue if it breaks? Is defined in combination with
#' 'steps' argument that defines that actuall function called for each break point.
#' Increasing break points makes the pipeline run faster at the cost of more
#' resource usage.
#' @param flag_steps list, mapping of functions to ids. The names per list element is a function.
#' The character elements inside each list element are the flag name ids for all steps that will
#' be "marked as done" inside that function. Example: In default preset,
#' the pipe_trim_collapse function marks as done both the trim and collapse flags.
#' @param pipeline_steps a list of the functions to actually run, the functions
#' must be named equal to names of the 'flag_step' argument.
#' @param mode = \code{c("online", "local")[1]}. "online" will assume project IDs for
#' online repository (SRA, ENA, PRJ etc). Local means local folders as accessions.
#' @param delete_raw_files logical, default: mode == "online". If online do delete
#' raw fastq files after trim step is done, for local samples do not delete.
#' Set only to TRUE for mode local, if you have backups!
#' @param delete_trimmed_files logical, default: mode == "online".
#' If TRUE deletes the trimmed fasta files.
#' @param delete_collapsed_files logical, default: mode == "online".
#' If TRUE deletes the collapsed fasta files.
#' @param keep_contaminants logical, default FALSE. Do not keep contaminant aligned reads,
#'  if TRUE they are saved in contamination dir. Ignored if "contam" is not in flags to use.
#' @param keep_unaligned_genome logical, default FALSE. Do not keep contaminant aligned reads,
#'  else saved in contamination dir.
#' @param compress_raw_data logical, default FALSE. If TRUE, will compress raw fastq files.
#' @param parallel_conf a bpoptions object, default:
#' \code{bpoptions(log =TRUE, stop.on.error = TRUE)}
#' Specific pipeline config for parallel settings and log directory
#' @param verbose logical, default TRUE, give start up message
#' @param logdir = file.path(project_dir, "log_pipeline"),
#' @param BPPARAM_TRIM BiocParallel::MulticoreParam(8),
#' number of cores/threads to use for trimming. Optimal is 8 for most data.
#' @param BPPARAM = bpparam(), number of cores/threads to use.
#' @return a list with a defined config
#' @export
pipeline_config <- function(project_dir = file.path(dirname(config)[1], "NGS_pipeline"),
                            config = ORFik::config(),
                            complete_metadata = file.path(project_dir, "FINAL_LIST.csv"),
                            backup_metadata = file.path(project_dir, "BACKUP_LIST.csv"),
                            temp_metadata = file.path(project_dir, "next_round_manual.csv"),
                            blacklist = file.path(project_dir, "BLACKLIST.csv"),
                            google_url = default_sheets(project_dir),
                            preset = "Ribo-seq",
                            flags = pipeline_flags(project_dir, mode, preset, contam),
                            flag_steps = flag_grouping(flags),
                            pipeline_steps = lapply(names(flag_steps),
                                                    function(x) get(x, mode = "function")),
                            mode = c("online", "local")[1],
                            contam = FALSE,
                            delete_raw_files = mode == "online",
                            delete_trimmed_files = mode == "online",
                            delete_collapsed_files = mode == "online",
                            keep_contaminants = FALSE,
                            keep_unaligned_genome = FALSE,
                            compress_raw_data = FALSE,
                            parallel_conf = bpoptions(log =TRUE,
                                                      stop.on.error = TRUE),
                            logdir = file.path(project_dir, "log_pipeline"),
                            verbose = TRUE,
                            BPPARAM_TRIM = BiocParallel::MulticoreParam(8),
                            BPPARAM = bpparam()) {
  if (verbose) {
    message("Setting up mNGSp config..")
    message("- Preset (mode): ", preset, " (", mode, ")")
    message("- Workers (threads): ", BPPARAM$workers)
    message("- Sync to google: ", !is.null(google_url))
    message("- Metadata dir: ", project_dir)
    message("- Output data dir: ", config["bam"])
  }

  stopifnot(mode %in% c("online", "local"))
  stopifnot(is.logical(delete_raw_files) & is.logical(delete_trimmed_files) &
            is.logical(delete_collapsed_files) & is.logical(keep_contaminants) &
            is.logical(keep_unaligned_genome) & is.logical(compress_raw_data))
  stopifnot(is(pipeline_steps, "list"))
  if (preset == "empty") message("Using empty preset, now add your steps: using add_step_to_pipeline()")
  if (preset != "empty") stopifnot(is(pipeline_steps[[1]], "function"))
  stopifnot(length(pipeline_steps) == length(flag_steps))

  delete_any_files <- c(delete_raw_files, delete_trimmed_files, delete_collapsed_files)
  names(delete_any_files) <- c("raw", "trimmed", "collapsed")
  if (mode == "local" & any(delete_any_files))
    message("You have run local fastq files with delete ",
    paste(names(delete_any_files)[delete_any_files], collapse = " & "), " files,",
            "are you sure this is what you want?")
  if (!is.null(parallel_conf) && !is.null(logdir)) {
    BPPARAM <- set_log_dir(BPPARAM, logdir)
  }


  return(list(project = project_dir, config = config, flag = flags,
              flag_steps = flag_steps,
              pipeline_steps = pipeline_steps,
              metadata = paste0(project_dir, "/metadata"),
              complete_metadata = complete_metadata,
              backup_metadata = backup_metadata,
              temp_metadata = temp_metadata,
              blacklist = blacklist,
              google_url = google_url, mode = mode,
              delete_raw_files = delete_raw_files,
              delete_trimmed_files = delete_trimmed_files,
              delete_collapsed_files = delete_collapsed_files,
              keep_contaminants = keep_contaminants,
              keep.unaligned.genome = keep_unaligned_genome,
              compress_raw_data = compress_raw_data,
              preset = preset, parallel_conf = parallel_conf,
              BPPARAM_TRIM = BPPARAM_TRIM, BPPARAM = BPPARAM))
}

#' Define which flags are done in each function
#' @param flags character vector, of paths with names as flag steps and
#' attribute called "grouping" with length same as flag with grouping information,
#' where group names are the group functions.
#' @param substeps a list of character vectors, where each sub list
#' is the steps that happen in that sub component of the pipeline
#' @return a list, returns the validated 'substeps' argument.
#' @export
flag_grouping <- function(flags, substeps = preset_grouping(flags)) {
  if (length(flags) == 0) return(list())
  stopifnot(all(unlist(substeps) %in% names(flags)))
  stopifnot(length(unlist(substeps)) > 0)
  stopifnot(length(flags) > 0)
  if (!all(names(flags) %in% unlist(substeps)))
    message("Running pipeline with not all flags set, is this intentional?")
  return(substeps)
}

preset_grouping <- function(flags) {
  grouping <- attr(flags, "grouping")
  if (is.null(names(flags))) stop("flags must have names!")
  if (is.null(grouping)) stop("Your flags does not have defined attr 'grouping',",
                              "it should contain the grouping of flags")
  if (length(flags) != length(grouping)) stop("The grouping attribute must be",
                                              "equal size to number of flags!")
  return(split(names(flags), grouping)[unique(grouping)])
}
