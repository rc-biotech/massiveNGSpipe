#' Initial path config
#'
#' Set up all paths and flag directories
#' @param project_dir where will specific pipeline outputs be put
#' @param config path, default \code{ORFik::config()}, where will
#' fastq, bam, references and ORFik experiments go
#' @param complete_metadata path, default: file.path(project_dir, "RFP_FINAL_LIST.csv")
#' Where should completed valid metadata be stored as csv?
#' @param google_url url or sheet object for google sheet to use. Set to NULL to
#' not use google sheet.
#' @param flags named character vector, with cut points, where can
#' the pipeline continue if it breaks? Is defined in combination with
#' 'steps' argument that defines that actuall function called for each break point.
#' Increasing break points makes the pipeline run faster at the cost of more
#' resource usage.
#' @param flag_steps a list of subsets of flags that maps which flags are set in
#' each function in 'step'.
#' @param mode = c("online", "local")[1]. "online" will assume project IDs for
#' online repository (SRA, ENA, PRJ etc). Local means local folders as accessions.
#' @param delete_raw_files logical, default: mode == "online". If online do delete
#' raw fastq files after trim step is done, for local samples do not delete.
#' Set only to TRUE for mode local, if you have backups!
#' @param steps a list of functions, the functions called during 'run_pipeline'
#' @param parallel_conf a bpoptions object, default:
#' \code{bpoptions(log =TRUE, stop.on.error = TRUE)}
#' Specific pipeline config for parallel settings and log directory
#' @param logdir = file.path(project_dir, "log_pipeline"),
#' @param BPPARAM = bpparam()
#' @return a list with a defined config
#' @export
pipeline_config <- function(project_dir = file.path(dirname(config)[1], "NGS_pipeline"),
                            config = ORFik::config(),
                            complete_metadata = file.path(project_dir, "RFP_FINAL_LIST.csv"),
                            backup_metadata = file.path(project_dir, "ALL_BACKUP_LIST.csv"),
                            temp_metadata = file.path(project_dir, "RFP_next_round_manual.csv"),
                            google_url = default_sheets(),
                            preset = "Ribo-seq",
                            flags = pipeline_flags(project_dir, mode),
                            flag_steps = config_substeps(names(flags), mode, preset),
                            pipeline_steps = preset_pipelines(preset, mode),
                            mode = c("online", "local")[1],
                            delete_raw_files = mode == "online",
                            parallel_conf = bpoptions(log =TRUE,
                                                      stop.on.error = TRUE),
                            logdir = file.path(project_dir, "log_pipeline"),
                            BPPARAM = bpparam()) {
  metadata <- paste0(project_dir, "/metadata")
  stopifnot(mode %in% c("online", "local"))
  stopifnot(is.logical(delete_raw_files))
  stopifnot(is(pipeline_steps, "list"))
  stopifnot(is(pipeline_steps[[1]], "function"))
  stopifnot(length(pipeline_steps) == length(flag_steps))
  if (mode == "local" & delete_raw_files)
    message("You have run local fastq files with delete raw files,",
            "are you sure this is what you want?")
  if (!is.null(parallel_conf)) {
    if (!is.null(logdir)) {
      dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
      bplog(BPPARAM) <- TRUE
      bplogdir(BPPARAM) <- logdir
    }
  }


  return(list(project = project_dir, config = config, flag = flags,
              flag_steps = flag_steps,
              pipeline_steps = pipeline_steps, metadata = metadata,
              complete_metadata = complete_metadata, backup_metadata = backup_metadata,
              temp_metadata = temp_metadata,
              google_url = google_url, mode = mode,
              delete_raw_files = delete_raw_files,
              parallel_conf = parallel_conf, BPPARAM = BPPARAM))
}

#' Configure which flags are done in each function
#' @inheritParams pipeline_config
#' @param flags character vector, names of flags
#' @param substeps a list of character vectors, where each sub list
#' is the steps that happen in that sub component of the pipeline
#' @return a list, returns the validated 'substeps' argument.
config_substeps <- function(flags, mode, preset = "Ribo-seq",
                            substeps = preset_substeps(preset, flags, mode)) {
  stopifnot(all(unlist(substeps) %in% flags))
  stopifnot(length(unlist(substeps)) > 0)
  stopifnot(length(flags) > 0)
  if (!all(flags %in% unlist(substeps)))
    message("Running pipeline with not all flags set, is this intentional?")
  return(substeps)
}

preset_substeps <- function(type = "Ribo-seq", flags, mode = c("online", "local")[1]) {
  valid_presets <- c("Ribo-seq")
  if (type == "Ribo-seq") {
    if (mode == "online") {
      list(pipe_fetch = flags[1:2],
           pipe_trim_col = flags[3:4],
           pipe_align = flags[5:6],
           pipe_exp = flags[7:8],
           pipe_pshift = flags[9:10],
           pipe_merge = flags[11],
           pipe_convert = flags[12:13],
           pipe_counts = flags[14])
    } else {
      list(pipe_trim_col = flags[1:2],
           pipe_align = flags[3:4],
           pipe_exp = flags[5:6],
           pipe_pshift = flags[7:8],
           pipe_merge = flags[9],
           pipe_convert = flags[10:11],
           pipe_counts = flags[12])
    }
  } else stop(paste("Currently valid preset pipelines are of types:",
                    valid_presets))
}


preset_pipelines <- function(type = "Ribo-seq", mode = c("online", "local")[1]) {
  valid_presets <- c("Ribo-seq")
  if (type == "Ribo-seq") {
    preset <- list(pipeline_fetch, pipeline_trim_collapse,
                   pipeline_align_clean, pipeline_exp_ofst,
                   pipeline_pshift_and_validate,
                   pipeline_merge_per_study,
                   pipeline_convert_psite_reads,
                   pipeline_counts_psites)
    if (mode == "local") preset <- preset[-1]
  } else stop(paste("Currently valid preset pipelines are of types:",
                    valid_presets))
  return(preset)
}
