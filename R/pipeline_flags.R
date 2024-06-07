#' Define the library specific steps
#'
#' This is called flags
#' @param preset character, default "Ribo-seq".
#'  Alternatives: c("Ribo-seq", "RNA-seq", "disome", "empty")
#' @param mode either of c("online", "local")
#' @return a character vector with names being the grouping in functions
libtype_flags <- function(preset, mode = "online") {
  valid_presets <- c("Ribo-seq", "RNA-seq", "disome", "empty")

  download_flags <- if (mode == "online") {
    c(fetch = "start", fetch = "fetch")
  } else {
    c()
  }
  trim_flags <- c(trim_collapse = "trim", trim_collapse = "collapsed")
  align_flags <- c(align_clean = "aligned", align_clean = "cleanbam")
  exp_flags <- c(exp_ofst = "exp", exp_ofst = "ofst")
  merge_flags <- c(merge_study = "merged_lib")

  lib_type_flags <- if (preset == "Ribo-seq") {
    c(pshift_and_validate = "pshifted", pshift_and_validate = "valid_pshift")
  } else if (preset %in% c("RNA-seq", "disome")) {
    trim_flags <- trim_flags[1]
    c(cigar_collapse = "cigar_collapse")
  } else stop(paste("Currently valid preset pipelines are of types:",
                    valid_presets))
  format_processing <- c(convert = "covrle", convert = "bigwig")
  count_tables <- c(counts = "pcounts")

  full_pipe <- c(download_flags, trim_flags, align_flags, exp_flags,
                 lib_type_flags, merge_flags, format_processing, count_tables)
  names(full_pipe) <- paste0("pipe_", names(full_pipe))
  return(full_pipe)
}

all_substeps_done <- function(config, steps, exp) {
  # Check which steps are done
  steps_done <- c()
  for (step in steps) {
    step_done <- TRUE
    for (e in exp) {
      if (!step_is_done(config, step, e)) {
        step_done <- FALSE
        break
      }
    }
    steps_done <- c(steps_done, step_done)
  }
  return(steps_done)
}

all_substeps_done_all <- function(config, steps, exps) {
  # Check which steps are done
  if (length(exps) == 0) stop("exps must have length > 0")
  steps_done <- c()
  for (exp in exps) {
    steps_done <- c(steps_done, all_substeps_done(config, steps, exp))
  }
  return(steps_done)
}


#' Get all pipeline flags and set dirs
#'
#' @param project_dir pipeline directory
#' @inheritParams pipeline_config
#' @param flag_names character vector of names of flag steps to include
#' @param create_dirs logical, default TRUE. For test purpose you can turn of
#' dir creation if needed.
#' @return a named character vector of directories in 'project_dir'
#'  named as 'flag_sub' with names 'flag_sub'
#' @export
pipeline_flags <- function(project_dir, mode = c("online", "local")[1],
                           preset,
                           flag_names = libtype_flags(preset, mode),
                           create_dirs = TRUE
                           ) {
  if (preset == "empty") return(character())
  flag_dir <- file.path(project_dir, "flags")
  flags <- file.path(flag_dir, flag_names)
  stopifnot(length(flags) > 0)
  names(flags) <- flag_names
  attr(flags, "grouping") <- names(flag_names)
  if (create_dirs)
    lapply(flags, function(f) dir.create(f, showWarnings = FALSE, recursive = TRUE))
  return(flags)
}

step_is_done <- function(config, step, experiment) {
  step_dir <- config[["flag"]][step]
  flag_rds <- file.path(step_dir, paste0(experiment, ".rds"))
  return(file.exists(flag_rds))
}

step_is_next_not_done <- function(config, step, experiment
                                  ) {
  if (step == "fetch") stop("fetch should use step_is_done directly!")
  all_steps <- names(config[["flag"]])
  all_flags <- file.exists(file.path(config[["flag"]], paste0(experiment, ".rds")))
  index <- which(all_steps %in% step)
  all_previous_done <- if (index == 1) {
    TRUE
  } else all(all_flags[seq(index-1)])
  this_one_not_done <- !all_flags[index]
  return(all_previous_done & this_one_not_done)
}

set_flag <- function(config, step, experiment, value = TRUE) {
  step_dir <- config[["flag"]][step]
  if (!dir.exists(step_dir)) stop(step, " is not a valid flag dir!")
  flag_rds <- file.path(step_dir, paste0(experiment, ".rds"))
  saveRDS(value, flag_rds)
  return(invisible(NULL))
}

#' Set flags for all pipelines specified
#'
#' Set all flags as done for all experiments in pipeline up to and included
#' the last step in steps
#' @inheritParams run_pipeline
#' @param steps, which flags to set to TRUE for given pipeline objects
#' @return invisible(NULL)
#' @export
set_flag_all <- function(config, steps = names(config$flag), pipelines) {
  if (length(steps) == 1) {
    if (steps == "all") {
      steps <- names(config$flag)
    }
  }
  if (!all(steps %in% names(config$flag))) stop("Some steps given are not valid steps!")
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
  for (e in exp) {
    for (step in steps) set_flag(config, step, e, value = TRUE)
  }
  message("Updated flags for pipeline subset:")
  progress_report(pipelines, config, show_stats = FALSE)
  message("Flag set: Done")
  return(invisible(NULL))
}

remove_flag <- function(config, step, experiment, warning = FALSE) {
  step_dir <- config[["flag"]][step]
  flag_rds <- file.path(step_dir, paste0(experiment, ".rds"))
  if (warning) {
    file.remove(flag_rds)
  } else suppressWarnings(file.remove(flag_rds))
}

#' Remove flags for all pipelines specified
#'
#' Set all flags specified as not done for all experiments in
#' for all pipeline objects specified.
#' @inheritParams run_pipeline
#' @param steps, which flags to set to FALSE (not done)
#'  for given pipeline objects
#' @return invisible(NULL)
#' @export
remove_flag_all <- function(config, steps = names(config$flag), pipelines) {
  if (length(steps) == 1) {
    if (steps == "all") {
      steps <- names(config$flag)
    }
  }
  if (!all(steps %in% names(config$flag))) stop("Some steps given are not valid steps!")
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
  for (e in exp) {
    for (step in steps) remove_flag(config, step, e)
  }
  message("Removed flags for pipeline subset:")
  progress_report(pipelines, config, show_stats = FALSE)
  message("Flag remove: Done")
  return(invisible(NULL))
}

#' Add additional step to pipeline
#'
#' @param config list, a NGS pipeline object
#' @param short_name character, the flag id, a short identifier (differential expression could be: difexp etc)
#' @param FUN a function, the to be added
#' @param group_name character, name of the function FUN.
#' @param tail_or_head character, c("tail", "head")[1]. Default is to append last, can also front append.
#' @return list, the updated config object
#' @export
add_step_to_pipeline <- function(config, short_name, FUN, group_name = name_of_function(FUN), tail_or_head = "tail") {
  stopifnot(is.list(config) & !is.null(config$project))
  stopifnot(is.character(short_name))
  stopifnot(is.function(FUN))
  stopifnot(is.character(group_name))
  stopifnot(tail_or_head %in% c("tail", "head"))

  flags <- config$flag
  grouping <- attr(flags, "grouping") # Grouping attribute

  if (short_name %in% names(flags)) stop("Step short_name already exists in pipeline config!")
  if (group_name %in% grouping) stop("Function name already exist in config!")

  if (tail_or_head == "tail") {
    flags <- c(flags, file.path(dirname(flags[1]), short_name))
    names(flags)[length(flags)] <- short_name
    grouping <- c(grouping, group_name) # <---- name of new function
  } else {
    flags <- c(file.path(dirname(flags[1]), short_name), flags)
    names(flags)[1] <- short_name
    grouping <- c(group_name, grouping) # <---- name of new function
  }

  dir.create(flags[short_name], showWarnings = FALSE, recursive = TRUE)

  attr(flags, "grouping") <- grouping
  flag_steps <- flag_grouping(flags)

  pipeline_steps <- lapply(names(flag_steps), function(x) get(x, mode = "function"))
  config$flag <- flags
  config$flag_steps <- flag_steps
  config$pipeline_steps <- pipeline_steps
  return(config)
}
