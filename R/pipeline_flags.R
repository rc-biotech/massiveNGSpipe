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
#' @param flag_sub character vector of names of flag steps to include
#' @param flags_to_use numeric, vector of indices of 'flag_sub' to use.
#' c(1,2), would mean only download.
#' @return a named character vector of directories in 'project_dir'
#'  named as 'flag_sub' with names 'flag_sub'
pipeline_flags <- function(project_dir, mode = c("online", "local")[1],
                           flag_sub = c("start","fetch", "trim", "collapsed",
                                        "aligned", "cleanbam","exp","ofst",
                                        "pshifted", "valid_pshift",
                                        "merged_lib",
                                        "covrle", "bigwig",
                                        "pcounts"),
                           flags_to_use = if(mode == "online")
                             {seq(length(flag_sub))} else seq(length(flag_sub))[-c(1,2)]) {
  flag_dir <- file.path(project_dir, "flags")
  flag_sub <- flag_sub[flags_to_use]
  flags <- file.path(flag_dir, flag_sub)
  stopifnot(length(flags) > 0)
  names(flags) <- flag_sub
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
#' Set all flags as not done for all experiments in pipeline up to
#' and included the last step in steps
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
