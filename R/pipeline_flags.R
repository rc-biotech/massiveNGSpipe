all_substeps_done <- function(config, steps, exp) {
  # Check which steps are done
  steps_not_done <- c()
  for (step in steps) {
    step_not_done <- FALSE
    for (e in exp) {
      if (!step_is_done(config, step, e)) {
        step_not_done <- TRUE
        break
      }
    }
    steps_not_done <- c(steps_not_done, step_not_done)
  }
  return(steps_not_done)
}

pipeline_flags <- function(project_dir) {
  flag_dir <- file.path(project_dir, "flags")
  flag_sub <- c("fetch", "trim", "collapsed", "aligned", "cleanbam","exp","ofst", 
                "pshifted", "valid_pshift", 
                "merged_lib")
  flags <- file.path(flag_dir, flag_sub)
  names(flags) <- flag_sub
  lapply(flags, function(f) dir.create(f, showWarnings = FALSE, recursive = TRUE))
  return(flags)
}

step_is_done <- function(config, step, experiment) {
  step_dir <- config[["flag"]][step]
  flag_rds <- file.path(step_dir, paste0(experiment, ".rds"))
  return(file.exists(flag_rds))
}

step_is_next_not_done <- function(config, step, experiment) {
  if (step == "fetch") stop("fetch should use step_is_done directly!")
  all_steps <- names(config[["flag"]])
  all_flags <- file.exists(file.path(config[["flag"]], paste0(experiment, ".rds")))
  index <- which(all_steps %in% step)
  return(all(all_flags[seq(index-1)]) & !all_flags[index])
}

set_flag <- function(config, step, experiment, value = TRUE) {
  step_dir <- config[["flag"]][step]
  if (!dir.exists(step_dir)) stop(step, " is not a valid flag dir!")
  flag_rds <- file.path(step_dir, paste0(experiment, ".rds"))
  saveRDS(value, flag_rds)
  return(invisible(NULL))
}

#' Set flags for all pipelines given
#' 
#' Set all flags as done for all experiments in pipeline up to and included
#' the last step in steps
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
  return(invisible(NULL))
}

remove_flag <- function(config, step, experiment) {
  step_dir <- config[["flag"]][step]
  flag_rds <- file.path(step_dir, paste0(experiment, ".rds"))
  file.remove(flag_rds)
}