set_log_dir <- function(BPPARAM, logdir) {
  stopifnot(is.character(logdir))
  dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
  bplog(BPPARAM) <- TRUE
  bplogdir(BPPARAM) <- logdir
  return(BPPARAM)
}

#' Get log for specific pipeline run
#'
#' Get the n lasts line of a specific job for pipeline
#' @param config a config object
#' @param task numeric, default 1. The task number wanted to check.
#' @param lines numeric, default 100. Number of lines to output (
#' starting from end of fail)
#' @param logdir character, default: bplogdir(config$BPPARAM),
#' path to log directory you want to check.
#' @return NULL (text output to screen)
#' @export
progress_log <- function(config, task = 1, lines = 100,
                         logdir = bplogdir(config$BPPARAM)) {
  stopifnot(dir.exists(logdir))
  log_job <- paste0("BPJOB.task", task, ".log")
  log_file <- file.path(logdir, log_job)
  system(paste("tail -n", lines, log_file))
}

#' Get latest time of pipeline update flag
#'
#' @param config a pipeline_config object
#' @param flagdir character, default: \code{dirname(config$flag[1])},
#' path to flag directory
#' @return character, a date string
last_update <- function(config, flagdir = dirname(config$flag[1])) {
  stopifnot(dir.exists(flagdir))
  info <- dir_info(flagdir)
  mod_time <- info$modification_time
  return(max(mod_time))
}

#' Get last pipeline update, in days relative to now
#'
#' @inheritParams last_update
#' @param units character, default "days"
#' @return character, a date string
last_update_diff <- function(config, units = "days") {
  last_update <- last_update(config)
  now <- Sys.time()
  return(difftime(now, last_update, units = units))
}

#' Get last update time message
#' @param config a list, the pipeline config
#' @param units = "hours"
#' @return invisible(NULL), only cat prints the message
#' @export
last_update_message <- function(config, units = "hours") {
  cat("Last update:")
  cat("(", round(last_update_diff(config, units = "hours"), 1), " hours ago): ",
      format(last_update(config), usetz=TRUE), "\n", sep = "")
  return(invisible(NULL))
}

#' Check if anything happened since given time
#'
#' @inheritParams last_update
#' @param max_time_since_update numeric, in days, default 1.
#' 1 Hour is 0.041
#' @return logical, FALSE if nothing happend since max time
#' @export
is_progressing <- function(config, max_time_since_update = 1) {
  last_update_diff(config) < max_time_since_update
}
