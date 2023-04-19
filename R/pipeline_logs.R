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
#' @param a config object
#' @param flagdir character, default: dirname(config$flag[1]),
#' path to flag directory
#' @return character, a date string
last_update <- function(config, flagdir = dirname(config$flag[1])) {
  stopifnot(dir.exists(flagdir))
  info <- dir_info(flagdir)
  mod_time <- info$modification_time
  return(max(mod_time))
}

#' Check if anything happened since given time
#'
#' @param a config object
#' @param max_time_since_update numeric, in days, default 1.
#' @return logical, FALSE if nothing happend since max time
#' @export
is_progressing <- function(config, max_time_since_update = 1) {
  last_update <- last_update(config)
  now <- Sys.time()
  diff <- difftime(now, last_update, units = "days")
  ifelse(diff > max_time_since_update,
         FALSE, TRUE)
}
