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
  cat("Last update: ")
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

report_failed_pipe <- function(try, config, step, exp) {

  if (is(try, "try-error")) {
    warning("Failed at step, ", step, "for exp: ", exp)
    warning(as.character(try))
    if (!is.null(config$error_dir)) {
      error_call <- paste(exp, "(", step, ")", as.character(try), paste(Sys.time()))
      error_path <- report_failed_pipe_path(config, exp)
      dir.create(config$error_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(error_call, error_path)
    }
    return(FALSE)
  }
  return(TRUE)
}

report_failed_pipe_path <- function(config, exp) {
  if (is.null(config$error_dir)) {
    warning("error_dir not defined in config, returning nonsense directory!")
    return(tempfile())
  }
  file.path(config$error_dir, paste0(exp, ".rds"))
}

session_error_dirs <- function(config) {
  sort(list.dirs(file.path(config$project, "error_logs"), recursive = F), decreasing = TRUE)
}

session_error_dirs_count <- function(config) {
  res <- sapply(session_error_dirs(config), function(x) length(list.files(x)))
  names(res) <- basename(names(res))
  return(res)
}
#' Get last active session error
#' @param config a list from a massiveNGSpipe::pipeline_config call
#' @return a character vector of errors
#' @export
last_session_errors <- function(config, index = 1, regex = NULL) {
  all_sessions <- session_error_dirs(config)
  if (length(all_sessions) < index) stop("You selected session ", index, " , but there is only ",
                                         length(all_sessions), " existing sessions")
  errors <- sapply(list.files(all_sessions[index], full.names = TRUE), readRDS)
  names(errors) <- gsub(".*error_logs|\\.rds$", "", names(errors))
  if (!is.null(regex)) errors <- grep(regex, errors, value = TRUE)
  return(errors)
}
