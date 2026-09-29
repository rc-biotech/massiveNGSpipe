# Graceful drain-and-abort for a running pipeline session.
#
# Signalling is file-based, same convention as checklist.txt
# (pipeline_checklist.R) -- BiocParallel's log=TRUE/logdir= buffering
# means a running worker can't be reached by any other in-process
# mechanism, so a plain file under the session directory, polled between
# studies, is what every other cross-process signal in this package
# already uses (checklist.txt, sample_flags/*, flags/*).

#' Path to this session's stop-request marker
#' @inheritParams pipeline_log_base
#' @noRd
stop_signal_path <- function(config) file.path(pipeline_log_base(config), "stop_requested.rds")

#' Has a graceful stop been requested for this session?
#' @inheritParams pipeline_log_base
#' @noRd
stop_requested <- function(config) file.exists(stop_signal_path(config))

#' Request a graceful stop for this session (internal -- see
#' abort_session_when_next_study_done_per_step() for the actual
#' user-facing entry point, which also waits for the drain and
#' terminates the mother process).
#' @inheritParams pipeline_log_base
#' @noRd
request_stop <- function(config) {
  d <- dirname(stop_signal_path(config))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  saveRDS(TRUE, stop_signal_path(config))
  invisible(NULL)
}

#' Mark one stage-group's parallel_wrap() loop as exited -- for EITHER
#' reason (genuine completion or an honored stop request), always, not
#' just the stopped case. This matters: a stage-group that finishes
#' naturally (the common case for most stages, especially when
#' BPPARAM_MAIN is SerialParam and stage-groups run one after another
#' rather than concurrently) exits long before any stop is ever
#' requested, and would otherwise never get a marker at all --
#' abort_session_when_next_study_done_per_step() would then wait forever
#' for a marker that stage was never going to produce. Recording the
#' reason in the stored value lets a caller distinguish the two if it
#' cares, but "has this stage-group's loop exited at all" (either reason)
#' is the only thing abort_session_when_next_study_done_per_step() itself
#' needs.
#' @inheritParams pipeline_log_base
#' @param stage_name character, the pipe_* group name (e.g. "pipe_trim_collapse")
#' @param stopped logical, TRUE if this exit was due to an honored stop
#' request, FALSE if the stage-group had simply finished all its work.
#' @noRd
mark_stage_exited <- function(config, stage_name, stopped) {
  if (is.null(stage_name)) return(invisible(NULL))
  d <- file.path(pipeline_log_base(config), "stage_stopped")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  saveRDS(stopped, file.path(d, paste0(stage_name, ".rds")))
  invisible(NULL)
}

#' Gracefully halt a running pipeline session
#'
#' Stops every pipeline stage from picking up any *new* study -- whichever
#' study each stage-group is currently mid-way through still finishes
#' normally -- and once every stage has drained, terminates the main
#' \code{\link{run_pipeline}} process for that session. Intended for
#' planned maintenance/restarts without corrupting in-flight work.
#'
#' Callable from a separate R session than the one actually running
#' \code{run_pipeline()} (same convention already used by
#' \code{mNGSp_app()}/\code{watch_pipeline_checklist()} to observe another
#' running session by its timestamp) -- only \code{config$project} and
#' \code{config$flag_steps} are used from \code{config} here, not any
#' live session state.
#'
#' @param config the mNGSp config object for this project (as returned by
#' \code{\link{pipeline_config}}) -- used for \code{config$project} and
#' \code{config$flag_steps} (the set of stage-group names to wait on),
#' not for any session-specific fields.
#' @param session character, the session timestamp naming
#' \code{file.path(config$project, "session_logs", session)} -- i.e.
#' \code{config$session_dir} as set once per \code{run_pipeline()} call
#' by \code{run_pipeline_set_up_session()}. Find it via
#' \code{session_info_table(config)}.
#' @param poll_interval numeric, seconds between drain-status checks,
#' default 10.
#' @param timeout numeric, seconds to wait for all stages to drain before
#' giving up (the stop request is left in place either way -- no new
#' study will start even if this times out), default 3600.
#' @return invisible(TRUE) if the mother process was signalled to
#' terminate, invisible(FALSE) if the session had already finished on its
#' own (or timed out) before that could happen.
#' @export
abort_session_when_next_study_done_per_step <- function(config, session,
                                                        poll_interval = 10,
                                                        timeout = 3600) {
  session_dir <- file.path(config$project, "session_logs", session)
  if (!dir.exists(session_dir)) stop("No such session directory: ", session_dir)
  session_config <- list(session_dir = session_dir, project = config$project)

  info <- session_info_read(session_config)
  if (info$status %in% c("Completed", "Failed")) {
    message("Session already finished (status: ", info$status, "), nothing to abort.")
    return(invisible(FALSE))
  }

  message("Requesting graceful stop for session: ", session)
  request_stop(session_config)

  stage_names <- names(config$flag_steps)
  stopped_dir <- file.path(session_dir, "stage_stopped")
  start_time <- Sys.time()
  repeat {
    info <- session_info_read(session_config)
    if (info$status %in% c("Completed", "Failed")) {
      message("Session finished on its own while waiting for it to drain.")
      return(invisible(FALSE))
    }
    stopped <- if (dir.exists(stopped_dir)) sub("\\.rds$", "", list.files(stopped_dir)) else character()
    if (all(stage_names %in% stopped)) break
    if (as.numeric(difftime(Sys.time(), start_time, units = "secs")) > timeout) {
      warning("Timed out waiting for all stages to drain; the stop request is ",
              "still in place (no new study will start), but the mother process ",
              "was NOT terminated. Stages still active: ",
              paste(setdiff(stage_names, stopped), collapse = ", "))
      return(invisible(FALSE))
    }
    Sys.sleep(poll_interval)
  }

  # All stage-groups have exited by this point, so the mother process's
  # bplapply() call has nothing left to wait on and may well have already
  # returned and exited on its own (run_pipeline_end_session() already
  # wrote "Completed") in the moment between the drain check above and
  # here -- that race is expected and not an error, so re-read status
  # once more and only send a signal / overwrite it if the process (and
  # its "started" status) are still actually there.
  info <- session_info_read(session_config)
  if (info$status %in% c("Completed", "Failed")) {
    message("Session finished on its own right as it finished draining -- nothing left to terminate.")
    return(invisible(FALSE))
  }

  message("All stages drained. Terminating session process group (pgid ", info$pgid, ")...")
  still_alive <- system(paste0("kill -0 -", info$pgid), ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  if (still_alive) system(paste0("kill -TERM -", info$pgid), ignore.stdout = TRUE, ignore.stderr = TRUE)
  info$status <- "Aborted"
  info$end_time <- Sys.time()
  saveRDS(info, file.path(session_dir, "session_info.rds"))
  invisible(TRUE)
}
