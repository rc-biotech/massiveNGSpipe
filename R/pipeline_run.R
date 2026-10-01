#' Run massive_NGS_pipe
#'
#' Leaves a final status note in \code{checklist.txt}'s title line no
#' matter how the run ends: \code{"(done after N.N hours)"},
#' \code{"(stopped gracefully after N.N hours)"} (a deliberate
#' \code{request_stop()}, R/pipeline_abort.R), or \code{"(aborted after
#' N.N hours)"} for anything else -- an uncaught R error propagating out
#' of \code{bplapply()} (e.g. a worker error re-raised in the master
#' under \code{stop.on.error = TRUE}), or an interrupt (Ctrl+C/Esc)
#' delivered to this R session. This is an \code{on.exit()} handler, so
#' it genuinely cannot catch everything: a hard kill of the whole R
#' process (\code{kill -9}, an OOM-killer, a container/machine crash)
#' gives no R code any chance to run at all -- in that case
#' checklist.txt just goes stale (its own timestamp stops advancing),
#' which is itself the signal, just not a labeled one.
#' @param pipelines list, output of pipeline_init_all
#' @param config list, output from pipeline_config(), the global config for your
#' NGS pipeline
#' @param wait numeric, default 100 (in seconds). How long should each
#' partial pipeline wait if it done before it check for new results
#' ready to continue with.
#' @return invisible(NULL)
#' @export
run_pipeline <- function(pipelines, config, wait = 100) {

  config <- run_pipeline_set_up_session(pipelines, config)

  completed <- FALSE
  on.exit({
    if (!completed) {
      try(pipeline_checklist(pipelines, config, print = FALSE,
                             run_status = paste0("aborted after ", format_elapsed_hours(config$init_time))),
         silent = TRUE)
    }
  }, add = TRUE)

  # Run pipeline
  BiocParallel::bplapply(seq_along(config$pipeline_steps),
                         function(i, config, pipelines, wait)
      parallel_wrap(config$pipeline_steps[[i]], pipelines, config,
                    config$flag_steps[[i]], wait,
                    stage_name = names(config$flag_steps)[i]),
    pipelines = pipelines, wait = wait, config = config,
    BPPARAM = config$BPPARAM_MAIN)

  completed <- TRUE
  all_done <- all(pipeline_checklist(pipelines, config, print = FALSE)$state == "done")
  final_label <- if (all_done) "done" else if (stop_requested(config)) "stopped gracefully" else "finished early"
  try(pipeline_checklist(pipelines, config, print = FALSE,
                         run_status = paste0(final_label, " after ", format_elapsed_hours(config$init_time))),
     silent = TRUE)

  return(run_pipeline_end_session(pipelines, config))
}

#' One stage-group's idle-poll loop
#'
#' Repeatedly calls `function_call(pipelines, config)` until every step in
#' `steps` is done for every experiment, or a graceful stop is requested
#' (`stop_requested()`), sleeping `wait` seconds and re-rendering the
#' checklist between rounds.
#' @param function_call function(pipelines, config), a `pipe_*()` wrapper
#' @param pipelines the pipelines list
#' @param config the mNGSp config object
#' @param steps character vector, this stage-group's step ids (one
#' element of `config$flag_steps`)
#' @param wait numeric, seconds to sleep between idle rounds
#' @param stage_name character, this stage-group's name (one name of
#' `config$flag_steps`), used to record its exit via `mark_stage_exited()`
#' @return invisible(NULL)
#' @noRd
parallel_wrap <- function(function_call, pipelines, config, steps, wait = 100,
                          stage_name = NULL) {
  steps_merged <- paste(steps, collapse = ", ", sep = ", ")
  message("Start step pipeline:\n", steps_merged)
  exps <- pipelines_names(pipelines)
  steps_done <- all_substeps_done_all(config, steps, exps)
  while(!all(steps_done) && !stop_requested(config)) {
    function_call(pipelines, config)
    steps_done <- all_substeps_done_all(config, steps, exps)
    if (!all(steps_done) && !stop_requested(config)) {
      Sys.sleep(wait)
      pipeline_checklist(pipelines, config)
    }
  }
  # Always mark this stage-group's loop as exited, whether it finished
  # naturally or honored a stop request -- see mark_stage_exited()'s own
  # doc for why "only mark on stop" would leave abort_session_when_next_study_done_per_step()
  # waiting forever on stages that simply finished before any stop was requested.
  mark_stage_exited(config, stage_name, stopped = !all(steps_done))
  if (!all(steps_done)) {
    message("Stopped (graceful shutdown requested) for step pipeline:\n", steps_merged)
  } else {
    message("Done for step pipeline:\n", steps_merged)
  }
}

#' Validate inputs and initialize a run_pipeline() session
#'
#' Sets `config$BPPARAM_MAIN`/`init_time`/`error_dir`/`session_dir`,
#' creates the session directory, and writes its initial info file.
#' @inheritParams run_pipeline
#' @return `config`, with the session fields above set
#' @noRd
run_pipeline_set_up_session <- function(pipelines, config) {
  stopifnot(length(pipelines) > 0 & is(pipelines, "list"))
  stopifnot(!anyNA(names(pipelines)) & all(lengths(pipelines) == 3))
  message("---- Starting pipline:")
  config$BPPARAM_MAIN <- bpparam_from_config(config, "main")
  message("Number of workers: ", config$threads$default)
  message("Number of studies to run: ", length(pipelines))
  message("Steps to run: ", paste(names(config$flag), collapse = ", "))
  init_time <- Sys.time()
  init_time_char <- as.character.Date(init_time)
  config$init_time <- init_time
  config$error_dir <- file.path(config$project, "error_logs", init_time_char)
  config$session_dir <- file.path(config$project, "session_logs", init_time_char)
  dir.create(config$session_dir, recursive = TRUE, showWarnings = FALSE)
  session_info_save(config, pipelines)
  return(config)
}

#' Finalize a run_pipeline() session
#'
#' Determines success from whether `config$error_dir` has any recorded
#' errors, prints a progress report on success, updates the session's
#' info file status to "Completed"/"Failed", and optionally sends a
#' Discord notification.
#' @inheritParams run_pipeline
#' @return logical, TRUE if the session completed with no recorded errors
#' @noRd
run_pipeline_end_session <- function(pipelines, config) {
  # Done
  no_errors <- ifelse(!is.null(config$error_dir) && dir.exists(config$error_dir),
                      length(list.files(config$error_dir)) == 0,
                      TRUE)
  if (no_errors) {
    title_message <- "Pipeline is done, without any errors"
    progress_report(pipelines, config, show_stats = TRUE)
    info <- session_info_read(config)
    info$status <- "Completed"
    info$end_time <- Sys.time()
    saveRDS(info, file.path(config$session_dir, "session_info.rds"))

  } else {
    title_message <- paste("Pipeline is done, but had errors, skipping report.",
    "See the directory: ", config$error_dir, "for more information!")
    info <- session_info_read(config)
    info$status <- "Failed"
    info$end_time <- Sys.time()
    saveRDS(info, file.path(config$session_dir, "session_info.rds"))
  }
  message(title_message)
  if (!is.null(config$discord_webhook)) {
    exps <- pipelines_names(pipelines[1])
    exps_main <- gsub("-.*", "", exps)
    title_message <- paste0(config$preset, " ", title_message)
    message <- paste(title_message, "\n",
                     exps_main, "is done:\n" ,
                     "- experiments are named: ", paste(exps, collapse = " and "),
                     collapse = "\n")

    discord_connection_default_cached()
    discordr::send_webhook_message(message)
  }

  if (!is.null(config$init_time)) {
    cat("Total run time: "); print(round(Sys.time() - config$init_time, 2))
  }
  return(no_errors)
}

#' Pipeline experiment names
#'
#' Get all pipeline experiment names from pipeline
#' Remember all species per pipeline object is combined, so this
#' function unlists the whole when recursive is TRUE
#' @param pipelines a list, the pipelines object
#' @param recursive logical, default TRUE, If false return as list
#' @return character vector of experiment names, recursive TRUE gives list.
#' @export
pipelines_names <- function(pipelines, recursive = TRUE) {
  unlist(lapply(pipelines,
                function(x) lapply(x$organisms,
                                   function(o) o$conf["exp"])),
         recursive = recursive, use.names = FALSE)
}
