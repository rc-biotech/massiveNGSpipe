#' Run massive_NGS_pipe
#' @param pipelines list, output of pipeline_init_all
#' @param config list, output from pipeline_config(), the global config for your
#' NGS pipeline
#' @param wait numeric, default 100 (in seconds). How long should each
#' partial pipeline wait if it done before it check for new results
#' ready to continue with.
#' @param BPPARAM a BiocParallel Param object, default is the user standard,
#' namely bpparam(). To quickly check how many threads you use, do:
#' \code{bpparam()$workers}. To adjust number of threads do for instance:
#' \code{MulticoreParam(3)}, which gives 3 threads.
#' @return invisible(NULL)
#' @export
run_pipeline <- function(pipelines, config, wait = 100, BPPARAM = config$BPPARAM) {

  config <- run_pipeline_set_up_session(pipelines, config, BPPARAM)

  # Run pipeline
  BiocParallel::bplapply(seq_along(config$pipeline_steps),
                         function(i, config, pipelines, wait)
      parallel_wrap(config$pipeline_steps[[i]], pipelines, config,
                    config$flag_steps[[i]], wait),
    pipelines = pipelines, wait = wait, config = config,
    BPOPTIONS = config$parallel_conf, BPPARAM = config$BPPARAM_MAIN)

  return(run_pipeline_end_session(pipelines, config))
}

parallel_wrap <- function(function_call, pipelines, config, steps, wait = 100) {
  steps_merged <- paste(steps, collapse = " ,", sep = " ,")
  message("Start step pipeline:\n", steps_merged)
  exps <- pipelines_names(pipelines)
  idle_round <- 0
  steps_done <- all_substeps_done_all(config, steps, exps)
  while(!all(steps_done)) {
    function_call(pipelines, config)
    steps_done <- all_substeps_done_all(config, steps, exps)
    if (!all(steps_done)) {
      message("Sleep (", steps_merged,")")
      Sys.sleep(wait)
      idle_round <- idle_round + 1;
      message("Stopped sleeping (", steps_merged,")", " - ", idle_round)
    }
  }
  message("Done for step pipeline:\n", steps_merged)
}

run_pipeline_set_up_session <- function(pipelines, config, BPPARAM = config$BPPARAM) {
  stopifnot(length(pipelines) > 0 & is(pipelines, "list"))
  stopifnot(!anyNA(names(pipelines)) & all(lengths(pipelines) == 3))
  message("---- Starting pipline:")
  message("Number of workers: ", BPPARAM$workers)
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

    massiveNGSpipe:::discord_connection_default_cached()
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
