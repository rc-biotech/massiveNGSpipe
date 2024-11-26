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
  stopifnot(length(pipelines) > 0 & is(pipelines, "list"))
  message("---- Starting pipline:")
  message("Number of workers: ", BPPARAM$workers)
  substeps <- config$flag_steps
  pipeline_steps <- config$pipeline_steps
  # Run pipeline
  BiocParallel::bplapply(seq_along(pipeline_steps),
                         function(i, substeps, pipelines, pipeline_steps, wait)
    parallel_wrap(pipeline_steps[[i]], pipelines, config, substeps[[i]], wait),
    substeps = substeps, pipelines = pipelines, wait = wait,
    pipeline_steps = pipeline_steps, BPOPTIONS = config$parallel_conf, BPPARAM = BPPARAM)
  # Done
  message("Pipeline is done, creating report")
  progress_report(pipelines, config, show_stats = TRUE)
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
      message("Stopped sleeping (", steps_merged,")")
      idle_round <- idle_round + 1; message("- ", idle_round)
    }
  }
  message("Done for step pipeline:\n", steps_merged)
}

pipelines_names <- function(pipelines, recursive = TRUE) {
  unlist(lapply(pipelines,
                function(x) lapply(x$organisms,
                                   function(o) o$conf["exp"])),
         recursive = recursive, use.names = FALSE)
}
