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
run_pipeline <- function(pipelines, config, wait = 100, BPPARAM = bpparam()) {
  message("---- Starting pipline:")
  message("Number of workers: ", BPPARAM$workers)
  # Prepare sub steps of pipeline with flags
  pipeline_steps <- list(pipeline_fetch, pipeline_trim_align,
                         pipeline_ofst_pshift, pipeline_merge)
  flags <- names(config$flag)
  substeps <- list(pipe_fetch = flags[1],
                   pipe_align = flags[2:5],
                   pipe_pshift = flags[6:9],
                   pipe_merge = flags[10])
  # Run pipeline
  BiocParallel::bplapply(seq_along(pipeline_steps),
                         function(i, substeps, pipelines, pipeline_steps)
    parallel_wrap(pipeline_steps[[i]], pipelines, config, substeps[[i]], wait = 100),
    BPPARAM = BPPARAM, substeps = substeps, pipelines = pipelines,
    pipeline_steps = pipeline_steps)
  # Done
  message("Pipeline is done, creating report")
  progress_report(pipelines, config)
}

parallel_wrap <- function(function_call, pipelines, config, steps, wait = 100) {
  steps_not_done <- TRUE
  exp <-
  idle_round <- 0
  while(any(steps_not_done)) {
    function_call(pipelines, config)
    steps_not_done <- all_substeps_done(config, steps, exp)

    if (any(steps_not_done)) {
      message("Sleep")
      Sys.sleep(wait)
      message("Stopped sleeping")
      idle_round <- idle_round + 1; message("- ", idle_round)
    }
  }
  message("Done for step pipeline:")
  cat(steps, sep = ","); cat("\n")
}

pipelines_names <- function(pipelines, recursive = TRUE) {
  unlist(lapply(pipelines,
                function(x) lapply(x$organisms,
                                   function(o) o$conf["exp"])),
         recursive = recursive, use.names = FALSE)
}
