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
  exp <- unlist(lapply(pipelines, 
                       function(x) lapply(x$organisms, 
                                          function(o) o$conf["exp"])),
                recursive = FALSE)
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