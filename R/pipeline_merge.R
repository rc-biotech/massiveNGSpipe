# Pipeline, Merge libraries: per study -> all merged per organism
pipe_merge_study <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "merged_lib", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipeline_merge_study_single(df_list, config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, merge_exp, study: ", experiments[1])
  }
}

pipeline_merge_study_single <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "merged_lib", name(df))) next
    ORFik::mergeLibs(df, file.path(libFolder(df), "pshifted_merged"),
                     "lib", "pshifted", FALSE)
    set_flag(config, "merged_lib", name(df))
  }
}
