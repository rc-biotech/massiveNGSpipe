# All metadata
pipeline_init_all <- function(config, complete_metadata = 
                                file.path(config$project, "RFP_FINAL_LIST.csv")) {
  if (!file.exists(complete_metadata)) stop("You have not create a successful metadata table yet!")
  final_list <- fread(complete_metadata)
  accessions <- unique(final_list$study_accession)
  pipelines <- lapply(accessions, function(accession)
    pipeline_init(final_list[final_list$study_accession == accession,], accession,  config))
  return(pipelines)
}


# Pipeline 1: Download
pipeline_fetch <- function(pipelines, config) {
  for (pipeline in pipelines) {
    try <- try(
      pipeline_download(pipeline, config)  
    )
    if (is(try, "try-error")) 
      warning("Failed at step, fetch, study: ", pipeline$accession)
  }
}

# Pipeline 1: trim -> bam files
pipeline_trim_align <- function(pipelines, config) {
  for (pipeline in pipelines) {
    try <- try({
      pipeline_trim(pipeline,     config)
      pipeline_collapse(pipeline, config)
      pipeline_align(pipeline,    config)
      pipeline_cleanup(pipeline,  config)
    })
    if (is(try, "try-error")) 
      warning("Failed at step, trim_align, study: ", pipeline$accession)
  }
}

# Pipeline 3: OFST -> validated pshifts
pipeline_ofst_pshift <- function(pipelines, config) {
  for (pipeline in pipelines) {
    try <- try({
      df_list <- pipeline_create_experiment(pipeline, config)
      if (is.null(df_list)) next
      pipeline_create_ofst(df_list,                   config)
      pipeline_pshift(df_list,                        config)
      pipeline_validate_shifts(df_list,               config)
    })
    if (is(try, "try-error")) 
      warning("Failed at step, ofst_psfhit, study: ", pipeline$accession)
  }
}
#TODO: Add possibility to fix wrong pshifting

# Pipeline 4: Merge libraries: per study -> all merged per organism
pipeline_merge <- function(pipelines, config) {
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "merged_lib", e)))
  
  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e) read.experiment(e, validate = F))
      pipeline_merge_per_study(df_list, config)
    })
    if (is(try, "try-error")) 
      warning("Failed at step, ofst_psfhit, study: ", pipeline$accession)
  }
}