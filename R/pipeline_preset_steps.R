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

pipeline_trim_collapse <- function(pipelines, config) {
  for (pipeline in pipelines) {
    try <- try({
      pipeline_trim(pipeline,     config)
      pipeline_collapse(pipeline, config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, trim_collapse, study: ", pipeline$accession)
  }
}

# Pipeline: trim -> bam files
pipeline_align_clean <- function(pipelines, config) {
  for (pipeline in pipelines) {
    try <- try({
      pipeline_align(pipeline,    config)
      pipeline_cleanup(pipeline,  config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, align, study: ", pipeline$accession)
  }
}

# Pipeline: exp -> ofst
pipeline_exp_ofst <- function(pipelines, config) {
  for (pipeline in pipelines) {
    try <- try({
      df_list <- pipeline_create_experiment(pipeline, config)
      if (is.null(df_list)) next
      pipeline_create_ofst(df_list,                   config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, exp_ofst, study: ", pipeline$accession)
  }
}

# Pipeline: pshifts -> validated pshifts
pipeline_pshift_and_validate <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e)
    step_is_next_not_done(config, "pshifted", e) |
      step_is_next_not_done(config, "valid_pshift", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipeline_pshift(df_list, accepted.length = c(20, 21, 25:33), config)
      pipeline_validate_shifts(df_list,               config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, pshift_validate, study: ",
              experiments)
  }
}
#TODO: Add possibility to fix wrong pshifting

# Pipeline, Merge libraries: per study -> all merged per organism
pipeline_merge_per_study <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "merged_lib", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipeline_merge_study(df_list, config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, merge_exp, study: ", experiments[1])
  }
}


# Pipeline 4: Merge libraries: per study -> all merged per organism
pipeline_convert_psite_reads <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "covrle", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e) read.experiment(e, validate = F))
      pipeline_convert_covRLE(df_list, config)
      pipeline_convert_bigwig(df_list, config)
    })
    if (is(try, "try-error")) {
      warning("Failed at step, pshift_convert, study: ", experiments[1])
      print(try)
    }
  }
}

#' Merge all studies per organism
#' @inheritParams run_pipeline
#' @return invisible(NULL)
#' @export
pipeline_merge_org <- function(config, pipelines = pipeline_init_all(config)) {
  names(pipelines) <- NULL
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_done(config, "merged_lib", e)))

  # Merge all per organism
  done_exp_list <- exp[done_exp]
  done_organisms <- unique(names(done_exp_list))
  for (org in done_organisms) {
    message("-- Organism: ", org)
    df_list <- lapply(done_exp_list[names(done_exp_list) == org], function(e)
      read.experiment(e, validate = F)[1,])
    df <- do.call(rbind, df_list)
    # Overwrite default paths to merged
    df@listData$filepath <- file.path(dirname(filepath(df, "default")), "pshifted_merged", "RFP.ofst")
    df@listData$rep <- seq(nrow(df))
    exp_name <- organism_merged_exp_name(org)
    out_dir <- file.path(config$config["bam"], exp_name)
    ORFik::mergeLibs(df, out_dir, "all", "default", FALSE)
    create.experiment(out_dir, exper = exp_name,
                      txdb = df@txdb,
                      libtype = "RFP",  fa = df@fafile, organism = org)
    df <- read.experiment(exp_name, output.env = new.env())
    convert_to_covRleList(df)
    convert_to_bigWig(df)
    ORFik::countTable_regions(df, lib.type = "pshifted", forceRemake = TRUE)
  }
  return(invisible(NULL))
}

pipeline_counts_psites <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "pcounts", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e) read.experiment(e, validate = F))
      pipeline_count_table_psites(df_list, config)
    })
    if (is(try, "try-error")) {
      warning("Failed at step, count tables psites, study: ", experiments[1])
      print(try)
    }
  }
  return(invisible(NULL))
}