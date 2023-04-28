#' Init all studies to pipeline objects
#' @inheritParams curate_metadata
#' @param simple_progress_report logical, default TRUE. Display current progress.
#' @return a list of pipelines
#' @export
pipeline_init_all <- function(config, complete_metadata = config$complete_metadata,
                              simple_progress_report = TRUE) {
  if (!file.exists(complete_metadata)) stop("You have not create a successful metadata table yet!")
  final_list <- fread(complete_metadata)[KEEP == TRUE,]
  if (nrow(final_list) == 0) stop("complete metadata table has 0 rows, ",
                         "did you forget to set the 'KEEP' column to TRUE in",
                         " 'config$temp_metadata'?")
  accessions <- unique(final_list$study_accession)
  pipelines <- lapply(accessions, function(accession)
    tryCatch(
      pipeline_init(final_list[final_list$study_accession == accession,], accession,  config),
            error=function(e) {warning(e); return(NULL)}))
  names(pipelines) <- accessions

  crashed_studies <- pipelines %in% list(NULL)
  if (any(crashed_studies)) {
    pipelines <- pipelines[!crashed_studies]
  }


  if (simple_progress_report) progress_report(pipelines, config, FALSE)
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
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
  done_exp <- unlist(lapply(exp, function(e)
    step_is_next_not_done(config, "pshifted", e) |
      step_is_next_not_done(config, "valid_pshift", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipeline_pshift(df_list,                        config)
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
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
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
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
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
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
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
    df <- read.experiment(exp_name)
    convert_to_covRleList(df)
    convert_to_bigWig(df)
    ORFik::countTable_regions(df, lib.type = "pshifted", forceRemake = TRUE)
  }
  return(invisible(NULL))
}

pipeline_counts_psites <- function(pipelines, config) {
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
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

organism_merged_exp_name <- function(organisms) {
  paste0("all_merged-", gsub(" ", "_", organisms))
}

