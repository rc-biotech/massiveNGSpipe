# Pipeline 1: Download
pipe_fetch <- function(pipelines, config) {
  for (pipeline in pipelines) {
    try <- try(
      pipeline_download(pipeline, config)
    )
    if (is(try, "try-error"))
      warning("Failed at step, fetch, study: ", pipeline$accession)
  }
}

pipe_trim_collapse <- function(pipelines, config) {
  do_collapse <- "collapsed" %in% names(config$flag)
  for (pipeline in pipelines) {
    try <- try({
      pipeline_trim(pipeline,     config)
      if (do_collapse)
        pipeline_collapse(pipeline, config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, trim_collapse, study: ", pipeline$accession)
  }
}

# Pipeline: trim -> bam files
pipe_align_clean <- function(pipelines, config) {
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
pipe_exp_ofst <- function(pipelines, config) {
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
pipe_pshift_and_validate <- function(pipelines, config) {
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


# Pipeline 4: Merge libraries: per study -> all merged per organism
pipe_convert <- function(pipelines, config) {
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

#' Create the superset collection of all samples per organism
#' @inheritParams run_pipeline
#' @return invisible(NULL)
#' @export
pipeline_collection_org <- function(config, pipelines = pipelines_init_all(config)) {
  message("- Collection of all samples per organism")
  names(pipelines) <- NULL
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_done(config, "pcounts", e)))

  # Merge all per organism
  done_exp_list <- exp[done_exp]
  done_organisms <- unique(names(done_exp_list))

  for (org in done_organisms) {
    message("-- Organism: ", org)
    df_list <- lapply(done_exp_list[names(done_exp_list) == org], function(e)
      read.experiment(e, validate = F))
    df <- do.call(rbind, df_list)
    exp_name <- organism_collection_exp_name(org)
    out_dir <- file.path(config$config["bam"], exp_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fraction <- df$fraction
    #fraction <- gsub("auto_.*", "", df$fraction)
    study_size <- unlist(lapply(df_list, function(e) nrow(e)))
    studies <- rep(unlist(done_exp_list[names(done_exp_list) == org], use.names = F), study_size)
    studies <- gsub("-.*", "", studies)
    fraction <- paste(fraction, studies, sep = "_")
    fraction <- gsub("^_", "",fraction)
    fraction <- gsub("^NA_", "",fraction)
    create.experiment(out_dir, exper = exp_name,
                      txdb = df@txdb, fa = df@fafile, organism = org,
                      libtype = df$libtype,
                      condition = df$condition, rep = df$rep,
                      stage = df$stage, fraction = fraction,
                      files = df$filepath, result_folder = out_dir)

    df <- read.experiment(exp_name, output.env = new.env())
    rel_dir <- "QC_STATS"
    count_folder <- file.path(out_dir, rel_dir)
    dir.create(count_folder, recursive = TRUE, showWarnings = FALSE)
    for (region in c("mrna", "cds", "leaders", "trailers")) {
      message("--- ", region)
      count_lists <- bplapply(df_list,
                             function(e, region) suppressMessages(countTable(e, region,
                                                           type = "summarized")),
                             region = region)
      count_list <- do.call(BiocGenerics::cbind, count_lists)
      saveName <- file.path(count_folder, paste0("countTable_", region, ".rds"))
      saveRDS(count_list, file = saveName)
      saveName <- file.path(count_folder, paste0("totalCounts_", region, ".rds"))
      saveRDS(colSums(assay(count_list)), file = saveName)
    }
  }
}

#' Merge all studies per organism
#'
#' Collected
#' @inheritParams run_pipeline
#' @return invisible(NULL)
#' @export
pipeline_merge_org <- function(config, pipelines = pipeline_init_all(config)) {
  message("- Merge all samples per organism")
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

pipe_counts <- function(pipelines, config) {
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
