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
  do_trim <- "trim" %in% names(config$flag)
  do_collapse <- "collapsed" %in% names(config$flag)
  for (pipeline in pipelines) {
    try <- try({
      if (do_trim)
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
  do_contamint_removal <- "contam" %in% names(config$flag)
  for (pipeline in pipelines) {
    try <- try({
      if (do_contamint_removal)
        pipeline_align_contaminants(pipeline, config)

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
#' @param path_suffix character, default "". Add suffix to saved collection name, like a date to seperate version.
#' @return invisible(NULL)
#' @export
pipeline_collection_org <- function(config, pipelines = pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE),
                                    path_suffix = "") {
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

    libtype_df <- libraryTypes(df)
    if (length(libtype_df) != 1) stop("Only single libtype experiments supported for merging")
    if (libtype_df == "") stop("Libtype of experiment must be defined!")
    exp_name <- organism_collection_exp_name(org)
    if (libtype_df != "RFP") exp_name <- paste0(exp_name, "_", libtype_df)
    exp_name <- paste0(exp_name, path_suffix)
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
    count_folder <- QCfolder(df)
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
#' @param done_experiments named character vector, default:
#' step_is_done_pipelines(config, "merged_lib", pipelines). Set of experiments to merge
#' per species.
#' @param done_organisms which organisms to run merging for, default
#' unique(names(done_experiments))
#' @return invisible(NULL)
#' @export
pipeline_merge_org <- function(config, pipelines = pipeline_init_all(config, only_complete_genomes = TRUE,
                                                                     gene_symbols = FALSE),
                               done_experiments = step_is_done_pipelines(config, "merged_lib", pipelines),
                               done_organisms = unique(names(done_experiments)),
                               libtype_to_merge = libtype_long_to_short(config$preset),
                               out_dir_base = config$config["bam"]) {
  message("- Merge all samples per organism")
  done_organisms <- unique(done_organisms)
  stopifnot(length(names(done_experiments)) > 0)
  for (org in done_organisms) {
    message("-- Organism: ", org)
    df_list <- lapply(done_experiments[names(done_experiments) == org], function(e)
      read.experiment(e, validate = F)[1,])
    stopifnot(length(df_list) > 0)
    df <- do.call(rbind, df_list)
    df <- df[df$libtype == libtype_to_merge]

    # Overwrite default paths to merged
    libtype_df <- libraryTypes(df)

    df <- update_path_per_sample(df, libtype, libtype_to_merge = libtype_to_merge)

    exp_name <- organism_merged_exp_name(org)
    if (libtype_df != "RFP") exp_name <- paste0(exp_name, "_", libtype_df)
    out_dir <- file.path(out_dir_base, exp_name)
    ORFik::mergeLibs(df, out_dir, "all", "default", FALSE)

    make_additional_formats(df[1,], exp_name, out_dir)
  }
  return(invisible(NULL))
}

libtype_long_to_short <- function(long) {
  short <- long
  short[long == "Ribo-seq"] <- "RFP"
  short[long == "RNA-seq"] <- "RNA"
  return(short)
}

#' Make additional formats for an experiment
#'
#' Creates a new experiment and new file formats.
#' This is used for custom merge experiments mostly like all samples
#' per organism.
#' @param df_ref an ORFik experiment with only 1 libtype, giving info for new.
#' @param exp_name name of experiment to be made
#' @param out_dir_exp folder path, where study files are found.
#' @return invisible(NULL)
#' @export
make_additional_formats <- function(df_ref, exp_name, out_dir_exp) {
  message(exp_name)
  stopifnot(is(df_ref, "experiment"))

  libtype_df_ref <- libraryTypes(df_ref)
  stopifnot(length(libtype_df_ref) == 1)
  create.experiment(out_dir_exp, exper = exp_name,
                    txdb = df_ref@txdb,
                    libtype = libtype_df_ref,  fa = df_ref@fafile,
                    organism = organism(df_ref))
  df <- read.experiment(exp_name, output.env = new.env())

  convert_to_covRleList(df)
  if (libtype_df_ref == "RFP") {
    message("- Bigwig method: ORFik (5' ends)")
    convert_to_bigWig(df, in_files = filepath(df, "cov"))
  } else {
    message("- Bigwig method: full read")
    dir <- file.path(libFolder(df), "bigwig")
    dir.create(dir, showWarnings = FALSE)
    bw_files <- file.path(dir, c("all_forward.bigWig", "all_reverse.bigWig"))
    covrle <- fimport(filepath(df, "cov"))
    message("-- Bigwig forward")
    rtracklayer::export.bw(object = f(covrle), bw_files[1])
    message("-- Bigwig reverse")
    rtracklayer::export.bw(object = r(covrle), bw_files[2])
  }

  ORFik::countTable_regions(df, lib.type = "pshifted", forceRemake = TRUE)
  return(invisible(NULL))
}

#' Gather all modalities per organism into exp
#'
#' Using libtype modalities like Ribo-seq, RNA-seq, disome-seq, create
#' a new experiment. Only applies to organisms with > 1 modality completed.
#' @param all_exp data.table of all done experiment,
#'  default: list.experiments(validate = FALSE)
#' @return invisible(NULL)
#' @export
pipeline_merge_org_modalities <- function(all_exp = list.experiments(validate = FALSE)) {
  message("- All modalities per organism")
  candidates <- organism_merged_exp_name(unique(all_exp$organism))

  all_done <- all_exp[grep(paste(candidates, collapse = "|"), all_exp$name),]
  all_done <- all_done[!grepl("_modalities$", name),]
  modalities <- table(all_done$organism)
  done_organisms <- names(modalities[modalities > 1])
  for (org in done_organisms) {
    message("-- Organism: ", org)
    all_modalities_org <- all_done[organism == org,]
    df_list <- lapply(all_modalities_org$name, function(e)
      read.experiment(e, validate = F)[1,])
    df <- do.call(rbind, df_list)
    # Overwrite default paths to merged
    libtype_df <- unique(df$libtype)
    if (length(libtype_df) == 1) stop("Minimum 2 libtypes required for merge")
    if (any(libtype_df == "")) stop("Libtype of experiment must be defined!")
    exp_name <- organism_merged_modalities_exp_name(org)
    create.experiment(NULL, exper = exp_name,
                      txdb = df@txdb,
                      libtype = libtype_df, fa = df@fafile,
                      organism = org, files = df$filepath)
    df <- read.experiment(exp_name, output.env = new.env())
  }
  return(invisible(NULL))
}

pipeline_merge_exp_modalities <- function(all_exp = list.experiments(validate = FALSE)) {
  message("- All modalities per organism")
  exp_name <- sub("-.*", "", all_exp$name)
  candidates <- !grepl("^all_|_modalities$", all_exp$name)

  all_done <- all_exp[candidates,]
  exp_name <- exp_name[candidates]
  all_done[, exp_short := exp_name]
  libtypes_per_exp <- all_done[, table(exp_short, organism)]
  libtypes_per_exp <- rowSums(libtypes_per_exp > 1)
  libtypes_per_exp <- libtypes_per_exp[libtypes_per_exp > 0]

  multi_modals <- names(libtypes_per_exp)
  all_done <- all_done[exp_short %in% multi_modals]
  for (exp in multi_modals) {
    message("-- Experiment: ", exp)
    all_modalities <- all_done[exp_short == exp,]
    df_list <- lapply(all_modalities$name, function(e)
      read.experiment(e, validate = F))
    df <- do.call(rbind, df_list)
    # Overwrite default paths to merged
    libtype_df <- unique(df$libtype)
    if (length(libtype_df) == 1) stop("Minimum 2 libtypes required for merge")
    if (any(libtype_df == "")) stop("Libtype of experiment must be defined!")
    author <- ifelse(is.null(df$author), "", df$author)
    org <- organism(df)
    exp_name <- paste0(exp, "_", gsub(" ", "_", tolower(org)), "_modalities")
    create.experiment("", exper = exp_name, author = author,
                      txdb = df@txdb, libtype = df$libtype, rep = df$rep,
                      condition = df$condition, stage = df$stage,
                      fraction = df$fraction, runIDs = runIDs(df),
                      fa = df@fafile, organism = org, files = df$filepath)
    df <- read.experiment(exp_name, output.env = new.env())
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
