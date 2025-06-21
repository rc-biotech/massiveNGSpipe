# Pipeline 1: Download
pipe_fetch <- function(pipelines, config) {
  drive_cap <- config$stop_downloading_new_data_at_drive_usage
  disc_is_full <- as.numeric(gsub("%", "", get_system_usage()$Drive_Usage_Percent)) >= drive_cap
  if (disc_is_full) {
    message("Disc is nearly full, will pause downloading data for a while..")
    return(invisible(NULL))
  }

  flag_step <- which(names(config$flag) == "fetch") # Step after fetch
  unprocessed_downloads <- sum(progress_report(pipelines, config,
                                          show_status_per_exp = FALSE,
                                          show_done = FALSE,
                                          return_progress_vector = T) == flag_step)
  if (unprocessed_downloads > config$max_unprocessed_downloads) {
    message("Max unprocessed downloads reached, will pause downloading data for a while..")
    return(invisible(NULL))
  }

  for (pipeline in pipelines) {
    exp <- pipeline$accession
    if (file.exists(report_failed_pipe_path(config, exp))) next
    try <- try(
      pipeline_download(pipeline, config)
    )
    report_failed_pipe(try, config, "fetch", exp)
  }
}

pipe_trim_collapse <- function(pipelines, config) {
  do_trim <- "trim" %in% names(config$flag)
  do_collapse <- "collapsed" %in% names(config$flag)
  for (pipeline in pipelines) {
    exp <- pipeline$accession
    if (file.exists(report_failed_pipe_path(config, exp))) next
    try <- try({
      if (do_trim)
        pipeline_trim(pipeline,     config)
    })
    status <- report_failed_pipe(try, config, "trim", pipeline$accession)
    if (status) {
      try <- try({
        if (do_collapse)
          pipeline_collapse(pipeline, config)
      })
      report_failed_pipe(try, config, "collapse", pipeline$accession)
    }
  }
}

# Pipeline: trim -> bam files
pipe_align_clean <- function(pipelines, config) {
  do_contamint_removal <- "contam" %in% names(config$flag)
  for (pipeline in pipelines) {
    exp <- pipeline$accession
    if (file.exists(report_failed_pipe_path(config, exp))) next
    try <- try({
      if (do_contamint_removal)
        pipeline_align_contaminants(pipeline, config)

      pipeline_align(pipeline,    config)

    })
    status <- report_failed_pipe(try, config, "align", pipeline$accession)
    if (status) {
      try <- try({
        pipeline_cleanup(pipeline,  config)
      })
      report_failed_pipe(try, config, "clean", pipeline$accession)
    }
  }
}

# Pipeline: exp -> ofst
pipe_exp_ofst <- function(pipelines, config) {
  for (pipeline in pipelines) {
    exp <- pipeline$accession
    if (file.exists(report_failed_pipe_path(config, exp))) next
    try <- try({
      df_list <- pipeline_create_experiment(pipeline, config)
    })
    status <- report_failed_pipe(try, config, "exp", pipeline$accession)
    if (is.null(df_list)) next
    if (status) {
      try <- try({
        pipeline_create_ofst(df_list,                   config)
      })
      report_failed_pipe(try, config, "ofst", pipeline$accession)
    }
  }
}

# Pipeline: pshifts -> validated pshifts
pipe_pshift_and_validate <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e)
    step_is_next_not_done(config, "pshifted", e) |
      step_is_next_not_done(config, "valid_pshift", e)))

  for (experiments in exp[done_exp]) {
    if (file.exists(report_failed_pipe_path(config, experiments[1]))) next

    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipeline_pshift(df_list,          config)
    })
    status <- report_failed_pipe(try, config, "pshift", experiments[1])
    if (status) {
      try <- try({
        pipeline_validate_shifts(df_list, config)
      })
      report_failed_pipe(try, config, "validate", experiments[1])
    }

  }
}
#TODO: Add possibility to fix wrong pshifting


# Pipeline 4: Merge libraries: per study -> all merged per organism
pipe_convert <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "covrle", e)))

  for (experiments in exp[done_exp]) {
    if (file.exists(report_failed_pipe_path(config, experiments[1]))) next
    try <- try({
      df_list <- lapply(experiments, function(e) read.experiment(e, validate = F))
      pipeline_convert_covRLE(df_list, config)
      pipeline_convert_bigwig(df_list, config)
    })
    report_failed_pipe(try, config, "bigwig_covRLE", experiments[1])
  }
}

pipe_counts <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "pcounts", e)))

  for (experiments in exp[done_exp]) {
    if (file.exists(report_failed_pipe_path(config, experiments[1]))) next
    try <- try({
      df_list <- lapply(experiments, function(e) read.experiment(e, validate = F))
      pipeline_count_table_psites(df_list, config)
    })
    report_failed_pipe(try, config, "counts", experiments[1])
  }
  return(invisible(NULL))
}

#' Create the superset collection of all samples of all organisms
#' Also saves a file in pipeline dir called ./metadata_done_samples.csv.
#' @inheritParams pipeline_merge_org
#' @param path_suffix character, default "". Add suffix to saved collection name, like a date to seperate version.
#' @return invisible(NULL)
#' @export
#' @examples
#' pipeline_collection_org(config, done_organisms = "Homo sapiens", path_suffix = "10_02_2025")
pipeline_collection_master <- function(config, pipelines = pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE),
                                    done_experiments = step_is_done_pipelines(config, "pcounts", pipelines),
                                    done_organisms = unique(names(done_experiments)),
                                    path_suffix = "", BPPARAM = MulticoreParam(workers = bpparam()$workers,
                                                                               exportglobals = FALSE,
                                                                               log = FALSE, force.GC = FALSE, fallback = T)) {
  message("- Collection of all samples of all organism")
  done_organisms <- unique(done_organisms)
  stopifnot(length(names(done_experiments)) > 0)

  exps_species <- done_experiments[names(done_experiments) %in% done_organisms]
  message("-- Organism: ", length(done_organisms), " (", length(exps_species), " studies)")
  df_list <- bplapply(exps_species, function(e)
    read.experiment(e, validate = FALSE), BPPARAM = BPPARAM)
  df <- do.call(rbind, df_list)

  libtype_df <- libraryTypes(df)
  if (length(libtype_df) != 1) stop("Only single libtype experiments supported for merging")
  if (libtype_df == "") stop("Libtype of experiment must be defined!")
  exp_name <- "all_samples-all_species"
  if (libtype_df != "RFP") exp_name <- paste0(exp_name, "_", libtype_df)
  message("Total samples done: ", nrow(df))

  out_dir <- file.path(config$config["bam"], exp_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fraction <- df$fraction
  #fraction <- gsub("auto_.*", "", df$fraction)
  study_size <- unlist(lapply(df_list, function(e) nrow(e)))
  studies <- rep(unlist(exps_species, use.names = FALSE), study_size)
  studies <- gsub("-.*", "", studies)
  fraction <- paste(fraction, studies, sep = "_")
  fraction <- gsub("^_", "",fraction)
  fraction <- gsub("^NA_", "",fraction)
  fraction <- paste0(fraction, "_", seq(length(fraction)))

  create.experiment(out_dir, exper = exp_name,
                    txdb = df@txdb, fa = df@fafile, organism = "all_species",
                    libtype = df$libtype,
                    condition = df$condition, rep = df$rep,
                    stage = df$stage, fraction = fraction,
                    files = df$filepath, result_folder = out_dir)
  message("Saved collection as: ", exp_name)
  message("--- Trying to load created collection")
  df <- read.experiment(exp_name, output.env = new.env())

  dt <- fread(config$complete_metadata)
  dt_complete <- dt[Run %in% runIDs(df),]
  dt_complete <- dt_complete[chmatch(Run, runIDs(df))]
  if (libtype_df == "RFP") {
    files <- filepath(df, "pshifted", base_folders = libFolder(df, mode = "all"))
    dt_complete[, is_pshifted := grepl("_pshifted\\.ofst$", files)]
    if (any(!dt_complete$is_pshifted)) {
      warning("Some samples are not pshifted even though they are flagged as 'complete', they are listed below:")
      print(dt_complete[is_pshifted == FALSE,])
    }
  }
  fwrite(dt_complete, file.path(config$project, "metadata_done_samples.csv"))
  message("Saved all complete studies metadata to: ", file.path(config$project, "metadata_done_samples.csv"))

  return(dt_complete)
}

#' Create the superset collection of all samples per organism
#' @inheritParams pipeline_merge_org
#' @param path_suffix character, default "". Add suffix to saved collection name, like a date to seperate version.
#' @return invisible(NULL)
#' @export
#' @examples
#' config <- pipeline_config()
#' pipelines <- pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE)
#' dirs <- lapply(pipelines, function(p) {
#' dirs_pipe <- massiveNGSpipe:::list_dirs_of_pipeline(pipeline = p)
#' names(dirs_pipe) <- tolower(names(dirs_pipe))
#' lapply(dirs_pipe, function(exp) {
#'   all(dir.exists(exp[-which(names(exp) %in% c("fastq", "trim", "exp"))]))
#' })
#' })
#' exists <- unlist(dirs, use.names = TRUE)
#' names(exists) <- sub("\\.", "-", names(exists))
#' names(exists) <- sub(" ", "_", names(exists))
#' sum(exists / length(exists))
#' all_exp_names <- unlist(massiveNGSpipe:::get_experiment_names(pipelines), use.names = FALSE)
#' stopifnot(all(names(exists) %in% all_exp_names))
#' done_experiments <- names(exists)[exists]
#' length(done_experiments)
#' names(done_experiments) <- stringr::str_to_sentence(sub("_", " ", sub(".*-", "", done_experiments)))
#' pipeline_collection_org(config, pipelines, done_organisms = "Homo sapiens", done_experiments, path_suffix = "_13_06_2025")
#' #pipeline_collection_org(config, done_organisms = "Homo sapiens", path_suffix = "_10_02_2025")
pipeline_collection_org <- function(config, pipelines = pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE),
                                    done_experiments = step_is_done_pipelines(config, "pcounts", pipelines),
                                    done_organisms = unique(names(done_experiments)),
                                    path_suffix = "", BPPARAM = MulticoreParam(workers = bpparam()$workers,
                                                                               exportglobals = FALSE,
                                                                               log = FALSE, force.GC = FALSE, fallback = T)) {
  message("- Collection of all samples per organism")
  done_organisms <- unique(done_organisms)
  stopifnot(length(names(done_experiments)) > 0)
  for (org in done_organisms) {
    exps_species <- done_experiments[names(done_experiments) == org]
    message("-- Organism: ", org, " (", length(exps_species), " studies)")
    df_list <- bplapply(exps_species, function(e)
      read.experiment(e, validate = FALSE), BPPARAM = BPPARAM)
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
    studies <- rep(unlist(exps_species, use.names = FALSE), study_size)
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
    message("--- Validating created collection")
    df <- read.experiment(exp_name, output.env = new.env())
    message("--- Merging count tables from studies")
    count_folder <- QCfolder(df)
    dir.create(count_folder, recursive = TRUE, showWarnings = FALSE)
    for (region in c("mrna", "cds", "leaders", "trailers")) {
      message("---- ", region)
      count_lists <- bplapply(df_list,
                             function(e, region) suppressMessages(countTable(e, region,
                                                           type = "summarized")),
                             region = region, BPPARAM = BPPARAM)
      count_list <- do.call(BiocGenerics::cbind, count_lists)
      saveName <- file.path(count_folder, paste0("countTable_", region, ".qs"))
      qs::qsave(count_list, file = saveName, nthreads = 5)
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
#' @param done_organisms Scientific name (e.g. Homo sapiens), which organisms to run merging for, default
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
    df <- df[df$libtype == libtype_to_merge,]

    # Overwrite default paths to merged
    libtype_df <- libraryTypes(df)

    df <- update_path_per_sample(df, libtype_df, libtype_to_merge = libtype_to_merge)

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
#' @param out_dir_base = ORFik::config()["bam"], will make directory internally
#' called name of exp (e.g. all_merged-Homo_sapiens_modalities)
#' @return invisible(NULL)
#' @export
pipeline_merge_org_modalities <- function(all_exp = list.experiments(validate = FALSE),
                                          out_dir_base = ORFik::config()["bam"]) {
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
    if (length(libtype_df) != length(df$libtype)) {
      message("Here are the duplicates:")
      tab <- table(df$libtype)
      print(tab[tab > 1])
      stop("Organism have > 1 all_merged track of the same type")
    }
    exp_name <- organism_merged_modalities_exp_name(org)
    create.experiment(NULL, exper = exp_name,
                      txdb = df@txdb,
                      libtype = libtype_df, fa = df@fafile,
                      organism = org, files = df$filepath,
                      result_folder = file.path(out_dir_base, exp_name))
    df <- read.experiment(exp_name)
    message("- Merging count tables")
    qc_folder <- QCfolder(df)
    regions <- c("mrna", "leaders", "cds", "trailers")
    dir.create(qc_folder, recursive = TRUE, showWarnings = FALSE)
    counts <- lapply(regions, function(region) {
      message("-- ", region)
      sum_exp <- do.call(cbind, lapply(seq(nrow(df)), function(i) {
        df_sub <- df[i,]
        df_sub@resultFolder <- libFolder(df[i,])
        countTable(df_sub, region = region, type = "summarized")
      }))
      save_RDSQS(sum_exp, file.path(qc_folder, paste0("countTable_", region, ".qs")))
    })
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
