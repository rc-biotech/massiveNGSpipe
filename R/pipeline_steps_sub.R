.datatable.aware <- TRUE # nolint


#' Create a pipeline for the specified study.
#'
#' @param study a data.table of output from ORFik::download.SRA.metadata
#' @param study_accession any accession accepted by
#' \code{ORFik::download.SRA.metadata}. It is needed to know
#' if you used GEO, PRJ, SRP etc as the search query.
#' @param config Configured directories for pipeline as a list
#' @return a list of pipeline objects, one for each study,
#' subsetted by organism per study.
pipeline_init <- function(study, study_accession, config) {
    # For each organism in the study, create an ORFik experiment config,
    # fetch reference genome and contaminants, and create an index.
    organisms <- list()
    for (organism in unique(study$ScientificName)) {
        message("---- ", organism)
        assembly_name <- gsub(" ", "_", trimws(tolower(organism)))
        experiment <- paste(study_accession, assembly_name, sep = "-")

        # Create ORFik experiment config manually, without separate folders
        # for different library strategies.
        conf <- config.exper(experiment, assembly_name, "", config[["config"]])
        # Fix bad naming from config.exper
        conf <- gsub("//", "/", conf); conf <- gsub("_$", "", conf)
        names(conf) <- gsub(" $", "", names(conf))
        sapply(conf[1:3], fs::dir_create)

        # Will now swap to toplevel if needed by itself
        # Will swap to refseq if ensembl fails
        ensembl_db <- try(
          annotation <- ORFik::getGenomeAndAnnotation(
            organism = organism,
            genome = TRUE, GTF = TRUE,
            phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
            output.dir = conf["ref"],
            assembly_type = "primary_assembly", optimize = TRUE,
            pseudo_5UTRS_if_needed = 100, notify_load_existing = FALSE,
            gene_symbols = TRUE
          )
        )
        if (is(ensembl_db, "try-error")) {
          message("-- Genome download failed for: ", organism, " on Ensembl")
          message("Trying refseq")
          annotation <- ORFik::getGenomeAndAnnotation(
            organism = organism,
            genome = TRUE, GTF = TRUE,
            phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
            output.dir = conf["ref"],
            assembly_type = "primary_assembly", optimize = TRUE,
            pseudo_5UTRS_if_needed = 100, db = "refseq",
            notify_load_existing = FALSE, gene_symbols = TRUE
          )
        }
        index <- ORFik::STAR.index(annotation, notify_load_existing = FALSE)
        organisms[[organism]] <- list(
            conf = conf, annotation = annotation, index = index,
            experiment = experiment
        )
    }
    return(list(
        accession = study_accession, study = study, organisms = organisms
    ))
}

#' Download all SRA files for all studies
#'
#' Extract them into '<'accession'>'.fastq.gz or '<'accession'>'_{1,2}.fastq.gz
#' for SE/PE reads respectively.
#' @param pipeline a pipeline object
pipeline_download <- function(pipeline, config) {
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (step_is_done(config, "fetch", conf["exp"])) next
        set_flag(config, "start", conf["exp"])
        download_sra(
            study[ScientificName == organism],
            conf["fastq"],
            compress = FALSE
        )
        set_flag(config, "fetch", conf["exp"])
    }
}

#' Trim all .fastq.gz files in a given pipeline. Paired end reads
#' ("SRR<...>_1.fastq.gz", "SRR<...>_2.fastq.gz") are trimmed separately
#' from each other to allow for differing adapters.
#' @param pipeline a pipeline object
pipeline_trim <- function(pipeline, config) {
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "trim", conf["exp"])) next
        index <- pipeline$organisms[[organism]]$index
        source_dir <- conf["fastq"]
        target_dir <- conf["bam"]
        runs <- study[ScientificName == organism]
        for (i in seq_len(nrow(runs))) {
            message(runs[i]$Run)
            filenames <-
                if (runs[i]$LibraryLayout == "PAIRED") {
                  paste0(runs[i]$Run, c("_1", "_2"), ".fastq")
                } else {
                  paste0(runs[i]$Run, ".fastq")
                }
            filenames <- fs::path(source_dir, filenames)
            if (!all(file.exists(filenames))) {
              filenames <- paste0(filenames, ".gz")
              if (!all(file.exists(filenames))) {
                stop("File does not exist: ", filenames[1])
              }
            }
            file <- filenames[1]
            file2 <- if(is.na(filenames[2])) {NULL} else filenames[2]

            ORFik::STAR.align.single(
              file, file2,
              output.dir = target_dir,
              adapter.sequence = fastqc_adapters_info(file),
              index.dir = index, steps = "tr"
            )
            if (config$delete_raw_files) fs::file_delete(file)
        }
        set_flag(config, "trim", conf["exp"])
    }
}

#' Collapse trimmed single-end reads and move them into "trim/SINGLE"
#' subdirectory. PE reads are moved into "trim/PAIRED" without modification.
#' @param pipeline a pipeline object
pipeline_collapse <- function(pipeline, config) {
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "collapsed", conf["exp"])) next
        trimmed_dir <- fs::path(conf["bam"], "trim")
        runs_paired <- study[ScientificName == organism &
            LibraryLayout == "PAIRED"]
        runs_single <- study[ScientificName == organism &
            LibraryLayout != "PAIRED"]
        for (i in seq_len(nrow(runs_paired))) {
            outdir <- fs::path(trimmed_dir, runs_paired[i]$LibraryLayout)
            fs::dir_create(outdir)
            filenames <- paste0(
                "trimmed_", runs_paired[i]$Run, c("_1", "_2"), ".fastq"
            )
            for (filename in filenames) {
                fs::file_move(fs::path(trimmed_dir, filename), outdir)
            }
        }

        for (i in seq_len(nrow(runs_paired))) {
          outdir <- fs::path(trimmed_dir, runs_paired[i]$LibraryLayout)
          fs::dir_create(outdir)
          filenames <- paste0(
            "trimmed_", runs_paired[i]$Run, c("_1", "_2"), ".fastq"
          )
          for (filename in filenames) {
            fs::file_move(fs::path(trimmed_dir, filename), outdir)
          }
        }
        outdir <- fs::path(trimmed_dir, "SINGLE")
        fs::dir_create(outdir)
        BiocParallel::bplapply(runs_single$Run, function(srr) {
          filename <- paste0("trimmed_", srr, ".fastq")
          ORFik::collapse.fastq(
            fs::path(trimmed_dir, filename), outdir,
            compress = TRUE
          )
          fs::file_delete(fs::path(trimmed_dir, filename))
        }, BPPARAM = BiocParallel::MulticoreParam(16))
        set_flag(config, "collapsed", conf["exp"])
    }
}

#' Remove contaminants and align the reads to genome. Single and paired end
#' reads are handled separately, so the resulting logs are renamed with a
#' "_SINGLE" and/or "_PAIRED" suffix.
#' @param pipeline a pipeline object
pipeline_align <- function(pipeline, config) {
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "aligned", conf["exp"])) next
        index <- pipeline$organisms[[organism]]$index
        runs <- study[ScientificName == organism]
        trimmed_dir <- fs::path(conf["bam"], "trim")
        output_dir <- conf["bam"]
        if (any(runs$LibraryLayout == "SINGLE")) {
            ORFik::STAR.align.folder(
                input.dir = fs::path(trimmed_dir, "SINGLE"),
                output.dir = output_dir,
                index.dir = index, steps = "co-ge", paired.end = FALSE,
            )
            for (stage in c("contaminants_depletion", "aligned")) {
                fs::file_move(
                    fs::path(output_dir, stage, "LOGS"),
                    fs::path(output_dir, stage, "LOGS_SINGLE")
                )
            }
            for (filename in c("full_process.csv", "runCommand.log")) {
                fs::file_move(
                    fs::path(output_dir, filename),
                    fs::path(output_dir, paste0(
                        fs::path_ext_remove(filename), "_SINGLE.",
                        fs::path_ext(filename)
                    ))
                )
            }
            fs::dir_delete(fs::path(trimmed_dir, "SINGLE"))
        }
        if (any(runs$LibraryLayout == "PAIRED")) {
            ORFik::STAR.align.folder(
                input.dir = fs::path(trimmed_dir, "PAIRED"),
                output.dir = output_dir,
                index.dir = index, steps = "co-ge", paired.end = TRUE,
            )
            for (stage in c("contaminants_depletion", "aligned")) {
                fs::file_move(
                    fs::path(output_dir, stage, "LOGS"),
                    fs::path(output_dir, stage, "LOGS_PAIRED")
                )
            }
            for (filename in c("full_process.csv", "runCommand.log")) {
                fs::file_move(
                    fs::path(output_dir, filename),
                    fs::path(output_dir, paste0(
                        fs::path_ext_remove(filename), "_PAIRED.",
                        fs::path_ext(filename)
                    ))
                )
            }
            fs::dir_delete(fs::path(trimmed_dir, "PAIRED"))
        }
        set_flag(config, "aligned", conf["exp"])
    }
}

#' Remove all files apart from logs and final aligned BAMs.
#' Rename the BAMs into <run_accession>.bam format.
#' @param pipeline a pipeline object
pipeline_cleanup <- function(pipeline, config) {
    accession <- pipeline$accession
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "cleanbam", conf["exp"])) next
        runs <- study[ScientificName == organism]$Run
        fs::file_delete(fs::dir_ls(
            fs::path(conf["bam"], "contaminants_depletion"),
            glob = "**/*.out.*"
        ))
        for (run in runs) {
            fs::file_move(
                fs::dir_ls(
                    fs::path(conf["bam"], "aligned"),
                    glob = paste0("**/*", run, "*.out.bam")
                ),
                fs::path(conf["bam"], "aligned", run, ext = "bam")
            )
        }
        set_flag(config, "cleanbam", conf["exp"])
    }
}

#' Create ORFik experiment from config study and bam files
pipeline_create_experiment <- function(pipeline, config) {
    df_list <- list()
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "exp", conf["exp"])) {
          if (!step_is_done(config, "exp", conf["exp"])) return(NULL)
          df_list <- c(df_list, read.experiment(conf["exp"],
                                                output.env = new.env()))
          next
        }
        annotation <- pipeline$organisms[[organism]]$annotation
        assembly_name <- gsub(" ", "_", trimws(tolower(organism)))
        accession <- pipeline$accession
        study <- pipeline$study
        stopifnot(nrow(study) > 0)
        study <- study[ScientificName == organism,]
        if (nrow(study) == 0)
          stop("No samples for organism wanted in study!")
        experiment <- paste(accession, assembly_name, sep = "-")
        # Do some small correction to info and merge
        remove <- "^_|_$|^NA_|_NA$|^NA$|^_$|^__$|^___$"
        # Stage
        stage <- paste0(study$CELL_LINE, "_",study$TISSUE)
        stage <- gsub(paste0(remove, "|NONE_|_NONE"), "", stage)
        stage <- gsub(paste0(remove, "|NONE_|_NONE"), "", stage) # Twice
        # Condition
        condition <- paste0(study$CONDITION)
        condition <- gsub(remove, "", condition)
        condition <- gsub(remove, "", condition)
        # Fraction
        fraction <- paste(study$FRACTION,study$TIMEPOINT, study$BATCH, sep = "_")
        if (!is.null(study$GENE)) fraction <- paste(fraction, study$GENE, sep = "_")
        fraction <- gsub(remove, "", fraction)
        study$INHIBITOR[is.na(study$INHIBITOR)] <- ""
        add_inhibitor_to_fraction <-
          !all(study$INHIBITOR %in% c("chx", "CHX"))
        if (add_inhibitor_to_fraction) {
          fraction <- paste0(fraction, "_",study$INHIBITOR)
        }
        fraction <- gsub(remove, "", fraction); fraction <- gsub(remove, "", fraction)
        fraction <- gsub(remove, "", fraction); fraction <- gsub(remove, "", fraction)
        # PAIRED END
        paired_end <- study$LibraryLayout == "PAIRED"
        if (any(paired_end)) stop("TODO: Fix for paired end, will not work for now")

        bam_dir <- fs::path(conf["bam"], "aligned")
        bam_files <- match_bam_to_metadata(bam_dir, study, paired_end)
        ORFik::create.experiment(
            dir = bam_dir,
            exper = experiment, txdb = paste0(annotation["gtf"], ".db"),
            libtype = "RFP",  fa = annotation["genome"], organism = organism,
            stage = stage, rep = study$REPLICATE,
            condition = condition, fraction = fraction,
            author = unique(study$AUTHOR),
            files = bam_files
        )
        df <- ORFik::read.experiment(experiment,
                                     output.env = new.env())
        # TODO: Copy experiment over
        # df[-(1:4),3] <- sapply(
        #    stringr::str_split(fs::path_file(df[-(1:4),6]), "_"),
        #    function(x) x[3]
        # )
        set_flag(config, "exp", conf["exp"])
        df_list <- c(df_list, df)
    }
  return(df_list)
}

pipeline_create_ofst <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "ofst", name(df))) next
    convertLibs(df)
    remove.experiments(df)
    set_flag(config, "ofst", name(df))
  }
}

pipeline_pshift <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "pshifted", name(df))) next
    res <- tryCatch(shiftFootprintsByExperiment(df, output_format = "ofst",
                                                accepted.lengths = c(20, 21, 25:33)),
                    error = function(e) {
                      message(e)
                      message("Fix manually (skipping to next project!)")
                      return(e)
                    })
    if(!inherits(res, "error")) {
      dir.create(QCfolder(df))
      set_flag(config, "pshifted", name(df))
    }
  }
}

pipeline_convert_covRLE <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "covrle", name(df))) next
    convert_to_covRleList(df)
    set_flag(config, "covrle", name(df))
  }
}

pipeline_convert_bigwig <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "bigwig", name(df))) next
    convert_to_bigWig(df)
    set_flag(config, "bigwig", name(df))
  }
}

pipeline_count_table_psites <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "pcounts", name(df))) next
    ORFik::countTable_regions(df, lib.type = "pshifted", forceRemake = TRUE)
    set_flag(config, "pcounts", name(df))
  }
}
