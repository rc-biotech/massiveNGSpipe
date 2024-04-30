.datatable.aware <- TRUE # nolint

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
        # Files to run (Single end / Paired end)
        all_files <- run_files_organizer(runs, source_dir)
        # Trim
        BiocParallel::bplapply(seq_len(nrow(runs)), function(i, all_files, runs) {
            message(runs[i]$Run)
            filenames <- all_files[[i]]
            file <- filenames[1]
            file2 <- if(is.na(filenames[2])) {NULL} else filenames[2]

            if (!grepl("\\.fasta$|\\.fasta\\.gz$", file)) {
              adapter <- try(fastqc_adapters_info(file))
              if (is(adapter, "try-error")) {
                message("This is a fasta file, fastqc adapter detection disabled")
                adapter <- "disable"
              }
            } else adapter <- "disable"

            ORFik::STAR.align.single(
              file, file2,
              output.dir = target_dir,
              adapter.sequence = adapter,
              index.dir = index, steps = "tr"
            )

        }, all_files = all_files, runs = runs,
        BPPARAM = BiocParallel::MulticoreParam(8))

        set_flag(config, "trim", conf["exp"])
        if (config$delete_raw_files) fs::file_delete(unlist(all_files))
    }
}



#' Remove contaminants and align the reads to genome. Single and paired end
#' reads are handled separately, so the resulting logs are renamed with a
#' "_SINGLE" and/or "_PAIRED" suffix.
#' @param pipeline a pipeline object
#' @param config the mNGSp config object
pipeline_align <- function(pipeline, config) {
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "aligned", conf["exp"])) next
        index <- pipeline$organisms[[organism]]$index
        runs <- study[ScientificName == organism]
        trimmed_dir <- fs::path(conf["bam"], "trim")
        output_dir <- conf["bam"]
        did_collapse <- "collapsed" %in% names(config$flag)
        if (any(runs$LibraryLayout == "SINGLE")) {
          input_dir <- ifelse(did_collapse,
                              fs::path(trimmed_dir, "SINGLE"),
                              trimmed_dir)
            ORFik::STAR.align.folder(
                input.dir = input_dir,
                output.dir = output_dir,
                index.dir = index, steps = "co-ge", paired.end = FALSE
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
            if (config$delete_collapsed_files)
              fs::dir_delete(input_dir)
        }
        if (any(runs$LibraryLayout == "PAIRED")) {
          input_dir <- ifelse(did_collapse,
                              fs::path(trimmed_dir, "PAIRED"),
                              trimmed_dir)
            # message("Paired end ignored for now, running collapsed pair mode only!")
            collapsed_paired_end_mode <- TRUE
            ORFik::STAR.align.folder(
                input.dir = input_dir,
                output.dir = output_dir,
                index.dir = index, steps = "tr-co-ge",
                resume = "co",
                paired.end = collapsed_paired_end_mode
            )
            # for (stage in c("contaminants_depletion", "aligned")) {
            #     fs::file_move(
            #         fs::path(output_dir, stage, "LOGS"),
            #         fs::path(output_dir, stage, "LOGS_PAIRED")
            #     )
            # }
            # for (filename in c("full_process.csv", "runCommand.log")) {
            #     fs::file_move(
            #         fs::path(output_dir, filename),
            #         fs::path(output_dir, paste0(
            #             fs::path_ext_remove(filename), "_PAIRED.",
            #             fs::path_ext(filename)
            #         ))
            #     )
            # }
            if (config$delete_collapsed_files)
              fs::dir_delete(input_dir)
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
        study_org <- study[ScientificName == organism,]
        bam_dir <- fs::path(conf["bam"], "aligned")

        fs::file_delete(fs::dir_ls(
            fs::path(conf["bam"], "contaminants_depletion"),
            glob = "**/*.out.*"
        ))
        old_file_names <- match_bam_to_metadata(bam_dir, study_org, FALSE)
        new_file_names <- fs::path(bam_dir, study_org$Run, ext = "bam")
        stopifnot(length(old_file_names) == length(new_file_names))
        for (i in seq_along(old_file_names)) {
          fs::file_move(old_file_names[i], new_file_names[i])
        }
        set_flag(config, "cleanbam", conf["exp"])
    }
}

#' Create ORFik experiment from config study and bam files
#' @noRd
pipeline_create_experiment <- function(pipeline, config) {
    df_list <- list()
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        experiment <- conf["exp"]
        if (!step_is_next_not_done(config, "exp", experiment)) {
          if (!step_is_done(config, "exp", conf["exp"])) return(NULL)
          df_list <- c(df_list, read.experiment(conf["exp"],
                                                output.env = new.env()))
          next
        }
        annotation <- pipeline$organisms[[organism]]$annotation
        study <- pipeline$study
        stopifnot(nrow(study) > 0)
        study <- study[ScientificName == organism,]
        if (nrow(study) == 0)
          stop("No samples for organism wanted in study!")


        metadata_clean <- cleanup_metadata_for_exp(study)
        bam_dir <- fs::path(conf["bam"], "aligned")
        bam_files <- match_bam_to_metadata(bam_dir, study, metadata_clean$paired_end)
        ORFik::create.experiment(
            dir = bam_dir,
            exper = experiment, txdb = paste0(annotation["gtf"], ".db"),
            libtype = study$LIBRARYTYPE,
            fa = annotation["genome"], organism = organism,
            stage = metadata_clean$stage, rep = study$REPLICATE,
            condition = metadata_clean$condition,
            fraction = metadata_clean$fraction,
            pairedEndBam = metadata_clean$paired_end,
            author = unique(study$AUTHOR),
            files = bam_files, runIDs = study$Run
        )
        df <- ORFik::read.experiment(experiment,
                                     output.env = new.env())
        set_flag(config, "exp", conf["exp"])
        df_list <- c(df_list, df)
    }
  return(df_list)
}

pipeline_create_ofst <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "ofst", name(df))) next
    convert_bam_to_ofst(df)
    set_flag(config, "ofst", name(df))
  }
}

pipeline_pshift <- function(df_list, accepted.length = c(20, 21, 25:33),
                            config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "pshifted", name(df))) next
    res <- tryCatch(shiftFootprintsByExperiment(df, output_format = "ofst",
                                                accepted.lengths = accepted.length),
                    error = function(e) {
                      message(e)
                      message("Fix manually (skipping to next project!)")
                      return(e)
                    })
    if(!inherits(res, "error")) {
      dir.create(QCfolder(df))
      set_flag(config, "pshifted", name(df))
    } else {
      bad_pshifting_report(df)
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
