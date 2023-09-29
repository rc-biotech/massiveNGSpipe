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
        all_files <-
        lapply(seq_len(nrow(runs)), function(i) {
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
              stop("File does not exist to trim (both .gz and unzipped): ",
                   filenames[1])
            }
          }
          file <- filenames[1]
          file2 <- if(is.na(filenames[2])) {NULL} else filenames[2]
          return(c(file, file2))
        })
        # Trim
        BiocParallel::bplapply(seq_len(nrow(runs)), function(i, all_files, runs) {
            message(runs[i]$Run)
            filenames <- all_files[[i]]
            file <- filenames[1]
            file2 <- if(is.na(filenames[2])) {NULL} else filenames[2]

            ORFik::STAR.align.single(
              file, file2,
              output.dir = target_dir,
              adapter.sequence = fastqc_adapters_info(file),
              index.dir = index, steps = "tr"
            )

        }, all_files = all_files, runs = runs,
        BPPARAM = BiocParallel::MulticoreParam(8))

        set_flag(config, "trim", conf["exp"])
        if (config$delete_raw_files) fs::file_delete(unlist(all_files))
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
        any_paired_libs <- nrow(runs_paired) > 0
        if (any_paired_libs) {
            outdir <- fs::path(trimmed_dir, runs_paired[1]$LibraryLayout)
            fs::dir_create(outdir)
            filenames <- paste0(
                "trimmed" , c("1_", "2_"), runs_paired[i]$Run, "_1", ".fastq"
            )
            message("Paired end data is now collapsed into 1 file,
                    and 2nd file is reverse complimented before merging!")
            full_filenames <- file.path(trimmed_dir, filenames)
            BiocParallel::bplapply(full_filenames, function(filename) {
              ORFik::collapse.fastq(
                filename, outdir,
                compress = TRUE
              )
              fs::file_delete(filename)
            }, BPPARAM = BiocParallel::MulticoreParam(16))
            # Read in and reverse second file, then merge back into 1.
            full_filenames <- file.path(outdir, paste0("collapsed_", filenames))
            first_files <- filenames[seq_along(filenames) %% 2 == 1]
            second_files <- filenames[seq_along(filenames) %% 2 == 0]

            for (i in seq_along(first_files)) {
              first_file <- first_files[(i*2)-1]
              second_file <- second_files[(i*2)]
              a <- readDNAStringSet(first_file, format = "fasta", use.names = TRUE)
              b <- readDNAStringSet(second_file, format = "fasta", use.names = TRUE)
              b <- reverseComplement(b)
              writeXStringSet(c(a,b), first_file, format = "fasta")
              fs::file_delete(second_file)
            }
        }

        any_single_libs <- nrow(runs_single) > 0
        if (any_single_libs) {
          outdir <- fs::path(trimmed_dir, "SINGLE")
          fs::dir_create(outdir)
          BiocParallel::bplapply(runs_single$Run, function(srr) {
            filename <- list.files(trimmed_dir, paste0(srr, "\\."),
                                   full.names = TRUE)
            if (length(filename) != 1) {
              filename <- file.path(trimmed_dir,
                                    paste0("trimmed_", srr, ".fastq"))
            }
            ORFik::collapse.fastq(
              filename, outdir,
              compress = TRUE
            )
            fs::file_delete(filename)
          }, BPPARAM = BiocParallel::MulticoreParam(16))
        }

        set_flag(config, "collapsed", conf["exp"])
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
            if (config$delete_collapsed_files)
              fs::dir_delete(fs::path(trimmed_dir, "SINGLE"))
        }
        if (any(runs$LibraryLayout == "PAIRED")) {
            message("Paired end ignored for now, running collapsed pair mode only!")
            collapsed_paired_end_mode <- TRUE
            ORFik::STAR.align.folder(
                input.dir = fs::path(trimmed_dir, "PAIRED"),
                output.dir = output_dir,
                index.dir = index, steps = "co-ge",
                paired.end = collapsed_paired_end_mode,
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

        # Do some small correction to info and merge
        remove <- "^_|_$|^NA_|_NA$|^NA$|^_$|^__$|^___$"
        # Stage
        stage <- paste0(study$CELL_LINE, "_",study$TISSUE)
        stage <- gsub(paste0(remove, "|NONE_|_NONE"), "", stage)
        stage <- gsub(paste0(remove, "|NONE_|_NONE"), "", stage) # Twice
        # Condition
        condition <- paste0(study$CONDITION)
        if (!is.null(study$GENE)) condition <- paste(condition, study$GENE, sep = "_")
        condition <- gsub(remove, "", condition)
        condition <- gsub(remove, "", condition)
        condition <- gsub(remove, "", condition)
        condition <- gsub(remove, "", condition)
        # Fraction
        fraction <- paste(study$FRACTION,study$TIMEPOINT, study$BATCH, sep = "_")
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
        if (any(paired_end)) {
          message("Only running single end for now, make fix for this to work normally")
          paired_end <- FALSE
        }

        bam_dir <- fs::path(conf["bam"], "aligned")
        bam_files <- match_bam_to_metadata(bam_dir, study, paired_end)
        ORFik::create.experiment(
            dir = bam_dir,
            exper = experiment, txdb = paste0(annotation["gtf"], ".db"),
            libtype = "RFP",  fa = annotation["genome"], organism = organism,
            stage = stage, rep = study$REPLICATE,
            condition = condition, fraction = fraction,
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
