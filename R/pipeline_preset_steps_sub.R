.datatable.aware <- TRUE # nolint

#' Download all SRA files for all studies
#'
#' Extract them into '<'accession'>'.fastq.gz or '<'accession'>'_{1,2}.fastq.gz
#' for SE/PE reads respectively.
#' @param pipeline a pipeline object
#' @param config the mNGSp config object from [pipeline_config]
pipeline_download <- function(pipeline, config) {
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (step_is_done(config, "fetch", conf["exp"])) next
        set_flag(config, "start", conf["exp"])
        download_sra(
            study[ScientificName == organism],
            conf["fastq"],
            compress = config$compress_raw_data
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
    config$BPPARAM_TRIM <- SerialParam()
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "trim", conf["exp"])) next
        index <- pipeline$organisms[[organism]]$index
        source_dir <- conf["fastq"]
        process_dir <- target_dir <- conf["bam"]
        trimmed_dir <- fs::path(process_dir, "trim")

        runs <- study[ScientificName == organism]
        # Files to run (Single end / Paired end)
        all_files <- run_files_organizer(runs, source_dir)
        # Trim
        barcodes_dt <-
        BiocParallel::bplapply(seq_len(nrow(runs)), function(i, all_files, runs) {
          run <- runs[i]$Run
          message(run)
          filenames <- all_files[[i]]
          file <- filenames[1]
          single_end <- is.na(filenames[2])
          file2 <- if(single_end) {NULL} else filenames[2]

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
          barcode_dt <- data.table()
          if (single_end && runs[i]$LIBRARYTYPE == "RFP") {
            barcode_dt <- barcode_detector_single(runs[i], source_dir, target_dir, trimmed_dir,
                                                  redownload_raw_if_needed = TRUE)
            if (barcode_dt$barcode_detected)  {
              barcode_dir <- file.path(trimmed_dir, "before_barcode_removal")
              dir.create(barcode_dir, showWarnings = FALSE, recursive = TRUE)
              file_trim <- massiveNGSpipe:::run_files_organizer(runs[i], trimmed_dir)[[1]]
              json_file <- grep(run, dir(trimmed_dir, "\\.json$", full.names = TRUE), value = TRUE)
              html_file <- grep(run, dir(trimmed_dir, "\\.html$", full.names = TRUE), value = TRUE)
              all_3_old_files <- c(file_trim, json_file, html_file)
              all_3_new_files <- gsub("/trimmed_|/trimmed2_", "/", file.path(barcode_dir, basename(all_3_old_files)))
              fs::file_move(all_3_old_files, all_3_new_files)
              file <- all_3_new_files[1]
              ORFik::STAR.align.single(
                file,
                output.dir = target_dir,
                adapter.sequence = adapter,
                index.dir = index, steps = "tr",
                trim.front = barcode_dt$barcode5p_size,
                trim.tail = barcode_dt$barcode3p_size
              )
              fs::file_delete(file)
            }
          }

          return(barcode_dt)
        }, all_files = all_files, runs = runs,
        BPPARAM = config$BPPARAM_TRIM)
        barcodes_dt <- rbindlist(barcodes_dt)
        fwrite(barcodes_dt, file.path(trimmed_dir, "adapter_barcode_table.csv"))

        set_flag(config, "trim", conf["exp"])
        if (config$delete_raw_files) fs::file_delete(unlist(all_files))
    }
}



#' Remove contaminants and align the reads to genome. Single and paired end
#' reads are handled separately, so the resulting logs are renamed with a
#' "_SINGLE" and/or "_PAIRED" suffix.
#' @inheritParams pipeline_download
pipeline_align <- function(pipeline, config) {
  # TODO: fix multiqc error for trimmed
  study <- pipeline$study
  did_contamint_removal <- "contam" %in% names(config$flag)
  did_collapse <- "collapsed" %in% names(config$flag)
  did_trim <- "trim" %in% names(config$flag)
  can_use_raw <- did_trim & !did_collapse & !config$delete_raw_files
  steps <- if(can_use_raw) {"tr"} else NULL
  steps <- c(steps, if(did_contamint_removal) {"co"} else NULL)
  steps <- paste(c(steps, "ge"), collapse = "-")
  for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    if (!step_is_next_not_done(config, "aligned", conf["exp"])) next
    index <- pipeline$organisms[[organism]]$index
    runs <- study[ScientificName == organism]
    trimmed_dir <- fs::path(conf["bam"], "trim")
    output_dir <- conf["bam"]

    keep.unaligned.genome <- config$keep.unaligned.genome
    if (any(runs$LibraryLayout == "SINGLE")) {
      input_dir <- ifelse(did_collapse,
                          fs::path(trimmed_dir, "SINGLE"),
                          trimmed_dir)
      stopifnot(dir.exists(input_dir))

      ORFik::STAR.align.folder(
        input.dir = input_dir,
        output.dir = output_dir, multiQC = FALSE,
        keep.unaligned.genome = keep.unaligned.genome,
        index.dir = index, steps = "ge", paired.end = FALSE
      )
      STAR.allsteps.multiQC(output_dir, steps = steps)
      for (stage in c("aligned")) {
        new_dir <- fs::path(output_dir, stage, "LOGS_SINGLE")
        if (dir.exists(new_dir)) {
          fs::dir_delete(new_dir)
        }
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
    }
    if (any(runs$LibraryLayout == "PAIRED")) {
      input_dir <- ifelse(did_collapse,
                          fs::path(trimmed_dir, "PAIRED"),
                          ifelse(can_use_raw, conf["fastq"], trimmed_dir))
      stopifnot(dir.exists(input_dir))
      # message("Paired end ignored for now, running collapsed pair mode only!")
      collapsed_paired_end_mode <- TRUE
      ORFik::STAR.align.folder(
        input.dir = input_dir,
        output.dir = output_dir,
        index.dir = index, steps = steps,
        resume = "ge",
        keep.unaligned.genome = keep.unaligned.genome,
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
    }

    if (nrow(runs))
    if (config$delete_collapsed_files)
      fs::dir_delete(input_dir)
    set_flag(config, "aligned", conf["exp"])
  }
}

#' Remove contaminants and align the reads to genome. Single and paired end
#' reads are handled separately, so the resulting logs are renamed with a
#' "_SINGLE" and/or "_PAIRED" suffix.
#' @inheritParams pipeline_download
pipeline_align_contaminants <- function(pipeline, config) {
  study <- pipeline$study
  for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    if (!step_is_next_not_done(config, "contam", conf["exp"])) next
    index <- pipeline$organisms[[organism]]$index
    runs <- study[ScientificName == organism]
    trimmed_dir <- fs::path(conf["bam"], "trim")
    output_dir <- conf["bam"]
    did_collapse <- "collapsed" %in% names(config$flag)
    keep.contaminants <- config$keep_contaminants
    if (any(runs$LibraryLayout == "SINGLE")) {
      input_dir <- ifelse(did_collapse,
                          fs::path(trimmed_dir, "SINGLE"),
                          trimmed_dir)
      ORFik::STAR.align.folder(
        input.dir = input_dir,
        output.dir = output_dir, keep.contaminants = keep.contaminants,
        index.dir = index, steps = "co", paired.end = FALSE
      )
      for (stage in c("contaminants_depletion")) {
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
        index.dir = index, steps = "co",
        keep.contaminants = keep.contaminants,
        paired.end = collapsed_paired_end_mode
      )
    }
    set_flag(config, "contam", conf["exp"])
  }
}

#' Remove all files apart from logs and final aligned BAMs.
#' Rename the BAMs into <run_accession>.bam format.
#' @inheritParams pipeline_download
pipeline_cleanup <- function(pipeline, config) {
    accession <- pipeline$accession
    study <- pipeline$study
    did_contamint_removal <- "contam" %in% names(config$flag)
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "cleanbam", conf["exp"])) next
        study_org <- study[ScientificName == organism,]
        bam_dir <- fs::path(conf["bam"], "aligned")
        if (did_contamint_removal) {
          fs::file_delete(fs::dir_ls(
            fs::path(conf["bam"], "contaminants_depletion"),
            glob = "**/*.out.*"
          ))
        }
        new_file_names <- fs::path(bam_dir, study_org$Run, ext = "bam")
        old_file_names <- match_bam_to_metadata(bam_dir, study_org, FALSE)

        file_names_to_delete <- new_file_names[new_file_names != old_file_names]
        file_names_to_delete <- file_names_to_delete[file.exists(file_names_to_delete)]
        if (length(file_names_to_delete) > 0) try(file.remove(file_names_to_delete), silent = TRUE)
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


pipeline_pshift <- function(df_list, config, accepted_lengths = config$accepted_lengths_rpf) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "pshifted", name(df))) next
    res <- tryCatch(shiftFootprintsByExperiment(df, output_format = "ofst",
                                                accepted.lengths = accepted_lengths),
                    error = function(e) {
                      message(e)
                      message(name(df))
                      message("Fix manually (skipping to next project!)")
                      bad_pshifting_report(df)
                      return(e)
                    })
    if(!inherits(res, "error")) {
      dir.create(QCfolder(df))
      set_flag(config, "pshifted", name(df))
    }
  }
}

#' Validate pshifting
#'
#' It will do these two things: Plot all TIS regions and
#' make a frame distribution table for all libs.
#' It then save a rds file called either 'warning.rds' or 'good.rds'.
#' Depending on the results was accepted or not:
#' \code{any(zero_frame$percent_length < 25)}
#' @param df_list a list of ORFik experiments
#' @param config a pipeline config
#' @return invisible(NULL)
pipeline_validate_shifts <- function(df_list, config) {
  #
  # TODO: This function can be improved, especially how we detect bad shifts!
  for (df in df_list) {
    if (!step_is_next_not_done(config, "valid_pshift", name(df))) next
    # Plot max 39 libraries!
    subset <- if (nrow(df) >= 40) {seq(39)} else {seq(nrow(df))}
    invisible(shiftPlots(df[subset,], output = "auto", plot.ext = ".png"))
    # Check frame usage
    frameQC <- orfFrameDistributions(df)
    remove.experiments(df)
    zero_frame <- frameQC[frame == 0 & best_frame == FALSE,]
    # Store a flag that says good / bad shifting
    QCFolder <- QCfolder(df)
    data.table::fwrite(frameQC, file = file.path(QCFolder, "Ribo_frames_all.csv"))
    data.table::fwrite(zero_frame, file = file.path(QCFolder, "Ribo_frames_badzero.csv"))
    if (any(zero_frame$percent_length < 25)) {
      warning("Some libraries contain shift that is < 25% of CDS coverage")
      saveRDS(FALSE, file.path(QCFolder, "warning.rds"))
    } else saveRDS(TRUE, file.path(QCFolder, "good.rds"))
    set_flag(config, "valid_pshift", name(df))
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
