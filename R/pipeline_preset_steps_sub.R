
.datatable.aware <- TRUE # nolint

#' Download all SRA files for all studies
#'
#' Extract them into \verb{<accession>.fastq.gz} or \verb{<accession>_\{1,2\}.fastq.gz}
#' for SE/PE reads respectively.
#'
#' The per-run downloads run in a subprocess (aws s3 sync / ascp /
#' fasterq-dump's own console output does not otherwise get captured, see
#' run_experiment_subprocess()) that this function polls for progress,
#' rendering the checklist on each tick.
#' @param pipeline a pipeline object
#' @param config the mNGSp config object from [pipeline_config]
#' @param pipelines the full pipelines list, used only to render the
#' checklist (study-level context beyond this one pipeline/organism);
#' defaults to just this pipeline if called standalone.
#' @return invisible(NULL)
pipeline_download <- function(pipeline, config, pipelines = list(pipeline)) {
    study <- pipeline$study
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (step_is_done(config, "fetch", conf["exp"])) next
        set_flag(config, "start", conf["exp"])
        experiment <- conf["exp"]
        info <- study[ScientificName == organism]

        # Resume support: skip runs already marked done for this
        # experiment AND whose actual fastq output still validates on
        # disk -- the marker alone isn't trusted here, because
        # download_sra() has no idempotency of its own and can even
        # delete a stale-compression leftover file (see
        # fastq_output_exists_and_valid(), download_sra.R).
        done_runs <- samples_done(config, "fetch", experiment)
        already_done <- vapply(info$Run, function(run) {
          run %in% done_runs &&
            fastq_output_exists_and_valid(
              run, conf["fastq"],
              info$LibraryLayout[info$Run == run] == "PAIRED",
              config$compress_raw_data)
        }, logical(1))
        info <- info[!already_done]

        if (nrow(info) > 0) {
          run_experiment_subprocess(
            func = function(info, outdir, compress, config, experiment) {
              download_sra(
                info, outdir, compress = compress,
                after_run = function(run) set_sample_flag(config, "fetch", experiment, run)
              )
            },
            args = list(info = info, outdir = conf["fastq"],
                        compress = config$compress_raw_data, config = config, experiment = experiment),
            logfile_out = file.path(pipeline_log_base(config), "console", "fetch", paste0(experiment, ".out.log")),
            logfile_err = file.path(pipeline_log_base(config), "console", "fetch", paste0(experiment, ".err.log")),
            on_poll = function() pipeline_checklist(pipelines, config)
          )
        }
        set_flag(config, "fetch", conf["exp"])
    }
}

#' Trim all .fastq.gz files in a given pipeline. Paired end reads
#' ("SRR<...>_1.fastq.gz", "SRR<...>_2.fastq.gz") are trimmed separately
#' from each other to allow for differing adapters.
#'
#' The per-run work runs in a subprocess (fastp's own console output does
#' not otherwise get captured, see run_experiment_subprocess()) that this
#' function polls for progress, rendering the checklist on each tick.
#' @inheritParams pipeline_download
#' @param pipelines the full pipelines list, used only to render the
#' checklist (study-level context beyond this one pipeline/organism);
#' defaults to just this pipeline if called standalone.
pipeline_trim <- function(pipeline, config, pipelines = list(pipeline)) {
    study <- pipeline$study
    config$BPPARAM_TRIM <- bpparam_from_config(config, "trim")
    for (organism in names(pipeline$organisms)) {
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "trim", conf["exp"])) next
        index <- pipeline$organisms[[organism]]$index
        source_dir <- conf["fastq"]
        process_dir <- target_dir <- conf["bam"]
        trimmed_dir <- fs::path(process_dir, "trim")

        runs_full <- study[ScientificName == organism]
        experiment <- conf["exp"]

        # Resume support: skip runs already marked done for this
        # experiment. Trusted directly here (unlike fetch) -- the marker
        # is only ever written after detect_adapter_and_trim()/
        # run_barcode_detection_and_trim() return successfully for that
        # one run, and re-running a not-yet-done run naturally overwrites
        # its own output files, so no separate output-existence check is
        # needed.
        #
        # Filter to the not-yet-done subset BEFORE resolving files (not
        # after): run_files_organizer() requires every row it's given to
        # resolve to a real file, but an already-trimmed sibling's raw
        # fastq is typically already gone (config$delete_raw_files is the
        # online-mode default) -- resolving against the full study here
        # would error on that missing file even when only reprocessing
        # one already-done study's single fixed sample. Confirmed live:
        # this previously errored instantly, correctly got caught by
        # pipe_trim_collapse()'s try(), but because the "trim" flag could
        # then never be set, the step's own retry-skip-forever design
        # continued permanently and silently (warnings deferred in a
        # never-returning top-level call) -- hours of near-zero-CPU
        # "nothing happening", easy to misread as a BiocParallel/fork
        # deadlock when it is neither.
        done_runs <- samples_done(config, "trim", experiment)
        not_done <- !(runs_full$Run %in% done_runs)
        runs <- runs_full[not_done]
        all_files <- if (nrow(runs) > 0) run_files_organizer(runs, source_dir) else list()

        if (nrow(runs) > 0) {
          run_experiment_subprocess(
            func = function(all_files, runs, mode, target_dir, trimmed_dir, source_dir, config, experiment) {
              lapply(seq_len(nrow(runs)), function(i) {
                study_sample <- runs[i]
                run <- study_sample$Run
                message(run)
                filenames <- all_files[[i]]
                single_end <- is.na(filenames[2])
                file <- filenames[1]
                file2 <- if(!single_end) filenames[2]
                barcode_dt <- data.table()
                check_for_barcodes <- runs[i]$LIBRARYTYPE == "RFP"
                if (check_for_barcodes) {
                  # Detects barcode sizes from a cheap subsample, then
                  # does the ONE real full adapter+barcode trim pass --
                  # NOT detect_adapter_and_trim() first (that used to run
                  # a full adapter-only pass which this then reprocessed
                  # a second time in full; see run_barcode_detection_and_trim()).
                  barcode_dt <- run_barcode_detection_and_trim(study_sample, source_dir,
                                                               target_dir, trimmed_dir, mode,
                                                               file, file2)
                } else {
                  detect_adapter_and_trim(file, target_dir, file2)
                }
                # Store the actual barcode_dt row (not just TRUE) so it
                # can be reconstructed below even across a resumed run.
                set_sample_flag(config, "trim", experiment, run, value = barcode_dt)
                barcode_dt
              })
            },
            args = list(all_files = all_files, runs = runs, mode = config$mode,
                        target_dir = target_dir, trimmed_dir = trimmed_dir,
                        source_dir = source_dir, config = config, experiment = experiment),
            logfile_out = file.path(pipeline_log_base(config), "console", "trim", paste0(experiment, ".out.log")),
            logfile_err = file.path(pipeline_log_base(config), "console", "trim", paste0(experiment, ".err.log")),
            on_poll = function() pipeline_checklist(pipelines, config)
          )
        }

        # Reconstruct the full barcode table across both samples done in
        # this run and any completed in an earlier, interrupted attempt.
        # A handful of on-disk "trim" markers predate this reconstruction
        # existing at all and store plain TRUE rather than a barcode_dt
        # row (confirmed live, PRJNA750456-mus_musculus, 2026-10-06: 2 of
        # 20 markers) -- rbindlist() errors ("Item N of input is not a
        # data.frame...") on any such entry, not just ones this package's
        # own backfill writes (already fixed there separately, see
        # apply_trim_fix_to_sample()).
        #
        # A plain TRUE marker must NOT simply be dropped here (an earlier
        # version of this fix did exactly that, reasoning the row's real
        # detection detail was unrecoverable either way) -- this table is
        # read back by later code as a "what's already been
        # detected/trimmed" cache, so silently dropping a sample's row
        # makes it look never-processed and triggers a wasteful full
        # redo for it next time. Confirmed live, 2026-10-07,
        # PRJNA637713-zea_mays (296 samples) and
        # PRJEB36473-schizosaccharomyces_pombe (12 samples): after a
        # single-sample fix, the rebuilt table had ONLY that one sample's
        # row -- every other (legacy-marker) sibling's row vanished from
        # the file entirely. The real per-sample detail for a legacy
        # TRUE marker genuinely isn't recoverable from the marker alone,
        # but it typically still exists in THIS table from before this
        # rebuild overwrites it (the table was written wholesale by
        # older code, predating per-sample markers) -- reuse it from
        # there, matching apply_trim_fix_to_sample()'s own backfill
        # logic exactly. Only synthesize an id-only placeholder row (not
        # drop the sample) when no prior row exists anywhere.
        trim_marker_values <- sample_flag_values(config, "trim", experiment)
        existing_table_path <- file.path(trimmed_dir, "adapter_barcode_table.csv")
        existing_table <- if (file.exists(existing_table_path)) {
          tryCatch(fread(existing_table_path), error = function(e) NULL)
        } else NULL
        trim_marker_values <- stats::setNames(
          lapply(names(trim_marker_values), function(run) {
            v <- trim_marker_values[[run]]
            if (is.data.frame(v)) return(v)
            if (!is.null(existing_table) && run %in% existing_table$id)
              return(existing_table[id == run])
            data.table(id = run)
          }), names(trim_marker_values))
        barcodes_dt <- rbindlist(trim_marker_values, fill = TRUE)
        fwrite(barcodes_dt, file.path(trimmed_dir, "adapter_barcode_table.csv"))

        # Record-only P-shift-failure-cause signal (see
        # R/pshift_diagnostics.R): too-few-raw-reads is knowable as soon
        # as trimming finishes, well before align/pshift ever run.
        check_too_few_reads(conf["bam"], trimmed_dir, config$min_raw_reads_pshift)

        set_flag(config, "trim", conf["exp"])
        if (config$delete_raw_files && length(all_files) > 0) fs::file_delete(unlist(all_files))
    }
}



#' Try to use SSD scratch space for one sample's STAR alignment input +
#' output; fall back to the main drive on any space/copy problem
#'
#' Purely opportunistic -- see \code{config$ssd_scratch_dir}'s own
#' roxygen (\code{\link{pipeline_config}}) for the full rationale and
#' the benchmark this is based on. Never raises: any failure here
#' (insufficient space, or the copy itself erroring, e.g. another
#' process filled the SSD in the gap between the space check and the
#' copy) just means this one sample runs on the main drive instead,
#' exactly like it would if \code{ssd_scratch_dir} were NULL.
#' @param file1 character, the (already main-drive) trimmed fastq path
#' @param output_dir character, this organism's real, main-drive bam dir
#' @param run character, this sample's run id (used to scope the
#' scratch subdirectory and the copy-back glob so a problem with one
#' sample's scratch copy can never touch a sibling's files)
#' @param config the mNGSp config object
#' @return list(output.dir, file1, used_ssd, cleanup = function())
#' @noRd
resolve_align_scratch <- function(file1, output_dir, run, config) {
  main_drive_result <- list(output.dir = output_dir, file1 = file1,
                            used_ssd = FALSE, cleanup = function() invisible(NULL))
  if (is.null(config$ssd_scratch_dir)) return(main_drive_result)

  needed_gb <- 2 * (file.info(file1)$size / 1024^3) + config$ssd_min_free_gb
  free_gb <- tryCatch({
    drive <- ORFik::detect_drive(config$ssd_scratch_dir)
    as.numeric(sub("[a-zA-Z]+", "", ORFik::get_system_usage(drive = drive)$Drive_Free))
  }, error = function(e) NA_real_)
  if (is.na(free_gb) || free_gb < needed_gb) return(main_drive_result)

  scratch_dir <- file.path(config$ssd_scratch_dir, basename(output_dir), run)
  scratch_input <- file.path(scratch_dir, basename(file1))
  copied <- tryCatch({
    dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
    fs::file_copy(file1, scratch_input)
    TRUE
  }, error = function(e) {
    # Cleanup BEFORE the warning, not after: a caller further up the
    # stack that happens to catch warnings via tryCatch() (unlike
    # try(), which only intercepts errors) would otherwise unwind past
    # this handler the moment warning() is called, skipping unlink()
    # and leaving a partial scratch dir behind.
    unlink(scratch_dir, recursive = TRUE)
    warning("SSD scratch copy failed for ", run, ", falling back to main drive: ",
           conditionMessage(e), call. = FALSE)
    FALSE
  })
  if (!copied) return(main_drive_result)

  list(output.dir = scratch_dir, file1 = scratch_input, used_ssd = TRUE,
      cleanup = function() unlink(scratch_dir, recursive = TRUE))
}

#' Remove contaminants and align the reads to genome, for ONE organism of
#' ONE study. Single and paired end reads are handled separately, so the
#' resulting logs are renamed with a "_SINGLE" and/or "_PAIRED" suffix.
#'
#' The per-pair alignment loop plus alignment_final_checks() run in a
#' subprocess (STAR/multiQC's own console output does not otherwise get
#' captured, see run_experiment_subprocess()) that this function polls for
#' progress, rendering the checklist on each tick.
#'
#' STAR's shared-memory genome loading (keep.index.in.memory) is claimed
#' (\code{\link{star_index_claim}}) for the duration of this call and
#' released on exit regardless of success/failure, so a DIFFERENT
#' process aligning a different study of the SAME organism concurrently
#' (a separate run_pipeline() -- online vs local, Ribo-seq vs RNA-seq --
#' or an ad hoc manual script) can never have its still-in-use index
#' removed out from under it. This replaces an earlier, incorrect
#' assumption ("already fully loaded and unloaded within one
#' experiment's loop today, so wrapping at this per-experiment
#' granularity introduces no cross-process shared-memory risk") that a
#' real ~70-minute silent stall this session (PRJNA1051101-homo_sapiens,
#' 2026-10-08) directly contradicted: removing a shared-memory segment
#' out from under a live STAR process does not fail loudly, it hangs.
#'
#' Each sample's alignment input/output is routed through
#' \code{\link{resolve_align_scratch}} first (opportunistic SSD scratch
#' copy, transparently falls back to the main drive otherwise).
#' @inheritParams pipeline_download
#' @param organism character, one name from \code{names(pipeline$organisms)}
#' @param pipelines the full pipelines list, used only to render the
#' checklist; defaults to just this pipeline if called standalone.
#' @param keep_loaded_after logical, default FALSE. Set TRUE when the
#' caller already knows (from a sorted, flattened dispatch across many
#' studies, see \code{\link{pipe_align_clean}}) that another
#' (study, organism) unit for this SAME organism is coming up next --
#' skips even trying to remove the index, purely a scheduling
#' optimization layered on top of (never a substitute for) the
#' cross-process claim check above.
#' @noRd
pipeline_align_one_organism <- function(pipeline, organism, config, pipelines = list(pipeline),
                                        keep_loaded_after = FALSE) {
  # TODO: fix multiqc error for trimmed
  study <- pipeline$study
  did_contamint_removal <- "contam" %in% names(config$flag)
  did_collapse <- "collapsed" %in% names(config$flag)
  did_trim <- "trim" %in% names(config$flag)
  can_use_raw <- did_trim & !did_collapse & !config$delete_raw_files
  steps <- if(can_use_raw) {"tr"} else NULL
  steps <- c(steps, if(did_contamint_removal) {"co"} else NULL)
  steps <- paste(c(steps, "ge"), collapse = "-")

  star.path <- STAR.install()
  fastp.path <- install.fastp()

  conf <- pipeline$organisms[[organism]]$conf
  if (!step_is_next_not_done(config, "aligned", conf["exp"])) return(invisible(NULL))
  index <- pipeline$organisms[[organism]]$index
  runs_full <- study[ScientificName == organism]
  experiment <- conf["exp"]
  # browser()
  trimmed_dir <- fs::path(conf["bam"], "trim")
  raw_fastq_dir <- conf["fastq"]
  output_dir <- conf["bam"]

  keep.unaligned.genome <- config$keep.unaligned.genome

  # Alignment
  input_dir <- if (did_collapse) {
    c(fs::path(trimmed_dir, "SINGLE"), fs::path(trimmed_dir, "PAIRED"))
  } else ifelse(did_trim, trimmed_dir, raw_fastq_dir)
  input_dir <- input_dir[dir.exists(input_dir)]

  # Resume support: filter to the not-yet-done subset BEFORE resolving
  # files (not after): run_files_organizer() requires every row it's
  # given to resolve to a real file, but an already-aligned sibling's
  # trimmed/collapsed input can be gone (config$delete_trimmed_files/
  # delete_collapsed_files) -- resolving against the full study here
  # would error on that missing file even when only reprocessing one
  # already-done study's single fixed sample. Same class of bug as
  # pipeline_trim()/pipeline_collapse()'s own run_files_organizer()
  # calls.
  done_runs <- samples_done(config, "aligned", experiment)
  runs <- runs_full[!(Run %in% done_runs)]
  pairs <- if (nrow(runs) > 0) run_files_organizer(runs, input_dir) else list()

  first_pair_index <- which(lengths(pairs) == 2)[1]
  any_paired <- !is.na(first_pair_index)
  strandMode <- 1
  if (any_paired) {
    if (!file.exists(file.path(output_dir, "strandMode.rds"))) {
      message("Paired end ignored for now, detecting forward direction read file and
          using that only!")
      genomeDir <- file.path(index, "genomeDir")
      strandMode <- detect_strand_mode(pairs[[first_pair_index]][1],
                                       pairs[[first_pair_index]][2],
                                       genomeDir, star = star.path,
                                       keepGenomeLoaded = "LoadAndKeep")

      saveRDS(strandMode, file.path(output_dir, "strandMode.rds"))
    } else strandMode <- readRDS(file.path(output_dir, "strandMode.rds"))
  }

  if (nrow(runs) > 0) {
    owner_id <- star_index_claim(index)
    on.exit(star_index_release(index, owner_id), add = TRUE)
    run_experiment_subprocess(
      func = function(pairs, runs, strandMode, output_dir, index, steps,
                      keep.unaligned.genome, star.path, fastp.path,
                      input_dir, config, experiment, owner_id, keep_loaded_after) {
        cat("Total number of files are:\n")
        cat(length(pairs)); cat("\n")
        for (pair_index in seq_along(pairs)) {
          R1_R2 <- pairs[[pair_index]]
          cat("Single end mode\n")
          cat("Run ", pair_index, " / ", length(pairs), "\n")
          single_end <- lengths(pairs[pair_index]) == 1
          run <- runs[pair_index]$Run
          is_last_sample <- identical(R1_R2, tail(pairs, 1)[[1]])
          # Keep the shared-memory index loaded unless this is the last
          # sample of THIS call AND nothing else (neither a known
          # upcoming same-organism unit, nor -- checked fresh right
          # here, right before the only point that matters -- a
          # different process's still-active claim) needs it anymore.
          keep.index.in.memory <- !is_last_sample || keep_loaded_after ||
            !star_index_safe_to_remove(index, owner_id)
          file1 <- R1_R2[ifelse(single_end, 1, strandMode)]
          file2 <- if (!is.na(R1_R2[2]) & FALSE) {R1_R2[ifelse(strandMode == 1, 2, 1)]}

          scratch <- resolve_align_scratch(file1, output_dir, run, config)
          align_ok <- tryCatch({
            ORFik::STAR.align.single(scratch$file1, file2,
                                     output.dir = scratch$output.dir,
                                     index.dir = index, steps = steps,
                                     resume = "ge",
                                     keep.index.in.memory = keep.index.in.memory,
                                     keep.unaligned.genome = keep.unaligned.genome,
                                     star.path = star.path, fastp = fastp.path
            )
            TRUE
          }, error = function(e) {
            if (!scratch$used_ssd) stop(e) # a real alignment failure, not an SSD problem -- propagate
            scratch$cleanup() # before the warning -- see resolve_align_scratch()'s own comment on why
            warning("STAR run on SSD scratch failed for ", run, ", retrying on main drive: ",
                   conditionMessage(e), call. = FALSE)
            FALSE
          })
          if (identical(align_ok, FALSE)) {
            # Retry once, directly on the main drive -- covers the race
            # where free space looked fine at resolve_align_scratch()'s
            # own check but something else filled the SSD during the
            # copy/alignment itself.
            ORFik::STAR.align.single(file1, file2,
                                     output.dir = output_dir,
                                     index.dir = index, steps = steps,
                                     resume = "ge",
                                     keep.index.in.memory = keep.index.in.memory,
                                     keep.unaligned.genome = keep.unaligned.genome,
                                     star.path = star.path, fastp = fastp.path
            )
          } else if (scratch$used_ssd) {
            # Copy only THIS sample's own output back, by run-id match --
            # same "never touch a sibling's files" precision already used
            # by pipeline_cleanup()'s has_pending_native_output matching,
            # never a whole-directory move.
            produced <- list.files(file.path(scratch$output.dir, "aligned"),
                                   pattern = run, full.names = TRUE)
            dir.create(file.path(output_dir, "aligned"), recursive = TRUE, showWarnings = FALSE)
            fs::file_copy(produced, file.path(output_dir, "aligned", basename(produced)),
                         overwrite = TRUE)
            scratch$cleanup()
          }
          #TODO: Now _1 and _2 will be kept, but fixed in cleanup, do I want it like that ?
          set_sample_flag(config, "aligned", experiment, run)
        }
        alignment_final_checks(input_dir, output_dir, runs, pairs, config, steps)
      },
      args = list(pairs = pairs, runs = runs, strandMode = strandMode,
                  output_dir = output_dir, index = index, steps = steps,
                  keep.unaligned.genome = keep.unaligned.genome,
                  star.path = star.path, fastp.path = fastp.path,
                  input_dir = input_dir, config = config, experiment = experiment,
                  owner_id = owner_id, keep_loaded_after = keep_loaded_after),
      logfile_out = file.path(pipeline_log_base(config), "console", "align", paste0(experiment, ".out.log")),
      logfile_err = file.path(pipeline_log_base(config), "console", "align", paste0(experiment, ".err.log")),
      on_poll = function() pipeline_checklist(pipelines, config)
    )
  }
  set_flag(config, "aligned", conf["exp"])
}

#' Remove contaminants and align the reads to genome, for every organism
#' of one study, in whatever order \code{names(pipeline$organisms)}
#' gives. Thin per-study wrapper around
#' \code{\link{pipeline_align_one_organism}} for callers that want
#' whole-study semantics (e.g. standalone/interactive use); the live
#' pipeline's own dispatch (\code{\link{pipe_align_clean}}) calls
#' \code{pipeline_align_one_organism()} directly instead, in a
#' cross-study, organism-sorted order.
#' @inheritParams pipeline_download
#' @param pipelines the full pipelines list, used only to render the
#' checklist; defaults to just this pipeline if called standalone.
pipeline_align <- function(pipeline, config, pipelines = list(pipeline)) {
  for (organism in names(pipeline$organisms)) {
    pipeline_align_one_organism(pipeline, organism, config, pipelines)
  }
}

#' Post-alignment checks/cleanup for one pipeline_align() call
#' @param pairs the resolved input files (one list element per sample
#' just processed by this call) -- used ONLY to know exactly which
#' input files are safe to delete below, never the whole directory. A
#' sibling's own input file, already deleted after ITS OWN earlier,
#' separate pipeline_align() call, must never be touched again here.
#' @noRd
alignment_final_checks <- function(input_dir, output_dir, runs, pairs, config, steps) {
  cleanup_script <- system.file("STAR_Aligner", "cleanup_folders.sh",
                                package = "ORFik")
  system2("/bin/bash", c(cleanup_script, output_dir))
  STAR.allsteps.multiQC(output_dir, steps = steps)

  # Record-only P-shift-failure-cause signal (see R/pshift_diagnostics.R):
  # alignment rate is knowable right after align, well before pshift runs.
  check_alignment_rate(output_dir, config$min_alignment_rate_pshift)

  dir_info <- as.data.table(fs::dir_info(file.path(output_dir, "aligned"), type = "file"))
  dir_info <- dir_info[grep(paste(runs$Run, collapse = "|"), path)][grep("\\.bam$", path)]
  if (nrow(dir_info) < nrow(runs)) {
    stop("You have missing bam files in your aligned folder!")
  }
  if (any(dir_info$size == 0)) {
    stop("You have empty bam files in your aligned folder!")
  }
  if (config$delete_collapsed_files) {
    # Only the input files actually resolved/used for the samples JUST
    # processed here -- never the whole directory (fs::dir_delete(input_dir)
    # used to do exactly that), which would also destroy every sibling's
    # still-needed collapsed/trimmed input the moment any single sample
    # gets reprocessed on its own. Confirmed live this is a real,
    # reachable path, not just a theoretical risk -- see the identical
    # class of bug already fixed in pipeline_trim()/pipeline_collapse().
    used_files <- unlist(pairs, use.names = FALSE)
    used_files <- used_files[!is.na(used_files) & file.exists(used_files)]
    if (length(used_files) > 0) fs::file_delete(used_files)
  }
}

#' Remove contaminants and align the reads to genome. Single and paired end
#' reads are handled separately, so the resulting logs are renamed with a
#' "_SINGLE" and/or "_PAIRED" suffix.
#'
#' Runs in a subprocess for log hygiene (STAR's own console output does
#' not otherwise get captured, see run_experiment_subprocess()). This is a
#' whole-folder STAR call (ORFik::STAR.align.folder()), not a per-sample R
#' loop, so unlike pipeline_trim()/pipeline_align() there is no natural
#' per-sample point to hook a marker write into -- the checklist only ever
#' shows study-level counts for this stage, same limitation as
#' pshift/valid_pshift/pcounts (ORFik's own internal BiocParallel dispatch).
#'
#' Unlike \code{\link{pipeline_align_one_organism}}, this never passes
#' \code{keep.index.in.memory}/\code{keepGenomeLoaded} to
#' \code{ORFik::STAR.align.folder()} at all, so it always defaults to
#' FALSE -- the contaminant index is never put into shared memory here,
#' so it is not exposed to the cross-process race that function guards
#' against, and needs no claim/release of its own. Still extracted to
#' one-organism granularity (same as align/cleanup) purely so
#' \code{\link{pipe_align_clean}}'s flattened, organism-sorted dispatch
#' can call all three stages at matching granularity -- a study's
#' mouse-organism contaminant removal must never block its
#' already-ready human-organism alignment from starting.
#' @inheritParams pipeline_download
#' @param organism character, one name from \code{names(pipeline$organisms)}
#' @param pipelines the full pipelines list, used only to render the
#' checklist; defaults to just this pipeline if called standalone.
#' @noRd
pipeline_align_contaminants_one_organism <- function(pipeline, organism, config, pipelines = list(pipeline)) {
    study <- pipeline$study
    conf <- pipeline$organisms[[organism]]$conf
    if (!step_is_next_not_done(config, "contam", conf["exp"])) return(invisible(NULL))
    index <- pipeline$organisms[[organism]]$index
    runs <- study[ScientificName == organism]
    trimmed_dir <- fs::path(conf["bam"], "trim")
    output_dir <- conf["bam"]
    did_collapse <- "collapsed" %in% names(config$flag)
    keep.contaminants <- config$keep_contaminants
    experiment <- conf["exp"]
    run_experiment_subprocess(
      func = function(runs, trimmed_dir, output_dir, did_collapse, keep.contaminants, index) {
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
      },
      args = list(runs = runs, trimmed_dir = trimmed_dir, output_dir = output_dir,
                  did_collapse = did_collapse, keep.contaminants = keep.contaminants, index = index),
      logfile_out = file.path(pipeline_log_base(config), "console", "contam", paste0(experiment, ".out.log")),
      logfile_err = file.path(pipeline_log_base(config), "console", "contam", paste0(experiment, ".err.log")),
      on_poll = function() pipeline_checklist(pipelines, config)
    )
    set_flag(config, "contam", conf["exp"])
}

#' Remove contaminants for every organism of one study, in whatever
#' order \code{names(pipeline$organisms)} gives. Thin per-study wrapper
#' around \code{\link{pipeline_align_contaminants_one_organism}} for
#' callers that want whole-study semantics; \code{\link{pipe_align_clean}}
#' calls the one-organism version directly instead.
#' @inheritParams pipeline_download
#' @param pipelines the full pipelines list, used only to render the
#' checklist; defaults to just this pipeline if called standalone.
pipeline_align_contaminants <- function(pipeline, config, pipelines = list(pipeline)) {
  for (organism in names(pipeline$organisms)) {
    pipeline_align_contaminants_one_organism(pipeline, organism, config, pipelines)
  }
}

#' Remove all files apart from logs and final aligned BAMs.
#' Rename the BAMs into <run_accession>.bam format. As a final step, for
#' presets that pre-alignment-collapse reads (e.g. Ribo-seq's "collapsed"
#' flag), save the true (uncollapsed) per-read alignment metrics table +
#' plot -- see \code{\link{save_expanded_alignment_metrics}} for why the
#' native STAR/collapsed-level report can be badly misleading.
#' @inheritParams pipeline_download
#' @param organism character, one name from \code{names(pipeline$organisms)}
#' @noRd
pipeline_cleanup_one_organism <- function(pipeline, organism, config) {
    accession <- pipeline$accession
    study <- pipeline$study
    did_contamint_removal <- "contam" %in% names(config$flag)
    did_collapse <- "collapsed" %in% names(config$flag)
        conf <- pipeline$organisms[[organism]]$conf
        if (!step_is_next_not_done(config, "cleanbam", conf["exp"])) return(invisible(NULL))
        study_org <- study[ScientificName == organism,]
        bam_dir <- fs::path(conf["bam"], "aligned")
        if (did_contamint_removal) {
          fs::file_delete(fs::dir_ls(
            fs::path(conf["bam"], "contaminants_depletion"),
            glob = "**/*.out.*"
          ))
        }
        # Resume support: a sample needs renaming here only if a
        # pending, not-yet-renamed STAR-native output actually exists
        # for it -- NOT simply "does <Run>.bam not exist yet". A sample
        # being freshly reprocessed can have BOTH a stale <Run>.bam
        # left over from its previous run AND a fresh pending native
        # output waiting to replace it; checking only "does the target
        # already exist" (an earlier version of this fix) wrongly
        # skipped that sample, leaving the stale old BAM in place and
        # the fresh realignment never promoted. Confirmed live,
        # PRJNA770650-homo_sapiens, 2026-10-06. Checking for a pending
        # native output also still skips a truly already-done sibling
        # (no such file exists for it), avoiding the original bug this
        # guarded against: match_bam_to_metadata()/run_files_organizer()'s
        # loose substring match otherwise matching an already-renamed
        # file and fs::file_move() moving it onto itself, which
        # reports [ENOENT] "no such file" and DELETES it instead of
        # no-op'ing (PRJNA926112-homo_sapiens, 2026-10-05).
        has_pending_native_output <- vapply(study_org$Run, function(run) {
          length(list.files(bam_dir, pattern = paste0("^.*", run, ".*_Aligned\\.sortedByCoord\\.out\\.bam$"))) > 0
        }, logical(1))
        rename_study_org <- study_org[has_pending_native_output]
        new_file_names <- fs::path(bam_dir, rename_study_org$Run, ext = "bam")

        if (nrow(rename_study_org) > 0) {
          old_file_names <- match_bam_to_metadata(bam_dir, rename_study_org, FALSE,
                                                  format = c("_Aligned.sortedByCoord.out.bam"))

          file_names_to_delete <- new_file_names[new_file_names != old_file_names]
          file_names_to_delete <- file_names_to_delete[file.exists(file_names_to_delete)]
          if (length(file_names_to_delete) > 0) try(file.remove(file_names_to_delete), silent = TRUE)

          stopifnot(length(old_file_names) == length(new_file_names))
          for (i in seq_along(old_file_names)) {
            fs::file_move(old_file_names[i], new_file_names[i])
          }
        }

        if (did_collapse) {
          collapsed_dir <- fs::path(conf["bam"], "trim", "SINGLE")
          save_expanded_alignment_metrics(bam_dir, collapsed_dir, study_org,
                                          BPPARAM = bpparam_from_config(config, "trim"))
        }

        set_flag(config, "cleanbam", conf["exp"])
}

#' Remove all files apart from logs and final aligned BAMs, for every
#' organism of one study, in whatever order
#' \code{names(pipeline$organisms)} gives. Thin per-study wrapper around
#' \code{\link{pipeline_cleanup_one_organism}} for callers that want
#' whole-study semantics; \code{\link{pipe_align_clean}} calls the
#' one-organism version directly instead, immediately once that SAME
#' organism's alignment succeeds -- never waiting on sibling organisms
#' of the same study, which already proceed fully independently (every
#' later stage's own flag/resume state and I/O directories are keyed by
#' the per-organism experiment name, not the whole-study accession).
#' @inheritParams pipeline_download
pipeline_cleanup <- function(pipeline, config) {
    for (organism in names(pipeline$organisms)) {
        pipeline_cleanup_one_organism(pipeline, organism, config)
    }
}

#' Create ORFik experiment from config study and bam files
#' @inheritParams pipeline_download
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
        # ORFik::create.experiment()'s own `if (author != "")` check
        # errors ("missing value where TRUE/FALSE needed") on a literal
        # NA, not just an empty string -- confirmed live,
        # PRJEB50305-saccharomyces_cerevisiae, 2026-10-06: AUTHOR was
        # genuinely NA for every sample in that study's metadata at run
        # time (a real metadata gap, not something this call can fix),
        # and unique(study$AUTHOR) passed that NA straight through.
        # Collapse to the safe, already-expected "no author" value
        # instead, and tolerate more than one distinct author too
        # (create.experiment() expects a single scalar).
        author_value <- unique(study$AUTHOR[!is.na(study$AUTHOR) & study$AUTHOR != ""])
        author_value <- if (length(author_value) == 0) "" else paste(author_value, collapse = "; ")
        ORFik::create.experiment(
            dir = bam_dir,
            exper = experiment, txdb = paste0(annotation["gtf"], ".db"),
            libtype = study$LIBRARYTYPE,
            fa = annotation["genome"], organism = organism,
            stage = metadata_clean$stage, rep = study$REPLICATE,
            condition = metadata_clean$condition,
            fraction = metadata_clean$fraction,
            pairedEndBam = metadata_clean$paired_end,
            author = author_value,
            files = bam_files, runIDs = study$Run
        )
        df <- ORFik::read.experiment(experiment,
                                     output.env = new.env())
        set_flag(config, "exp", conf["exp"])
        df_list <- c(df_list, df)
    }
  return(df_list)
}

#' Run one format conversion per sample, with per-sample resume support
#'
#' Subsets \code{df} to one row at a time and calls \code{convert_fun} on
#' that single-row experiment instead of the whole experiment at once --
#' ORFik's own \code{convert_bam_to_ofst()}/\code{convert_to_covRleList()}/
#' \code{convert_to_bigWig()} are already plain serial per-file
#' \code{for} loops internally (verified directly: no \code{BPPARAM}, no
#' cross-sample state besides directory creation and \code{seqinfo()}),
#' so calling them one row at a time is functionally identical to calling
#' them once on the whole \code{df} -- just resumable: a sample already
#' marked done in an earlier, interrupted attempt is skipped instead of
#' redone, the same resume pattern \code{pipeline_trim()}/
#' \code{pipeline_align()} already use (see \code{R/pipeline_sample_flags.R}).
#'
#' Deliberately NOT extended to pshift/pcounts: pshift hands off to
#' ORFik's own internal per-experiment shifting (no natural per-sample
#' point to hook into without changing ORFik itself), and pcounts'
#' \code{countTable_regions()} needs every sample's data together to
#' build one \code{SummarizedExperiment} -- neither has a safe per-sample
#' subset-and-resume equivalent the way a plain per-file conversion loop
#' does.
#'
#' Dispatches across not-yet-done samples via \code{BPPARAM} rather than
#' a serial \code{for} loop -- each \code{convert_fun()} call is
#' independent (its own file in, its own file out), so there was no
#' reason this needed to be serial. Confirmed live, GSE151959-homo_sapiens,
#' 2026-10-07: bigwig conversion alone took ~7.5 minutes for 25 samples
#' run one at a time (~15-18s each); this and the other two conversions
#' (ofst, covrle) all shared the exact same unnecessary serial pattern.
#' \code{config$threads} has no per-step entry for "ofst"/"covrle"/
#' "bigwig" (only trim/collapse/pshifted/valid_pshift/pcounts are
#' defined there) -- uses \code{config$threads$default} (the full
#' configured worker count), same as any other step without its own
#' dedicated, narrower thread budget.
#' @param df an ORFik experiment, already at the correct
#' \code{uniqueMappers()} setting for this pass (see \code{step_id})
#' @param config the mNGSp config object
#' @param step_id character, the sample-flag namespace for this pass
#' (e.g. \code{"ofst"} or \code{"ofst_unique"} -- kept distinct per
#' mapper-mode so a resume mid-way through the unique-mappers pass never
#' re-skips or re-redoes the separately-tracked all-mappers pass)
#' @param convert_fun function(df_one_row), one of ORFik's own
#' \code{convert_bam_to_ofst()}/\code{convert_to_covRleList()}/
#' \code{convert_to_bigWig()} (the last one wrapped to also supply its
#' \code{in_files} argument per-row)
#' @param BPPARAM BiocParallel param, default \code{bpparam_from_config(config, "default")}
#' @return invisible(NULL)
#' @noRd
convert_per_sample <- function(df, config, step_id, convert_fun,
                               BPPARAM = bpparam_from_config(config, "default")) {
  experiment <- name(df)
  run_ids <- runIDs(df)
  done <- samples_done(config, step_id, experiment)
  idx <- which(!(run_ids %in% done))
  BiocParallel::bplapply(idx, function(i, df, convert_fun, config, step_id, experiment, run_ids) {
    convert_fun(df[i, ])
    set_sample_flag(config, step_id, experiment, run_ids[i])
  }, df = df, convert_fun = convert_fun, config = config, step_id = step_id,
     experiment = experiment, run_ids = run_ids, BPPARAM = BPPARAM)
  invisible(NULL)
}

#' convert_bam_to_ofst(), plus save that one sample's read-length
#' distribution right after
#'
#' Reads the just-written ofst file once (cheap: a compact, local,
#' already-converted file, not a BAM reparse) to save
#' \code{<ofst_dir>/read_length_distribution/<Run>.csv} alongside it --
#' see \code{\link{read_length_distribution}} for why this must be
#' score-weighted, not a plain row count. Wrapped in its own \code{try()}:
#' a bug in this diagnostic-only addition must never fail the real ofst
#' conversion it rides along with (same reasoning already used for
#' \code{\link{save_expanded_alignment_metrics}}).
#' @param df_one_row an ORFik experiment subset to exactly one row (see
#' \code{\link{convert_per_sample}})
#' @return invisible(NULL)
#' @noRd
convert_bam_to_ofst_with_length_dist <- function(df_one_row) {
  convert_bam_to_ofst(df_one_row)
  try({
    ofst_path <- filepath(df_one_row, "ofst")
    out_dir <- file.path(dirname(ofst_path), "read_length_distribution")
    save_read_length_distribution(ofst_path, out_dir, runIDs(df_one_row))
  }, silent = TRUE)
  invisible(NULL)
}

#' Create ORFik ofst files from bam files
#' @param df_list a list of ORFik experiments
#' @inheritParams pipeline_download
pipeline_create_ofst <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "ofst", name(df))) next
    if (config$all_mappers) {
      convert_per_sample(df, config, "ofst", convert_bam_to_ofst_with_length_dist)
      try(aggregate_read_length_distribution(file.path(dirname(filepath(df, "ofst")[1]),
                                                        "read_length_distribution")), silent = TRUE)
    }

    if (config$split_unique_mappers) {
      df_unique <- df
      uniqueMappers(df_unique) <- TRUE
      convert_per_sample(df_unique, config, "ofst_unique", convert_bam_to_ofst_with_length_dist)
      try(aggregate_read_length_distribution(file.path(dirname(filepath(df_unique, "ofst")[1]),
                                                        "read_length_distribution")), silent = TRUE)
    }
    set_flag(config, "ofst", name(df))
  }
}

#' Save per-sample + aggregated read-length distribution for every
#' already-shifted ofst in this experiment
#'
#' No per-sample loop exists in \code{pipeline_pshift()} itself
#' (\code{shiftFootprintsByExperimentSafe()} shifts the whole experiment
#' via ORFik's own internal dispatch -- same constraint already accepted
#' for the ofst/covRLE/bigwig per-sample-resume work this session, which
#' deliberately did not extend to pshift for the same reason). Called
#' once, right after a successful shift, looping over the shifted files
#' that already exist on disk at that point (cheap: no re-shift, just a
#' read of each small ofst). Comparing this against the pre-shift
#' distribution (see \code{convert_bam_to_ofst_with_length_dist()})
#' directly shows how much got dropped by the \code{accepted_lengths}
#' filter. Wrapped in its own \code{try()}, same reasoning as the ofst
#' side: a bug here must never fail the real shift it rides along with.
#' @param df an ORFik experiment, already at the correct
#' \code{uniqueMappers()} setting for this pass (shifted file paths come
#' from \code{filepath(df, "pshifted")})
#' @return invisible(NULL)
#' @noRd
save_pshifted_length_distributions <- function(df) {
  try({
    pshifted_paths <- filepath(df, "pshifted")
    run_ids <- runIDs(df)
    out_dir <- file.path(dirname(pshifted_paths[1]), "read_length_distribution")
    for (i in seq_along(pshifted_paths)) {
      if (file.exists(pshifted_paths[i])) {
        save_read_length_distribution(pshifted_paths[i], out_dir, run_ids[i])
      }
    }
    aggregate_read_length_distribution(out_dir)
  }, silent = TRUE)
  invisible(NULL)
}

#' Thin wrapper around \code{ORFik::filepath(df_one_row, "ofst")}
#'
#' Same rationale as \code{\link{pshifted_filepath}} (R/shift_qc_cache.R):
#' exists purely so tests can mock ONE plain massiveNGSpipe function
#' instead of needing \code{filepath()} itself to dispatch on a fake
#' experiment class -- \code{filepath()} is a plain function in ORFik
#' (not an S4 generic) with an internal \code{stopifnot(is(df, "experiment"))},
#' so it can never be given a method for a test stub class.
#' @param df an ORFik experiment subset, one or more rows
#' @return character vector, one path per row of \code{df}
#' @noRd
ofst_filepath <- function(df) filepath(df, "ofst")

#' Is this one sample's pshifted file already valid (present and at
#' least as new as its own raw ofst file)?
#'
#' The freshness check (not just existence) is what makes a redone
#' sample (e.g. a barcode fix followed by a reshift) correctly trigger
#' recomputation -- its ofst file gets a newer mtime than the old
#' pshifted file, which was shifted from the PREVIOUS version of that
#' ofst. Same pattern as \code{\link{shift_qc_cache_valid}}.
#' @param df_one_row an ORFik experiment subset to exactly one row,
#' already at the correct \code{uniqueMappers()} setting for this pass
#' @return logical
#' @noRd
pshift_sample_valid <- function(df_one_row) {
  ofst_path <- ofst_filepath(df_one_row)
  if (!file.exists(ofst_path)) return(FALSE)
  pshifted_path <- pshifted_filepath(df_one_row)
  if (!file.exists(pshifted_path)) return(FALSE)
  file.info(pshifted_path)$mtime >= file.info(ofst_path)$mtime
}

#' One mapper-mode pass of pshift, subsetting to only the samples that
#' actually need it
#'
#' \code{ORFik::shiftFootprintsByExperiment()} (called inside
#' \code{shiftFootprintsByExperimentSafe()}) reloads, reshifts, and
#' rewrites EVERY sample it's given, even the ones whose shift it
#' reuses from an existing \code{shift.list} entry instead of
#' re-detecting -- re-detection is skipped per-sample, but the
#' load/apply/write cost is not. Subsetting \code{df} to only the
#' not-yet-valid samples before calling it skips that cost entirely for
#' everyone else, with no scratch-dir/copy-back needed: confirmed
#' directly from ORFik's own source that this is already safe --
#' \code{libFolder()} (hence \code{out.dir}) resolves from
#' \code{dirname(x$filepath[1])}, identical for a subset or the full
#' experiment, and \code{shifts_save()} already merges a subset's
#' shifts into an existing \code{shifting_table.rds} by name rather
#' than overwriting it (\code{old_shifts[names(shifts)] <- shifts}) --
#' this is exactly the "pshift a subset, the existing shift file gets
#' the new entries merged in" mechanism, already built into ORFik,
#' not something added here. \code{shift.list = NULL} for the subset
#' call: every sample in it lacks a valid pshifted file by construction,
#' so none has anything worth reusing -- a few extra seconds of fresh
#' auto-detection per sample is a non-issue next to the reload/rewrite
#' cost this whole function exists to avoid.
#'
#' One real caveat, informational only: \code{shiftFootprintsByExperimentSafe()}'s
#' own FFT-strength marker file (\code{strong_periodicity.rds} etc.,
#' \code{\link{fft_strength_files}}) is experiment-level, not
#' per-sample -- calling it on a small subset overwrites that marker
#' with just this call's own outcome, which may not represent the
#' whole experiment's shift history anymore. Doesn't affect the actual
#' shifted data or the shift table, only that one diagnostic file.
#'
#' \code{reuse_shifts_if_existing} is still honored for the narrow edge
#' case where a sample's pshifted FILE is missing/stale but its
#' shift-table ENTRY is still present and usable (e.g. the file was
#' deleted by something else while the table itself wasn't touched):
#' that entry is reused (skips re-detection, matching the pre-subsetting
#' behavior) rather than always forcing a fresh auto-detect for
#' everything in the subset.
#' @param df an ORFik experiment, already at the correct
#' \code{uniqueMappers()} setting for this pass
#' @param accepted_lengths,allowed_hard12_species,BPPARAM,max_no_adapter_removed_pct
#' passed straight through to \code{shiftFootprintsByExperimentSafe()}
#' @param reuse_shifts_if_existing logical
#' @return whatever \code{shiftFootprintsByExperimentSafe()} returns
#' when there was anything to do; \code{1} (matching the no-op default
#' elsewhere in this file) if every sample was already valid
#' @noRd
pshift_needed_subset <- function(df, accepted_lengths, allowed_hard12_species, BPPARAM,
                                 max_no_adapter_removed_pct, reuse_shifts_if_existing) {
  needs <- !vapply(seq_len(nrow(df)), function(i) pshift_sample_valid(df[i, ]), logical(1))
  if (!any(needs)) return(1)
  df_subset <- df[needs, ]

  shifting_table <- NULL
  if (reuse_shifts_if_existing) {
    full_table <- try(shifts_load_safe(df, TRUE), silent = TRUE)
    if (!is(full_table, "try-error") && !is.null(full_table)) {
      subset_paths <- ofst_filepath(df_subset)
      candidate <- full_table[subset_paths]
      # shifts_load_safe() pads any sample the table doesn't cover with
      # a NULL/NA placeholder (so its own length check against nrow(df)
      # still passes) -- that placeholder must fall through to fresh
      # auto-detection inside shiftFootprintsByExperiment()'s own
      # per-file `if (is.null(shifts))` check, never be passed through
      # as if it were a real known shift.
      has_real_entry <- !vapply(candidate, function(x) is.null(x) || (length(x) == 1 && is.na(x)), logical(1))
      if (any(has_real_entry)) shifting_table <- candidate[has_real_entry]
    }
  }

  res <- shiftFootprintsByExperimentSafe(df_subset, shifting_table, accepted_lengths,
                                         allowed_hard12_species, BPPARAM,
                                         max_no_adapter_removed_pct)
  if (!inherits(res, "error")) save_pshifted_length_distributions(df_subset)
  res
}

#' @inheritParams pipeline_create_ofst
pipeline_pshift <- function(df_list, config, accepted_lengths = config$accepted_lengths_rpf,
                            reuse_shifts_if_existing = config$reuse_shifts_if_existing,
                            allowed_hard12_species = c("Escherichia coli", "Sars cov2"),
                            BPPARAM = bpparam_from_config(config, "pshifted")) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "pshifted", name(df))) next
    res <- 1
    if (config$all_mappers) {
      res <- pshift_needed_subset(df, accepted_lengths, allowed_hard12_species, BPPARAM,
                                  config$max_no_adapter_removed_pct, reuse_shifts_if_existing)
    }

    if(config$split_unique_mappers & !inherits(res, "error")) {
      uniqueMappers(df) <- TRUE
      res <- pshift_needed_subset(df, accepted_lengths, allowed_hard12_species, BPPARAM,
                                  config$max_no_adapter_removed_pct, reuse_shifts_if_existing)
    }

    if(!inherits(res, "error")) {
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
#' @inheritParams pipeline_create_ofst
#' @return invisible(NULL)
pipeline_validate_shifts <- function(df_list, config) {
  #
  # TODO: This function can be improved, especially how we detect bad shifts!
  for (df in df_list) {
    if (!step_is_next_not_done(config, "valid_pshift", name(df))) next
    BPPARAM <- bpparam_from_config(config, "valid_pshift")
    if (config$all_mappers)
      shift_qc_cached(df, BPPARAM, config$max_no_adapter_removed_pct)
    if (config$split_unique_mappers) {
      uniqueMappers(df) <- TRUE
      shift_qc_cached(df, BPPARAM, config$max_no_adapter_removed_pct)
    }
    set_flag(config, "valid_pshift", name(df))
  }
}

#' @inheritParams pipeline_create_ofst
pipeline_convert_covRLE <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "covrle", name(df))) next
    if (config$all_mappers)
      convert_per_sample(df, config, "covrle", convert_to_covRleList)
    if (config$split_unique_mappers) {
      df_unique <- df
      uniqueMappers(df_unique) <- TRUE
      convert_per_sample(df_unique, config, "covrle_unique", convert_to_covRleList)
    }
    set_flag(config, "covrle", name(df))
  }
}

#' @inheritParams pipeline_create_ofst
pipeline_convert_bigwig <- function(df_list, config) {
  # convert_to_bigWig()'s in_files must be supplied per-row too (it
  # defaults to "pshifted", but this pipeline builds bigwig from the
  # covRLE step's own output instead -- filepath type "cov" -- so the
  # override has to travel with the per-row wrapper, not be left at the
  # function default).
  convert_bigwig_one <- function(df_one_row) convert_to_bigWig(df_one_row, filepath(df_one_row, "cov"))
  for (df in df_list) {
    if (!step_is_next_not_done(config, "bigwig", name(df))) next
    if (config$all_mappers)
      convert_per_sample(df, config, "bigwig", convert_bigwig_one)
    if (config$split_unique_mappers) {
      df_unique <- df
      uniqueMappers(df_unique) <- TRUE
      convert_per_sample(df_unique, config, "bigwig_unique", convert_bigwig_one)
    }
    set_flag(config, "bigwig", name(df))
  }
}

#' @inheritParams pipeline_create_ofst
pipeline_count_table_psites <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "pcounts", name(df))) next
    BPPARAM <- bpparam_from_config(config, "pcounts")
    if (config$all_mappers)
      ORFik::countTable_regions(df, lib.type = "cov", forceRemake = TRUE,
                                BPPARAM = BPPARAM)
    if (config$split_unique_mappers) {
      uniqueMappers(df) <- TRUE
      if (config$all_mappers) remove.experiments(df)
      ORFik::countTable_regions(df, lib.type = "cov", forceRemake = TRUE,
                                BPPARAM = BPPARAM)
    }
    set_flag(config, "pcounts", name(df))
  }
}
