#' Safe worker count for processing a set of files in parallel, bounded
#' by available system memory
#'
#' Finds the largest number of files whose combined size (as a
#' conservative proxy for in-memory footprint -- the actual peak per
#' file depends on compression ratio and intermediate objects built
#' while processing, neither known here) can be held at once without
#' exceeding currently free memory, leaves a small safety margin, then
#' caps the result at \code{max_workers}. Pessimistic by design: this
#' estimates "no more workers than fit in free memory right now", not
#' an optimal count.
#' @param files character vector, paths to the files about to be
#' processed in parallel (one per worker)
#' @param max_workers numeric, the upper bound a caller would otherwise
#' use (e.g. from \code{config$threads$collapse}) -- never returns more
#' than this
#' @param safety_margin_workers integer, default 2. Subtracted from the
#' memory-derived worker count, so some headroom for other processes/
#' overhead always remains
#' @return integer >= 1
#' @noRd
memory_safe_worker_count <- function(files, max_workers, safety_margin_workers = 2) {
  if (length(files) == 0) return(max(max_workers, 1))
  file_sizes_cumsum_GB <- cumsum(file.size(files) / 1e9)
  system_usage <- get_system_usage()
  free_memory_GB <- system_usage$Memory_Total_GB - system_usage$Memory_Usage_GB
  first_over_budget <- head(which(file_sizes_cumsum_GB > free_memory_GB), 1)
  # Every file's cumulative size stays within budget -- memory is not
  # the limiting factor here, so defer entirely to max_workers.
  if (length(first_over_budget) == 0) return(max(max_workers, 1))
  memory_derived <- first_over_budget - safety_margin_workers
  max(min(max_workers, memory_derived), 1)
}

#' Collapse fastq/fasta
#' Collapse trimmed single-end fasta/fastq reads and move them into "trim/SINGLE"
#' subdirectory. PE reads are moved into "trim/PAIRED" without modification.
#'
#' Per-sample resume support, matching \code{\link{pipeline_trim}}/
#' \code{\link{pipeline_align}}'s pattern: a sample already marked done
#' for this experiment+step (\code{samples_done()}, R/pipeline_sample_flags.R)
#' is skipped even if the experiment-level "collapsed" flag itself is
#' still unset -- e.g. after a single sample's trim output was
#' regenerated (a barcode/adapter fix) and the experiment-level flags
#' from "trim" onward cleared, this redoes only that sample's collapse,
#' not the whole study's.
#'
#' Runs inside \code{\link{run_experiment_subprocess}} (like trim/align),
#' for log capture and checklist polling -- unlike align, the per-sample
#' work itself (\code{\link[ORFik]{collapse.fastq}}) is not internally
#' multi-threaded, so it still needs its own \code{BiocParallel::bplapply()}
#' dispatch *inside* that subprocess, bounded by
#' \code{\link{memory_safe_worker_count}()} (collapsing many large fastq
#' files at once is memory-heavy per worker) and \code{config$threads$collapse}
#' (both the paired and single-end branches now go through the same
#' config-aware, memory-checked path -- previously only the single-end
#' branch checked memory at all, and paired always used a flat 16 workers).
#' @param pipeline a pipeline object, subset of init_pipelines output
#' @param config a pipeline_config object
#' @param pipelines the full pipelines list, used only to render the
#' checklist; defaults to just this pipeline if called standalone.
pipeline_collapse <- function(pipeline, config, pipelines = list(pipeline)) {
  study <- pipeline$study
  for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    if (!step_is_next_not_done(config, "collapsed", conf["exp"])) next
    experiment <- conf["exp"]
    trimmed_dir <- fs::path(conf["bam"], "trim")
    runs_paired <- study[ScientificName == organism &
                           LibraryLayout == "PAIRED"]
    runs_single <- study[ScientificName == organism &
                           LibraryLayout != "PAIRED"]

    done_runs <- samples_done(config, "collapsed", experiment)
    todo_paired <- which(!(runs_paired$Run %in% done_runs))
    todo_single <- which(!(runs_single$Run %in% done_runs))

    all_files <- c()
    outdir_paired <- files_paired <- NULL
    if (nrow(runs_paired) > 0) {
      # Paired-end reads are collapsed using only read 1; read 2 is not
      # separately reverse-complemented/merged back in.
      outdir_paired <- fs::path(trimmed_dir, runs_paired[1]$LibraryLayout)
      files_paired <- run_files_organizer(runs_paired, trimmed_dir)
      all_files <- c(all_files, files_paired)
    }
    outdir_single <- files_single <- NULL
    if (nrow(runs_single) > 0) {
      outdir_single <- fs::path(trimmed_dir, "SINGLE")
      files_single <- run_files_organizer(runs_single, trimmed_dir)
      all_files <- c(all_files, files_single)
    }

    run_experiment_subprocess(
      func = function(runs_paired, files_paired, outdir_paired, todo_paired,
                      runs_single, files_single, outdir_single, todo_single,
                      config, experiment) {
        if (length(todo_paired) > 0) {
          fs::dir_create(outdir_paired)
          read1_paths <- heads(files_paired[todo_paired], 1)
          BPPARAM <- bpparam_from_config(config, "collapse",
            workers = memory_safe_worker_count(unlist(read1_paths), config$threads$collapse))
          BiocParallel::bplapply(seq_along(todo_paired), function(i, read1_paths, outdir_paired,
                                                                  runs_paired, todo_paired,
                                                                  config, experiment) {
            ORFik::collapse.fastq(read1_paths[[i]], outdir_paired, compress = TRUE)
            set_sample_flag(config, "collapsed", experiment, runs_paired$Run[todo_paired[i]])
          }, read1_paths = read1_paths, outdir_paired = outdir_paired, runs_paired = runs_paired,
             todo_paired = todo_paired, config = config, experiment = experiment, BPPARAM = BPPARAM)
        }
        if (length(todo_single) > 0) {
          fs::dir_create(outdir_single)
          files_todo <- files_single[todo_single]
          BPPARAM <- bpparam_from_config(config, "collapse",
            workers = memory_safe_worker_count(unlist(files_todo), config$threads$collapse))
          BiocParallel::bplapply(seq_along(todo_single), function(i, files_todo, outdir_single,
                                                                   runs_single, todo_single,
                                                                   config, experiment) {
            ORFik::collapse.fastq(files_todo[[i]], outdir_single, compress = TRUE)
            set_sample_flag(config, "collapsed", experiment, runs_single$Run[todo_single[i]])
          }, files_todo = files_todo, outdir_single = outdir_single, runs_single = runs_single,
             todo_single = todo_single, config = config, experiment = experiment, BPPARAM = BPPARAM)
        }
        invisible(NULL)
      },
      args = list(runs_paired = runs_paired, files_paired = files_paired,
                  outdir_paired = outdir_paired, todo_paired = todo_paired,
                  runs_single = runs_single, files_single = files_single,
                  outdir_single = outdir_single, todo_single = todo_single,
                  config = config, experiment = experiment),
      logfile_out = file.path(pipeline_log_base(config), "console", "collapse", paste0(experiment, ".out.log")),
      logfile_err = file.path(pipeline_log_base(config), "console", "collapse", paste0(experiment, ".err.log")),
      on_poll = function() pipeline_checklist(pipelines, config)
    )

    if (config$delete_trimmed_files) {
      existing <- unlist(all_files, use.names = FALSE, recursive = TRUE)
      existing <- existing[file.exists(existing)]
      if (length(existing) > 0) file.remove(existing)
    }
    set_flag(config, "collapsed", experiment)
  }
}

pipe_cigar_collapse_single <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "cigar_collapse", name(df))) next
    if (config$all_mappers) {
      files <- filepath(df, "ofst")
      for (f in files) {
        ofst <- fimport(f)
        ORFik::export.ofst(collapseDuplicatedReads(ofst), f)
      }
    }

    if (config$split_unique_mappers) {
      uniqueMappers(df) <- TRUE
      files <- filepath(df, "ofst")
      for (f in files) {
        ofst <- fimport(f)
        ORFik::export.ofst(collapseDuplicatedReads(ofst), f)
      }
    }

    set_flag(config, "cigar_collapse", name(df))
  }
}

pipe_cigar_collapse <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "cigar_collapse", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipe_cigar_collapse_single(df_list, config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, cigar_collapse, study: ", experiments[1])
  }
}
