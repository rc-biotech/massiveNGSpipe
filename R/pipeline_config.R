#' Initial config setup
#'
#' Set up all paths, functions to be run and flag directories.
#' Also adds parallel processing settings and google integration.
#' @inheritParams libtype_flags
#' @param project_dir where will specific pipeline outputs be put. Default:
#'  file.path(dirname(config)[1], "NGS_pipeline").
#'  If you need seperat pipelines, please use different locations.
#' @param config path, default \code{ORFik::config()}, where will
#' fastq, bam, references and ORFik experiments go
#' @param complete_metadata path, default: file.path(project_dir, "FINAL_LIST.csv")
#' Where should completed valid metadata be stored as csv?
#' @param backup_metadata path, default: file.path(project_dir, "BACKUP_LIST.csv").
#' The complete list of all unique runs checked in this pipelines,
#'  even ones deleted earlier. Useful for check of what has been done before.
#' @param temp_metadata  path, default: file.path(project_dir, "next_round_manual.csv").
#' The intermediate file used after curation, but before it is validated. This is
#' the where you update the current new metadata annotations for final approval into the
#' complete_metadata. This syncs automatically to the google_url sheet if included.
#' @param blacklist path, Which studies to ignore for this config,
#' default: file.path(project_dir, "BLACKLIST.csv")\cr
#' A csv with 1 column id, which gives BioProject IDs
#' @param google_url url or sheet object for google sheet to use. Set to NULL to
#' not use google sheet.
#' @param flags named character vector, with cut points, where can
#' the pipeline continue if it breaks? Is defined in combination with
#' 'steps' argument that defines that actuall function called for each break point.
#' Increasing break points makes the pipeline run faster at the cost of more
#' resource usage.
#' @param flag_steps list, mapping of functions to ids. The names per list element is a function.
#' The character elements inside each list element are the flag name ids for all steps that will
#' be "marked as done" inside that function. Example: In default preset,
#' the pipe_trim_collapse function marks as done both the trim and collapse flags.
#' @param pipeline_steps a list of the functions to actually run, the functions
#' must be named equal to names of the 'flag_step' argument.
#' @param mode = \code{c("online", "local")[1]}. "online" will assume project IDs for
#' online repository (SRA, ENA, PRJ etc). Local means local folders as accessions.
#' @param delete_raw_files logical, default: mode == "online". If online do delete
#' raw fastq files after trim step is done, for local samples do not delete.
#' Set only to TRUE for mode local, if you have backups!
#' @param delete_trimmed_files logical, default: mode == "online".
#' If TRUE deletes the trimmed fasta files.
#' @param delete_collapsed_files logical, default: mode == "online".
#' If TRUE deletes the collapsed fasta files.
#' @param keep_contaminants logical, default FALSE. Do not keep contaminant aligned reads,
#'  if TRUE they are saved in contamination dir. Ignored if "contam" is not in flags to use.
#' @param keep_unaligned_genome logical, default FALSE. Do not keep contaminant aligned reads,
#'  else saved in contamination dir.
#' @param compress_raw_data logical, default FALSE. If TRUE, will compress raw fastq files.
#' @param stop_downloading_new_data_at_drive_usage integer, default 92,
#' percentage value where the drive will stop downloading new data. Set to 101 to
#' disable a cap.
#' @param max_unprocessed_downloads numeric, default 30. Temporarily stop downloading more data
#' if > 30 studies are not done with the first post download processing step
#' (usually trimming)
#' @param accepted_lengths_rpf default c(20, 21, 25:33), which read lengths to pshift.
#' Default is the standard fractions of normal 80S ~28 and the smaller size of ~21.
#' @param reuse_shifts_if_existing for Ribo-seq, reuse shift table called shifting_table.rds
#' in pshifted folder if it is valid (equal number of sample shift tables in file
#' relative toexperiment)
#' @param split_unique_mappers logical, default FALSE.
#' Run for unique mappers only, split out into seperate directory.
#' @param all_mappers logical, default TRUE. Run for all mappers
#' @param min_raw_reads_pshift numeric, default 1e5. Raw read-count
#' threshold (summed across a study's runs) below which a study is
#' recorded (see \code{\link{qc_diagnostics_path}}) as likely
#' too-few-reads for a usable P-shift. Record-only for now: does not by
#' itself skip or block any step (see \code{skip_pshift_on_low_reads}).
#' @param min_alignment_rate_pshift numeric, percent, default 10. Unique
#' genome-mapping-rate threshold below which a study is recorded as
#' likely wrong-organism. Record-only for now (see
#' \code{skip_pshift_on_low_alignment}).
#' @param max_no_adapter_removed_pct numeric, percent, default 80. Share
#' of reads with no adapter removed at all (fastp), above which a study
#' is recorded as likely bad adapter/barcode/UMI trimming. Record-only:
#' this signal is never used to block a step today.
#' @param skip_pshift_on_low_reads logical, default FALSE. Placeholder
#' for a future pipeline change: does not currently skip anything, even
#' when \code{min_raw_reads_pshift} is tripped.
#' @param skip_pshift_on_low_alignment logical, default FALSE. Placeholder
#' for a future pipeline change: does not currently skip anything, even
#' when \code{min_alignment_rate_pshift} is tripped.
#' @param parallel_conf a bpoptions object, default:
#' \code{bpoptions(log =TRUE,
#'  jobname = "pipeline_step",
#'  logdir = file.path(project_dir, "log_pipeline"),
#'  stop.on.error = TRUE)}
#' Specific pipeline config for parallel settings and log directory for BPPARAM_MAIN
#' @param verbose logical, default TRUE, give start up message
#' @param threads named list, worker counts per step, read via
#' \code{\link{bpparam_from_config}(config, step)}. Every step started
#' with \code{bpparam_from_config()} must have a matching name here
#' (checked there, not here). \code{collapse} defaults to
#' \code{min(threads_default, 16)} -- 16 has worked well in practice as
#' an upper bound, but the real worker count used is further capped at
#' runtime by available memory (see \code{memory_safe_worker_count()},
#' R/pipeline_collapse.R) since collapsing large fastq files is memory-
#' heavy per worker. \code{pshifted}/\code{valid_pshift}/\code{pcounts}
#' default to \code{min(threads_default, threads_blas_cap)} for the
#' same kind of reason, see \code{threads_blas_cap}'s own roxygen.
#' @param threads_blas_cap numeric, default 32. Upper bound used for
#' \code{pshifted}/\code{valid_pshift}/\code{pcounts} in \code{threads}
#' below (the steps that otherwise default to the FULL
#' \code{threads_default} worker count, unlike \code{trim}/\code{collapse}
#' which already have their own smaller fixed caps). Confirmed live,
#' 2026-10-08: running \code{threads_default} (46) forked workers of a
#' stage-group that is ITSELF one of several workers forked by
#' \code{run_pipeline()}'s own main-level dispatch multiplies out to
#' hundreds of OS-level processes+threads very quickly -- each one
#' ALSO starting its own uncapped OpenBLAS thread pool (as many threads
#' as detected cores, by default) -- and can exhaust the container's
#' cgroup \code{pids.max} outright (\code{pthread_create() ->
#' EAGAIN "Resource temporarily unavailable"}, reproduced live; this
#' is a real, separate failure mode from the actual hang that prompted
#' the investigation, but a closely related one, and 46 is unsafe for
#' the same underlying reason either way). A real-data benchmark
#' (\code{shift_qc_cached()}, 30 samples, PRJNA637713-zea_mays) found
#' 32 workers fastest in practice anyway (180s vs 46 workers' 201-300s
#' across two repeated runs, and vs 371s at 8 workers) -- this is not a
#' safety-for-speed tradeoff, 32 wins on both counts. See also
#' \code{blas_set_num_threads()}/\code{omp_set_num_threads()}
#' (RhpcBLASctl) below, called once here with this SAME value: capping
#' the worker COUNT alone isn't sufficient on its own, since each
#' forked worker can ALSO independently try to spin up its own
#' full-width BLAS thread pool regardless of how many sibling workers
#' exist.
#' @param discord_webhook = discord_connection_default_cached()
#' @param BPPARAM_MAIN BiocParallel::MulticoreParam(length(pipeline_steps))
#' The main parallel backend for pipeline, specifying logging behavoir etc.
#' @param BPPARAM_TRIM BiocParallel::MulticoreParam(max(BiocParallel::bpworkers(), 8)),
#' number of cores/threads to use for trimming. Optimal is 8 for most data.
#' @param BPPARAM = bpparam(), number of cores/threads to use.
#' @return a list with a defined config
#' @export
pipeline_config <- function(project_dir = file.path(dirname(config)[1], "NGS_pipeline"),
                            config = ORFik::config(),
                            complete_metadata = file.path(project_dir, "FINAL_LIST.csv"),
                            backup_metadata = file.path(project_dir, "BACKUP_LIST.csv"),
                            temp_metadata = file.path(project_dir, "next_round_manual.csv"),
                            blacklist = file.path(project_dir, "BLACKLIST.csv"),
                            google_url = default_sheets(project_dir),
                            preset = "Ribo-seq",
                            flags = pipeline_flags(project_dir, mode, preset, contam),
                            flag_steps = flag_grouping(flags),
                            pipeline_steps = lapply(names(flag_steps),
                                                    function(x) get(x, mode = "function")),
                            mode = c("online", "local")[1],
                            contam = FALSE,
                            delete_raw_files = mode == "online",
                            delete_trimmed_files = mode == "online",
                            delete_collapsed_files = mode == "online",
                            keep_contaminants = FALSE,
                            keep_unaligned_genome = FALSE,
                            compress_raw_data = FALSE,
                            stop_downloading_new_data_at_drive_usage = 92,
                            max_unprocessed_downloads = 30,
                            accepted_lengths_rpf = c(20, 21, 25:33),
                            reuse_shifts_if_existing = TRUE,
                            split_unique_mappers = FALSE,
                            all_mappers = TRUE,
                            min_raw_reads_pshift = 1e5,
                            min_alignment_rate_pshift = 10,
                            max_no_adapter_removed_pct = 80,
                            skip_pshift_on_low_reads = FALSE,
                            skip_pshift_on_low_alignment = FALSE,
                            parallel_conf = bpoptions(log = TRUE,
                                                      jobname = "pipeline_step",
                                                      logdir = file.path(project_dir, "log_pipeline"),
                                                      stop.on.error = TRUE),
                            discord_webhook = discord_connection_default_cached(),
                            verbose = TRUE,
                            thread_type = BiocParallel::MulticoreParam,
                            threads_default = BiocParallel::bpworkers(),
                            threads_blas_cap = 32,
                            threads = list(main = length(pipeline_steps),
                                           default = threads_default,
                                           trim = min(threads_default, 8),
                                           collapse = min(threads_default, 16),
                                           pshifted = min(threads_default, threads_blas_cap),
                                           valid_pshift = min(threads_default, threads_blas_cap),
                                           pcounts = min(threads_default, threads_blas_cap)
                                           )) {
  # Caps each forked worker's OWN internal BLAS/OpenMP thread pool to
  # the same bound as threads_blas_cap above -- see that parameter's
  # own roxygen for why the worker COUNT cap alone isn't sufficient.
  # Set once, here, in the parent: confirmed live that this setting
  # survives fork() and is correctly honored by every forked
  # descendant without needing to be re-called inside each one
  # (RhpcBLASctl::blas_set_num_threads() sets a plain thread-count
  # value, not a live thread-pool handle, so copy-on-write fork
  # semantics carry it over cleanly -- unlike the thread pool itself,
  # which is exactly why OpenBLAS's own uncapped-by-default pool is
  # unsafe to inherit through a fork in the first place).
  RhpcBLASctl::blas_set_num_threads(threads_blas_cap)
  RhpcBLASctl::omp_set_num_threads(threads_blas_cap)
  if (verbose) {
    message("Setting up mNGSp config..")
    message("- Preset (mode): ", preset, " (", mode, ")")
    message("- Workers (threads): ", threads$default, " (", get_fun_name(thread_type, "BiocParallel"), ")")
    message("- Sync to google: ", !is.null(google_url))
    message("- Metadata dir: ", project_dir)
    message("- Output data dir: ", config["bam"])
  }

  stopifnot(mode %in% c("online", "local"))
  stopifnot(is.logical(delete_raw_files) & is.logical(delete_trimmed_files) &
            is.logical(delete_collapsed_files) & is.logical(keep_contaminants) &
            is.logical(keep_unaligned_genome) & is.logical(compress_raw_data) &
            is.logical(split_unique_mappers) & is.logical(all_mappers))
  stopifnot(all_mappers | split_unique_mappers)
  stopifnot(is(pipeline_steps, "list"))
  if (preset == "empty") message("Using empty preset, now add your steps: using add_step_to_pipeline()")
  if (preset != "empty") stopifnot(is(pipeline_steps[[1]], "function"))
  stopifnot(length(pipeline_steps) == length(flag_steps))

  delete_any_files <- c(delete_raw_files, delete_trimmed_files, delete_collapsed_files)
  names(delete_any_files) <- c("raw", "trimmed", "collapsed")
  if (mode == "local" & any(delete_any_files))
    message("You have run local fastq files with delete ",
    paste(names(delete_any_files)[delete_any_files], collapse = " & "), " files,",
            "are you sure this is what you want?")



  return(list(project = project_dir, config = config, flag = flags,
              flag_steps = flag_steps,
              pipeline_steps = pipeline_steps,
              metadata = paste0(project_dir, "/metadata"),
              complete_metadata = complete_metadata,
              backup_metadata = backup_metadata,
              temp_metadata = temp_metadata,
              blacklist = blacklist,
              google_url = google_url, mode = mode,
              delete_raw_files = delete_raw_files,
              delete_trimmed_files = delete_trimmed_files,
              delete_collapsed_files = delete_collapsed_files,
              keep_contaminants = keep_contaminants,
              keep.unaligned.genome = keep_unaligned_genome,
              compress_raw_data = compress_raw_data,
              stop_downloading_new_data_at_drive_usage = stop_downloading_new_data_at_drive_usage,
              max_unprocessed_downloads = max_unprocessed_downloads,
              accepted_lengths_rpf = accepted_lengths_rpf,
              reuse_shifts_if_existing = reuse_shifts_if_existing,
              split_unique_mappers = split_unique_mappers,
              all_mappers = all_mappers,
              min_raw_reads_pshift = min_raw_reads_pshift,
              min_alignment_rate_pshift = min_alignment_rate_pshift,
              max_no_adapter_removed_pct = max_no_adapter_removed_pct,
              skip_pshift_on_low_reads = skip_pshift_on_low_reads,
              skip_pshift_on_low_alignment = skip_pshift_on_low_alignment,
              preset = preset, parallel_conf = parallel_conf,
              discord_webhook = discord_webhook,
              thread_type = thread_type,
              threads = threads
              ))
}

#' Define which flags are done in each function
#' @param flags character vector, of paths with names as flag steps and
#' attribute called "grouping" with length same as flag with grouping information,
#' where group names are the group functions.
#' @param substeps a list of character vectors, where each sub list
#' is the steps that happen in that sub component of the pipeline
#' @return a list, returns the validated 'substeps' argument.
#' @export
flag_grouping <- function(flags, substeps = preset_grouping(flags)) {
  if (length(flags) == 0) return(list())
  stopifnot(all(unlist(substeps) %in% names(flags)))
  stopifnot(length(unlist(substeps)) > 0)
  stopifnot(length(flags) > 0)
  if (!all(names(flags) %in% unlist(substeps)))
    message("Running pipeline with not all flags set, is this intentional?")
  return(substeps)
}

#' Reconstruct the step-id groupings from a flags vector's grouping attribute
#'
#' `split(names(flags), grouping)`, preserving the grouping's own
#' first-occurrence order (not alphabetical `split()` order).
#' @param flags a flags vector as returned by `pipeline_flags()`: named by
#' step id, with a parallel `"grouping"` attribute naming each step's
#' owning pipe_*() function
#' @return named list, one element per group, each a character vector of
#' step ids
#' @noRd
preset_grouping <- function(flags) {
  grouping <- attr(flags, "grouping")
  if (is.null(names(flags))) stop("flags must have names!")
  if (is.null(grouping)) stop("Your flags does not have defined attr 'grouping',",
                              "it should contain the grouping of flags")
  if (length(flags) != length(grouping)) stop("The grouping attribute must be",
                                              "equal size to number of flags!")
  return(split(names(flags), grouping)[unique(grouping)])
}

#' Find the exported name(s) of a value in a package's namespace
#'
#' Reverse lookup via `identical()`, used only for a startup message
#' naming e.g. the resolved `BiocParallel` param type.
#' @param fun the value to search for
#' @param pkg character, package name to search
#' @return character vector of matching names, `character(0)` if none found
#' @noRd
get_fun_name <- function(fun, pkg) {
  nm <- ls(getNamespace(pkg))
  nm[vapply(nm, function(n)
    identical(get(n, envir = getNamespace(pkg)), fun),
    logical(1)
  )]
}
