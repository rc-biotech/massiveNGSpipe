# Optional first stage for pooled (multiplexed) libraries: pipeline_config(demultiplex = TRUE).
#
# Pools live in <raw_data>/<pool>/multiplexed/ with a controlfile.csv (see ?demultiplex). A pool is matched to the
# metadata by run name: every control-file final_sample (default: sample) that is a Run of the curated metadata
# (FINAL_LIST) is delivered to that Run's experiment, so one pool may feed several studies and organisms. The
# experiment's organism, study and labels therefore come from the metadata, as for any other run.
#
# Delivery per Run: <exp bam>/trim/SINGLE/collapsed_trimmed_<Run>.fasta.gz and <exp bam>/trim/report_<Run>.json,
# then the per-sample "trim" and "collapsed" markers, so the normal trim/collapse stage skips these runs and
# alignment proceeds as usual. Experiments without pooled runs get the "demultiplexed" flag immediately and are
# trimmed from their own raw files. Pools are demultiplexed once into a staging directory
# (<processed_data>/demultiplexed_pools/<pool>/), keyed by the control file checksum.

#' Pipeline stage: demultiplex pooled lanes and deliver runs to their experiments
#' @param pipelines list, output of pipeline_init_all
#' @param config list, output of pipeline_config(demultiplex = TRUE)
#' @return invisible(NULL)
#' @noRd
pipe_demultiplex <- function(pipelines, config) {
  pools <- demultiplex_pools(config)
  run_map <- demultiplex_run_map(pipelines)
  for (pipeline in pipelines) {
    if (file.exists(report_failed_pipe_path(config, pipeline$accession))) next
    if (stop_requested(config)) break
    try <- try({
      for (organism in names(pipeline$organisms)) {
        exp <- pipeline$organisms[[organism]]$conf["exp"]
        if (!step_is_next_not_done(config, "demultiplexed", exp)) next
        runs <- pipeline$study[ScientificName == organism]$Run
        fed <- pools[final_sample %in% runs]
        done <- samples_done(config, "collapsed", exp)
        for (p in unique(fed[!(final_sample %in% done)]$pool_dir)) demultiplex_pool_deliver(p, pools, run_map, config)
        missing <- setdiff(fed$final_sample, samples_done(config, "collapsed", exp))
        if (length(missing)) stop("Pooled runs not delivered: ", paste(missing, collapse = ", "))
        set_flag(config, "demultiplexed", exp)
      }
    })
    report_failed_pipe(try, config, "demultiplexed", pipeline$accession)
  }
  invisible(NULL)
}

#' All registered pools: <raw_data>/*/multiplexed/controlfile.csv
#' @return data.table(pool_dir, controlfile, final_sample)
#' @noRd
demultiplex_pools <- function(config) {
  raw_root <- unname(config$config["fastq"])
  cf <- Sys.glob(file.path(raw_root, "*", "multiplexed", "controlfile.csv"))
  if (!length(cf)) return(data.table(pool_dir = character(), controlfile = character(), final_sample = character()))
  rbindlist(lapply(cf, function(f) {
    x <- fread(f, colClasses = "character")
    fs <- if ("final_sample" %in% colnames(x)) ifelse(is.na(x$final_sample) | x$final_sample == "", x$sample, x$final_sample) else x$sample
    data.table(pool_dir = dirname(f), controlfile = f, final_sample = unique(fs))
  }))
}

#' Run -> experiment directories over all pipelines
#' @noRd
demultiplex_run_map <- function(pipelines) {
  rbindlist(lapply(pipelines, function(pipeline) rbindlist(lapply(names(pipeline$organisms), function(organism) {
    conf <- pipeline$organisms[[organism]]$conf
    data.table(Run = pipeline$study[ScientificName == organism]$Run, exp = unname(conf["exp"]), bam = unname(conf["bam"]))
  }))))
}

#' Demultiplex one pool (once per control-file version) and deliver its runs
#' @noRd
demultiplex_pool_deliver <- function(pool_dir, pools, run_map, config) {
  pool <- basename(dirname(pool_dir))
  staging <- file.path(unname(config$config["bam"]), "demultiplexed_pools", pool)
  dir.create(staging, showWarnings = FALSE, recursive = TRUE)
  lock <- file.path(staging, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    owner <- suppressWarnings(as.integer(readLines(file.path(lock, "pid"), warn = FALSE)))
    if (length(owner) && !is.na(owner) && dir.exists(file.path("/proc", owner))) stop("Pool ", pool, " is being demultiplexed by pid ", owner)
    unlink(lock, recursive = TRUE); dir.create(lock)
  }
  writeLines(as.character(Sys.getpid()), file.path(lock, "pid"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)

  ctrl_file <- file.path(pool_dir, "controlfile.csv")
  key <- unname(tools::md5sum(ctrl_file))
  done_file <- file.path(staging, "demultiplex_done.json")
  done <- file.exists(done_file) && identical(jsonlite::read_json(done_file)$controlfile_md5, key)
  if (!done) {
    message("Demultiplexing pool ", pool)
    demultiplex(pool_dir, out_dir = staging, tmp_dir = file.path(config$ssd_scratch_dir %||% tempdir(), paste0("demultiplex_", pool)),
                parallel_files = config$demultiplex_parallel_files %||% 4, collapse_workers = config$demultiplex_collapse_workers %||% 8)
    jsonlite::write_json(list(controlfile_md5 = key, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), done_file, auto_unbox = TRUE)
  }
  # deliver every curated run of this pool (runs absent from the metadata stay in staging)
  this_pool <- pool_dir
  for (fs in pools[pool_dir == this_pool]$final_sample) {
    target <- run_map[Run == fs]
    if (!nrow(target)) { message("- ", fs, ": not in the curated metadata, left in ", staging); next }
    if (nrow(target) > 1) stop("Run ", fs, " belongs to several experiments: ", paste(target$exp, collapse = ", "))
    src_fa <- file.path(staging, "SINGLE", paste0("collapsed_trimmed_", fs, ".fasta.gz"))
    src_js <- file.path(staging, paste0("report_", fs, ".json"))
    trim_dir <- file.path(target$bam, "trim")
    dst_fa <- file.path(trim_dir, "SINGLE", basename(src_fa))
    if (fs %in% samples_done(config, "collapsed", target$exp) && file.exists(dst_fa)) next
    if (!file.exists(src_fa)) stop("Demultiplexed output missing for ", fs, ": ", src_fa)
    dir.create(file.path(trim_dir, "SINGLE"), showWarnings = FALSE, recursive = TRUE)
    if (!file.copy(src_fa, dst_fa, overwrite = TRUE)) stop("Could not deliver ", src_fa)
    file.copy(src_js, file.path(trim_dir, basename(src_js)), overwrite = TRUE)
    set_sample_flag(config, "trim", target$exp, fs)
    set_sample_flag(config, "collapsed", target$exp, fs)
  }
  invisible(staging)
}
