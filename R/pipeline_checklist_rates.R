# Live processing-speed figures for pipeline_checklist() -- built from a
# real session's own stall (a ~70-minute silent hang with no visible
# symptom until investigated by hand, see R/star_index_lock.R) showing
# the cost of having a done/total count but no way to tell whether that
# count is actually still moving at a normal pace.

#' Package-private, in-memory "last seen" state for rate_since_first_seen().
#' Lives for the lifetime of whichever worker process calls into this --
#' each stage-group (pipe_fetch(), pipe_align_clean(), ...) runs as its
#' own long-lived forked worker across a run_pipeline() session, so this
#' persists correctly across that worker's own repeated on_poll() ticks
#' without needing a file or any cross-process coordination.
#' @noRd
.rate_state <- new.env(parent = emptyenv())

#' Average rate of work since `key` was first observed, in units/hour
#'
#' Deliberately an AVERAGE since first observation, not a delta since the
#' last call: run_experiment_subprocess()'s own poll_interval defaults to
#' 10s, far too short a window for any stage where one unit of work takes
#' minutes (trim, align, ...) -- a last-tick delta would mostly read 0,
#' then spike. An average since first-seen is stable, and still resets
#' itself correctly the moment the active identifier changes, since a new
#' identifier is simply a new, neverbefore-seen key.
#' @param key character, uniquely identifies what's being measured (e.g.
#' \code{paste(stage_name, active_experiment)}) -- a NEW key always
#' starts fresh at NA on its first call, by design.
#' @param amount_now numeric, the current cumulative amount of work done
#' @param now POSIXct, default Sys.time() (parameterized for tests)
#' @return numeric, units/hour, or NA_real_ on the first observation of
#' this key, or when elapsed time is <= 0. Never negative -- a decrease
#' (e.g. a flag or file got reset) is stale-baseline noise, not a real
#' negative rate, so it's clamped to 0.
#' @noRd
rate_since_first_seen <- function(key, amount_now, now = Sys.time()) {
  baseline <- .rate_state[[key]]
  if (is.null(baseline)) {
    .rate_state[[key]] <- list(amount = amount_now, time = now)
    return(NA_real_)
  }
  dt_hours <- as.numeric(difftime(now, baseline$time, units = "hours"))
  if (dt_hours <= 0) return(NA_real_)
  max(0, (amount_now - baseline$amount) / dt_hours)
}

#' One experiment's conf (named path vector: exp/bam/fastq/...), from
#' the live pipelines object
#' @param pipelines the pipelines list
#' @param experiment character, experiment id
#' @return named character vector, or NULL if not found
#' @noRd
experiment_conf <- function(pipelines, experiment) {
  for (p in pipelines) {
    for (organism in names(p$organisms)) {
      conf <- p$organisms[[organism]]$conf
      if (identical(unname(conf["exp"]), unname(experiment))) return(conf)
    }
  }
  NULL
}

#' Every experiment's full Run-id vector, from the live pipelines object
#'
#' Same iteration shape as \code{\link{experiment_sample_counts}}, just
#' returning the Run ids themselves instead of a count -- needed here to
#' resolve WHICH specific run is currently active for a marker step
#' (fetch/aligned), not just how many samples are done.
#' @param pipelines the pipelines list (see pipeline_init_all())
#' @return named list, names are experiment ids, values are character
#' vectors of Run ids
#' @noRd
experiment_run_ids <- function(pipelines) {
  runs <- list()
  for (p in pipelines) {
    for (organism in names(p$organisms)) {
      exp <- p$organisms[[organism]]$conf["exp"]
      runs[[exp]] <- p$study[ScientificName == organism]$Run
    }
  }
  runs
}

#' Which run is currently active for a marker step, if any
#' @param pipelines the pipelines list
#' @param config the mNGSp config object
#' @param marker_step character, e.g. "fetch" or "aligned"
#' @param experiment character, the active experiment id
#' @return character run id, or NA_character_ if every run is already
#' marked done (e.g. the done-flag for this whole step just hasn't been
#' written yet) or the experiment isn't found
#' @noRd
active_run_id <- function(pipelines, config, marker_step, experiment) {
  all_runs <- experiment_run_ids(pipelines)[[experiment]]
  if (is.null(all_runs)) return(NA_character_)
  done_runs <- samples_done(config, marker_step, experiment)
  remaining <- setdiff(all_runs, done_runs)
  if (length(remaining) == 0) return(NA_character_)
  remaining[1]
}

#' Fetch step's live download rate for one in-progress run, in MB/s
#'
#' Deliberately extension-agnostic: sums the size of every file matching
#' the run's accession under its fastq output dir, whatever its current
#' extension. Transparently covers whichever of download_sra()'s 3 paths
#' is in play (the raw .sra growing during download_sra_aws()/_ascp(), or
#' the final .fastq(.gz) growing during fasterq-dump's own conversion)
#' without needing to track which phase is currently active.
#' @param fastq_dir character, the experiment's raw fastq output dir
#' (\code{conf["fastq"]})
#' @param run character, the run id currently downloading
#' @return numeric, MB/s, or NA_real_ if nothing matching exists yet, or
#' this is the first observation (see rate_since_first_seen())
#' @noRd
fetch_progress_rate <- function(fastq_dir, run) {
  files <- list.files(fastq_dir, pattern = run, full.names = TRUE)
  if (length(files) == 0) return(NA_real_)
  bytes_now <- sum(file.info(files)$size, na.rm = TRUE)
  bytes_per_hour <- rate_since_first_seen(paste("fetch", run), bytes_now)
  if (is.na(bytes_per_hour)) return(NA_real_)
  bytes_per_hour / 3600 / 1e6
}

#' STAR's own reported alignment speed for one in-progress run, in M reads/hr
#'
#' STAR already computes and reports this directly in Log.progress.out
#' (column 2, "Speed", M reads/hr) -- read the last data row rather than
#' computing our own rate. Appends roughly once a minute, so a small/fast
#' input can have zero data rows yet; that's a normal "not available yet"
#' case, not an error (matches this file's own documented behavior).
#' @param bam_dir character, the experiment's bam dir (\code{conf["bam"]})
#' @param run character, the run id currently aligning
#' @return numeric, M reads/hr, or NA_real_ if the log doesn't exist yet,
#' or has no data rows yet
#' @noRd
align_progress_rate <- function(bam_dir, run) {
  candidates <- c(
    Sys.glob(file.path(bam_dir, "aligned", "LOGS_SINGLE", paste0("*", run, "*_Log.progress.out"))),
    Sys.glob(file.path(bam_dir, "aligned", "LOGS_PAIRED", paste0("*", run, "*_Log.progress.out")))
  )
  if (length(candidates) == 0) return(NA_real_)
  lines <- readLines(candidates[1], warn = FALSE)
  lines <- lines[!grepl("^ALL DONE!", lines)]
  data_lines <- lines[-seq_len(min(2, length(lines)))] # drop the 2 header lines
  if (length(data_lines) == 0) return(NA_real_)
  # Field 4, not 2: STAR's own "Time" column ("Jun 01 14:43:03") is
  # itself 3 whitespace-separated tokens, shifting every later column
  # right by 2 -- confirmed directly against a real Log.progress.out.
  last_fields <- strsplit(trimws(data_lines[length(data_lines)]), " +")[[1]]
  speed <- suppressWarnings(as.numeric(last_fields[4]))
  if (is.na(speed)) return(NA_real_)
  speed
}

#' Generic samples/hour (or studies/hour) rate for every other stage
#'
#' Covers both cases already computed by pipeline_checklist()'s own
#' per-stage loop: stages with per-sample marker detail (amount =
#' active_done, keyed by stage+experiment) and stages with only
#' study-level counts (amount = done, keyed by stage alone, since there
#' is no single "active experiment" concept for a stage whose work hands
#' off to ORFik's own internal BiocParallel dispatch or a whole-folder
#' call).
#' @param stage_name character
#' @param active_experiment character or NA -- NA selects the
#' study-level/studies-per-hour case
#' @param amount numeric, active_done (sample case) or done (study case)
#' @return numeric, units/hour, or NA_real_ on first observation
#' @noRd
generic_progress_rate <- function(stage_name, active_experiment, amount) {
  key <- if (is.na(active_experiment)) stage_name else paste(stage_name, active_experiment)
  rate_since_first_seen(key, amount)
}
