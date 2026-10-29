#' Pipeline error records, event stream and agent utilities
#'
#' Every per-study failure is recorded in three places:
#' \itemize{
#'   \item \code{<error_dir>/<exp>.rds}: the record \code{report_failed_pipe()} has always written (one character
#'   string). Its existence is the error flag: stage loops skip the study while it exists, the checklist shows it,
#'   and deleting it (\code{\link{clear_pipeline_error}}) lets a running pipeline retry the study.
#'   \item \code{<project>/error_logs/events.jsonl}: one JSON line per event (\code{"error"}, \code{"resolved"},
#'   \code{"crash"}), appended across sessions. Meant for watchers (\code{tail -F}, or
#'   \code{inst/scripts/watch_pipeline_errors.py}) so a person or an agent is notified instead of polling files.
#'   \item an optional hook: \code{config$error_hook} (an R function of one argument, the event list) or the
#'   environment variable \code{MNGSP_ERROR_HOOK} (a shell command run with \code{MNGSP_EVENT}, \code{MNGSP_EXP},
#'   \code{MNGSP_STEP}, \code{MNGSP_MESSAGE}, \code{MNGSP_HINT}, \code{MNGSP_PROJECT} and \code{MNGSP_ERROR_FILE}
#'   set). Hook failures are reported as warnings and never stop the pipeline.
#' }
#' @name pipeline_errors
NULL

#' Likely cause of a pipeline error message
#'
#' Pattern table for errors seen on this server; first match wins. Used in error records, events and
#' \code{\link{pipeline_errors}}.
#' @param msg character, error message(s)
#' @return character, a short hint, or \code{"unknown"}
#' @export
#' Disk-full hint from the drives themselves
#'
#' Error texts from a full disk are often empty or misleading (a worker killed
#' while writing, "cannot create R_TempDir"); checking free space at error time
#' is more reliable. Checks tempdir() and the raw/processed data drives.
#' @param config the mNGSp config object
#' @param min_gb numeric, default 2: below this a drive counts as full
#' @return character hint or NULL
#' @noRd
low_disk_hint <- function(config, min_gb = 2) {
  dirs <- c(temp = tempdir(), raw = unname(config$config["fastq"]), processed = unname(config$config["bam"]))
  dirs <- dirs[!is.na(dirs) & dir.exists(dirs)]
  free <- vapply(dirs, tmp_free_gb, numeric(1))
  full <- !is.na(free) & free < min_gb
  if (!any(full)) return(NULL)
  paste0("disk full: ", paste(sprintf("%s drive (%s) %.1f GB free", names(dirs)[full], dirs[full], free[full]),
                              collapse = "; "))
}

error_hint <- function(msg) {
  hints <- c(
    "could not open genome file|genomeParameters.txt" =
      "STAR index missing or incomplete: rebuild the index (reference folder STAR_index/genomeDir)",
    "STAR alignment step failed|EXITING because of FATAL ERROR" =
      "STAR failed: see the sample's aligned/LOGS/*Log.out (index, RAM, input FASTQ)",
    "storage exhausted|No space left on device|cannot write|disk full" =
      "disk full: check df for the raw/processed drive and the temp directory",
    "cannot allocate|std::bad_alloc|Killed|out of memory|vector memory exhausted" =
      "out of memory: lower threads for this step or wait for RAM",
    "no periodic or accepted read-lengths|Not enough reads in TIS region" =
      "no 3-nt periodicity: check adapter/barcode trimming, read lengths and whether the library is Ribo-seq",
    "detect_strand_mode" = "RNA strand detection failed: often a broken STAR index or too few mapped reads",
    "RNA orientation unresolved" =
      "RNA strand unresolved: see aligned/QC_STATS/rna_strand.csv. ~0 informative reads even after the whole-BAM fallback = few reads on annotated genes: check alignment rate and assay (e.g. shRNA/CRISPR screen tags, DNA; PRJNA418238 2026-10-20) before excluding; sense fraction between thresholds = moderately stranded -> per-run config$rna_strand_by_run",
    "dbGaP|Access denied|unauthorized while resolving" =
      "controlled-access data (dbGaP): cannot be downloaded; set KEEP = FALSE for the study's runs",
    "munmap_chunk|invalid pointer|double free|free\\(\\): invalid|Segmentation fault|core dumped" =
      "a compiled helper (trimcollapse/demultiplexer) crashed: report it with the input file",
    "fasterq-dump|prefetch|sra_to_fastq|aws.*s3" =
      "download failed: retry; check the accession, network and temp space",
    "lazy-load|lazy load|cannot open compressed file.*\\.rdb" =
      "package files changed under a running session (reinstall during run): restart the run",
    "Could not resolve host|Timeout was reached|HTTP error" = "network error: retry later",
    "gargle|token|oauth|googlesheets4|googledrive" = "Google authentication problem",
    "cannot create R_TempDir|TMPDIR" = "temporary directory not writable or full",
    "object '.*' not found|could not find function" =
      "code error (missing object or function): likely a bug or a version mismatch, report it",
    "subscript out of bounds|missing value where TRUE/FALSE needed|arguments imply differing number of rows" =
      "code or metadata edge case: inspect the study's metadata and report it",
    "line is too long|read_fastq_files|invalid FASTQ|FASTQ.*(truncated|malformed)" =
      "corrupt FASTQ (e.g. a block of NUL bytes from an interrupted write, PRJNA630198 2026-10-20): check the file, delete it and the run's fetch/trim markers, re-fetch",
    "unevaluated and other errors" =
      "a parallel worker died without an R error (often killed by the OOM killer or a signal)")
  msg <- paste(msg, collapse = "\n")
  for (p in names(hints)) if (grepl(p, msg, ignore.case = TRUE)) return(unname(hints[p]))
  "unknown"
}

#' Short one-line form of an error, with the first BiocParallel remote error pulled out
#' @param err a try-error, condition or character
#' @return character(1)
#' @noRd
error_one_line <- function(err) {
  txt <- if (inherits(err, "condition")) conditionMessage(err) else paste(as.character(err), collapse = "\n")
  lines <- trimws(strsplit(txt, "\n", fixed = TRUE)[[1]])
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return("(empty error message)")
  first <- sub("^Error( in [^:]*)? *: *", "", lines[1])
  remote <- grep("^Error", lines[-1], value = TRUE)
  if (grepl("BiocParallel errors", txt) && length(remote))
    first <- paste0(first, " | first remote error: ", sub("^Error( in [^:]*)? *: *", "", remote[1]))
  substr(first, 1, 500)
}

#' Path of a project's error event stream
#' @param config the mNGSp config object (needs \code{config$project})
#' @return character path \code{<project>/error_logs/events.jsonl}
#' @export
error_events_path <- function(config) file.path(config$project, "error_logs", "events.jsonl")

#' Append one event to the error event stream and run the error hook
#' @param config the mNGSp config object
#' @param event list with at least \code{event}, \code{exp}, \code{step}, \code{message}
#' @return invisible(event)
#' @noRd
error_event_emit <- function(config, event) {
  event <- c(list(time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                  project = basename(config$project %||% ""),
                  session = basename(config$error_dir %||% ""), pid = Sys.getpid()), event)
  event <- lapply(event, function(v) if (is.null(v) || length(v) == 0) "" else as.character(v)[1])
  if (!is.null(config$project)) tryCatch({
    path <- error_events_path(config)
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    line <- as.character(jsonlite::toJSON(event, auto_unbox = TRUE))
    # one write() of a single line: appends from parallel workers do not interleave within a line
    cat(line, "\n", sep = "", file = path, append = TRUE)
  }, error = function(e) warning("could not append error event: ", conditionMessage(e)))
  error_hook_run(config, event)
  invisible(event)
}

#' Run the configured error hook (R function or MNGSP_ERROR_HOOK shell command); never throws
#' @noRd
error_hook_run <- function(config, event) {
  hook <- config$error_hook
  if (is.function(hook)) {
    tryCatch(hook(event), error = function(e) warning("error_hook failed: ", conditionMessage(e)))
  }
  cmd <- if (is.character(hook) && length(hook) == 1 && nzchar(hook)) hook else Sys.getenv("MNGSP_ERROR_HOOK")
  if (nzchar(cmd)) {
    env <- c(MNGSP_EVENT = event$event, MNGSP_EXP = event$exp, MNGSP_STEP = event$step,
             MNGSP_MESSAGE = event$message, MNGSP_HINT = event$hint %||% "", MNGSP_PROJECT = event$project,
             MNGSP_ERROR_FILE = event$error_file %||% "")
    tryCatch({
      old <- Sys.getenv(names(env), unset = NA, names = TRUE)
      do.call(Sys.setenv, as.list(env))
      on.exit({ for (n in names(old)) if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, setNames(list(old[[n]]), n)) })
      # subshell: ignore.stdout appends ">/dev/null", which would override a redirect inside the hook command
      status <- system(paste0("(", cmd, ")"), wait = TRUE, ignore.stdout = TRUE, timeout = 30)
      if (!identical(status, 0L)) warning("MNGSP_ERROR_HOOK exited with status ", status)
    }, error = function(e) warning("MNGSP_ERROR_HOOK failed: ", conditionMessage(e)))
  }
  invisible(NULL)
}

#' Record a per-study failure: error flag file, event stream and hook
#'
#' Writes the same \code{<error_dir>/<exp>.rds} record as before (format unchanged:
#' \code{"<exp> ( <step> ) <error text> <time>"}, now with a \code{Likely cause:} line), then emits an
#' \code{"error"} event.
#' @param config the mNGSp config object
#' @param exp character, experiment or study id
#' @param step character, step name
#' @param err a try-error, condition or character
#' @return invisible(path of the record, or NA)
#' @noRd
error_record_write <- function(config, exp, step, err) {
  txt <- if (inherits(err, "condition")) conditionMessage(err) else paste(as.character(err), collapse = "\n")
  # Subprocess errors often only say "see error above": the real cause (e.g. STAR "could not open genome file")
  # is in this session's console logs, so the matching lines are added to the record and used for the hint.
  console <- tryCatch(console_error_excerpt(config, exp), error = function(e) character())
  # "0 remote errors, N unevaluated": the worker's own R error is often only in the BiocParallel task logs
  if (grepl("unevaluated and other errors", txt))
    console <- c(console, tryCatch(bp_task_error_excerpt(config, exp), error = function(e) character()))
  hint <- error_hint(c(txt, console))
  low <- tryCatch(low_disk_hint(config), error = function(e) NULL)
  if (!is.null(low)) hint <- low
  path <- NA_character_
  if (!is.null(config$error_dir)) {
    path <- report_failed_pipe_path(config, exp)
    dir.create(config$error_dir, showWarnings = FALSE, recursive = TRUE)
    console_txt <- if (length(console)) paste0("\nConsole log excerpt:\n", paste(console, collapse = "\n")) else ""
    record <- paste(exp, "(", step, ")", txt, console_txt, paste0("\nLikely cause: ", hint, "\n"), paste(Sys.time()))
    saveRDS(record, path)
  }
  error_event_emit(config, list(event = "error", exp = exp, step = step, message = error_one_line(err),
                                hint = hint, error_file = path))
  invisible(path)
}

#' Error lines from this session's console logs for one experiment or study
#'
#' Looks in \code{<session log base>/console/<stage>/<exp>*.{err,out}.log} (written by
#' \code{run_experiment_subprocess()}) and returns the last few distinct lines that look like errors.
#' @param config the mNGSp config object
#' @param exp character, experiment name or study accession (prefix match)
#' @param n integer, max lines
#' @return character vector (possibly empty)
#' @noRd
console_error_excerpt <- function(config, exp, n = 6) {
  if (is.null(config$session_dir) && is.null(config$project)) return(character())
  base <- file.path(pipeline_log_base(config), "console")
  if (!dir.exists(base)) return(character())
  files <- list.files(base, pattern = paste0("^", gsub("([.()+])", "\\\\\\1", exp), ".*\\.(err|out)\\.log$"),
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(character())
  lines <- unlist(lapply(files, function(f) utils::tail(readLines(f, warn = FALSE), 2000)), use.names = FALSE)
  hits <- grep("FATAL|ERROR|Error|error:|failed|Killed|exhausted|cannot|No such file", lines, value = TRUE)
  hits <- unique(trimws(hits)); hits <- hits[nzchar(hits)]
  utils::tail(substr(hits, 1, 300), n)
}

#' ERROR lines from this session's BiocParallel task logs that concern one experiment
#'
#' With \code{bpstopOnError} a worker's R error can be missing from the returned BiocParallel error ("0 remote
#' errors, N unevaluated", PRJNA630198 collapse 2026-10-20: ShortRead "line is too long" on a corrupt FASTQ, misread
#' as a killed worker). The task logs (\code{<parallel_conf$logdir>/<step>/<step>.task<i>.log}) keep it. The log
#' directory is shared by concurrent pipeline sessions, so only logs that mention \code{exp} and were written since
#' this session started are read, and only \code{ERROR [time]} lines from that time on are returned.
#' @param config the mNGSp config object
#' @param exp character, experiment name or study accession
#' @param n integer, max lines
#' @return character vector (possibly empty)
#' @noRd
bp_task_error_excerpt <- function(config, exp, n = 4) {
  logdir <- config$parallel_conf$logdir
  if (is.null(logdir) || is.na(logdir) || !dir.exists(logdir)) return(character())
  since <- suppressWarnings(as.POSIXct(basename(config$error_dir %||% ""), tz = ""))
  if (is.na(since)) since <- Sys.time() - 24 * 3600
  files <- list.files(logdir, pattern = "\\.task[0-9]+\\.log$", recursive = TRUE, full.names = TRUE)
  files <- files[file.mtime(files) >= since]
  hits <- unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    if (!any(grepl(exp, lines, fixed = TRUE))) return(character())
    err <- grep("^ERROR \\[", lines, value = TRUE)
    t <- suppressWarnings(as.POSIXct(sub("^ERROR \\[([^]]+)\\].*", "\\1", err), tz = ""))
    err[!is.na(t) & t >= since]
  }), use.names = FALSE)
  utils::tail(substr(unique(hits), 1, 300), n)
}

#' Record an unexpected crash of a stage loop or of run_pipeline() itself
#'
#' These errors are not tied to one study, so they are written as \code{<session_dir>/crash_<where>.rds}
#' (with the call stack) and emitted as a \code{"crash"} event. Called from inside the forked stage worker, so the
#' record survives even when the worker's console output is lost.
#' @param config the mNGSp config object
#' @param where character, stage name or "run_pipeline"
#' @param cond the condition
#' @param calls optional call stack (from \code{sys.calls()} inside a calling handler)
#' @return invisible(path)
#' @noRd
pipeline_crash_record <- function(config, where, cond, calls = NULL) {
  msg <- if (inherits(cond, "condition")) conditionMessage(cond) else as.character(cond)
  cls <- paste(class(cond), collapse = "/")
  stack <- if (length(calls)) vapply(calls, function(x) paste(deparse(x, nlines = 1L), collapse = ""), "") else character()
  path <- NA_character_
  base <- config$session_dir %||% config$error_dir
  if (!is.null(base)) {
    dir.create(base, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(base, paste0("crash_", gsub("[^A-Za-z0-9_.-]", "_", where), ".rds"))
    saveRDS(list(where = where, class = cls, message = msg, hint = error_hint(msg), calls = utils::tail(stack, 40),
                 time = Sys.time(), pid = Sys.getpid()), path)
  }
  message("massiveNGSpipe: unexpected ", cls, " in ", where, ": ", msg, "\nLikely cause: ", error_hint(msg),
          if (!is.na(path)) paste0("\nRecord: ", path) else "")
  error_event_emit(config, list(event = "crash", exp = "", step = where, message = error_one_line(cond),
                                hint = error_hint(msg), error_file = path))
  invisible(path)
}

#' Current pipeline errors (one row per error flag file)
#'
#' Reads the error flag files of one session. Only files that still exist are returned, so an error removed with
#' \code{\link{clear_pipeline_error}} disappears here and from the checklist.
#' @param config the mNGSp config object. With \code{config$error_dir} set (inside a run) that session is used,
#' otherwise the session selected by \code{index}.
#' @param index integer, session index, newest first (ignored when \code{config$error_dir} is set)
#' @return data.table with columns exp, step, message, hint, time, file (zero rows when there are no errors)
#' @export
pipeline_errors <- function(config, index = 1) {
  dir <- config$error_dir
  if (is.null(dir)) {
    dirs <- session_error_dirs(config)
    dir <- if (length(dirs) >= index) dirs[index] else NA_character_
  }
  empty <- data.table::data.table(exp = character(), step = character(), message = character(),
                                  hint = character(), time = as.POSIXct(character()), file = character())
  if (is.na(dir) || !dir.exists(dir)) return(empty)
  files <- list.files(dir, "\\.rds$", full.names = TRUE)
  if (!length(files)) return(empty)
  data.table::rbindlist(lapply(files, function(f) {
    raw <- tryCatch(readRDS(f), error = function(e) NA_character_)
    raw <- if (is.character(raw) && length(raw) == 1) raw else NA_character_
    step <- if (!is.na(raw) && regexpr("\\(\\s*[^)]*?\\s*\\)", raw) > 0)
      gsub("^\\(\\s*|\\s*\\)$", "", regmatches(raw, regexpr("\\(\\s*[^)]*?\\s*\\)", raw))) else NA_character_
    rest <- if (!is.na(raw)) sub("^[^)]*\\)", "", raw) else ""
    data.table::data.table(exp = sub("\\.rds$", "", basename(f)), step = step,
                           message = error_one_line(sub("\nLikely cause:.*", "", rest)),
                           hint = error_hint(rest), time = file.mtime(f), file = f)
  }))[order(time)]
}

#' One-line error summary for the checklist header
#' @param config the mNGSp config object
#' @return character(1), "" when there are no current errors
#' @noRd
pipeline_error_summary_line <- function(config) {
  if (is.null(config$error_dir)) return("")
  e <- tryCatch(pipeline_errors(config), error = function(err) NULL)
  crashes <- if (!is.null(config$session_dir)) length(list.files(config$session_dir, "^crash_.*\\.rds$")) else 0
  if (is.null(e) || (!nrow(e) && !crashes)) return("")
  items <- if (nrow(e)) paste0(e$exp, " (", e$step, ")") else character()
  shown <- paste(utils::head(items, 4), collapse = ", ")
  more <- if (length(items) > 4) sprintf(", +%d more", length(items) - 4) else ""
  crash_txt <- if (crashes) sprintf("; %d stage crash record(s) in session_logs", crashes) else ""
  sprintf("ERRORS this session: %d%s%s%s -- pipeline_errors(config) for details, clear_pipeline_error() to retry",
          nrow(e), if (nrow(e)) paste0(": ", shown) else "", more, crash_txt)
}

#' Remove an error flag so a running pipeline retries the study
#'
#' The record is moved to \code{<error_dir>/resolved/} (not deleted) and a \code{"resolved"} event is emitted.
#' Stage loops check the flag file before every study, so a running \code{run_pipeline()} picks the study up on its
#' next round, and the checklist stops showing it.
#' @param config the mNGSp config object
#' @param exp character, experiment/study ids (file names without .rds) to clear
#' @param index integer, session index (newest first) when \code{config$error_dir} is not set
#' @param note character, reason, stored in the event
#' @return invisible character vector of moved files
#' @export
clear_pipeline_error <- function(config, exp, index = 1, note = "") {
  dir <- config$error_dir %||% session_error_dirs(config)[index]
  files <- file.path(dir, paste0(exp, ".rds"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("No error record for: ", paste(basename(missing), collapse = ", "), " in ", dir)
  dest <- file.path(dir, "resolved"); dir.create(dest, showWarnings = FALSE)
  moved <- file.path(dest, paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"), basename(files)))
  ok <- file.rename(files, moved)
  if (!all(ok)) stop("Could not move: ", paste(basename(files[!ok]), collapse = ", "))
  cfg <- config; cfg$error_dir <- dir
  for (i in seq_along(exp)) error_event_emit(cfg, list(event = "resolved", exp = exp[i], step = "",
                                                         message = note, hint = "", error_file = moved[i]))
  invisible(moved)
}

#' Read the error event stream
#' @param config the mNGSp config object (or a list with \code{project})
#' @param since POSIXct or character time; only later events are returned (default: all)
#' @param n integer, return only the last n events
#' @return data.table of events (columns as written by the pipeline)
#' @export
pipeline_error_events <- function(config, since = NULL, n = Inf) {
  path <- error_events_path(config)
  if (!file.exists(path)) return(data.table::data.table())
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (is.finite(n)) lines <- utils::tail(lines, n)
  ev <- data.table::rbindlist(lapply(lines, function(l) tryCatch(jsonlite::fromJSON(l), error = function(e) NULL)),
                              fill = TRUE)
  if (!is.null(since) && nrow(ev)) ev <- ev[as.POSIXct(substr(time, 1, 19), format = "%Y-%m-%dT%H:%M:%S") >= as.POSIXct(since)]
  ev
}
