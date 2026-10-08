# "What is going on, right now" helpers -- built directly from this
# session's own repeated need to manually `ps aux | grep`,
# tail session_logs/error_logs by hand, and re-derive the same picture
# from scratch after every context reset while debugging a live hang
# (2026-10-07/08, PRJNA637713). The existing pshift_triage()
# (R/pipeline_qc_report.R) already established the right PATTERN for
# this (cheap, read-only, built purely from precomputed files) for one
# specific stage; everything here generalizes that same philosophy to
# "what's going on" as a whole, not just pshift.

#' Which studies (if any) currently have massiveNGSpipe-related processes running, system-wide
#'
#' Broader than \code{\link{sample_currently_processing}} (which only
#' checks for \code{fastp}/\code{STAR} processes mentioning ONE
#' specific experiment): scans every process on the machine for (a) a
#' recognizable study-accession-shaped token (\code{PRJNA}/\code{PRJEB}/
#' \code{PRJDB}/\code{GSE}/\code{ERP}/\code{SRP}/\code{DRP}-prefixed)
#' anywhere in its command line, and (b) \code{fastp}/\code{STAR}/
#' \code{massiveNGSpipe} activity generally, even when no accession is
#' visible in that specific process's own command line (e.g. an outer
#' \code{Rscript} wrapper or RStudio's own \code{sourceWithProgress}
#' launcher). Read-only: never signals, kills, or otherwise touches
#' anything it finds.
#' @return data.table(pid, accession, cmd_snippet) -- one row per
#' matching process. \code{accession} is \code{NA} when a process
#' looks massiveNGSpipe-related but no accession-shaped token was
#' found in its own command line (common for outer wrapper processes).
#' A 0-row data.table (not NULL) when nothing matches.
#' @export
studies_currently_processing <- function() {
  empty <- data.table::data.table(pid = integer(), accession = character(), cmd_snippet = character())
  procs <- tryCatch(system2("ps", c("-eo", "pid,args"), stdout = TRUE, stderr = FALSE),
                    warning = function(w) character(), error = function(e) character())
  if (length(procs) <= 1) return(empty)
  procs <- procs[-1] # drop the header line ("PID COMMAND")

  accession_pattern <- "(PRJNA|PRJEB|PRJDB|GSE|ERP|SRP|DRP)[0-9]+"
  relevant <- grepl("massiveNGSpipe|fastp|\\bSTAR\\b", procs, perl = TRUE) |
    grepl(accession_pattern, procs, perl = TRUE)
  procs <- procs[relevant]
  if (length(procs) == 0) return(empty)

  pid <- as.integer(trimws(sub("^\\s*([0-9]+)\\s.*", "\\1", procs)))
  m <- regexpr(accession_pattern, procs, perl = TRUE)
  # regmatches() on a regexpr()-derived m silently DROPS non-matching
  # elements rather than returning "" placeholders for them, so its
  # result can be shorter than procs -- indexing by has_match keeps
  # NA aligned to the right row instead of (as bare vapply(raw, ...)
  # did, before this fix) letting data.table() recycle a too-short
  # accession vector across every row.
  has_match <- m > 0
  accession <- rep(NA_character_, length(procs))
  accession[has_match] <- regmatches(procs, m)
  snippet <- substr(trimws(sub("^\\s*[0-9]+\\s*", "", procs)), 1, 120)

  data.table::data.table(pid = pid, accession = accession, cmd_snippet = snippet)
}

#' Is one session's checklist stale -- no update in a while, with no
#' final status label written?
#'
#' A clean end to \code{\link{run_pipeline}()} always labels
#' \code{checklist.txt}'s title line \code{"(done/stopped gracefully/
#' aborted after N.N hours)"} via its own \code{on.exit()} handler --
#' but a hard kill (\code{kill -9}, an OOM-killer, a container crash)
#' gives that handler no chance to run at all, so the file just goes
#' stale (its own mtime stops advancing) with NO label, which is the
#' only signal available. Confirmed live, 2026-10-07/08: exactly this
#' pattern (checklist present, unlabeled, mtime far in the past) was
#' the only externally-visible symptom of a ~13 hour nested-fork/
#' thread-exhaustion hang (PRJNA637713) -- discovering it took manual
#' process-state digging (\code{ps}, \code{/proc/<pid>/wchan}) because
#' nothing surfaced it automatically. This function is that automatic
#' surfacing.
#' @inheritParams session_checklist_path
#' @param stale_after_mins numeric, default 30. How long since the
#' checklist's own last modification before an UNLABELED title counts
#' as "possibly hung" rather than "still legitimately running".
#' @return a list: \code{path} (character), \code{exists} (logical),
#' \code{title} (character or NA), \code{labeled} (logical -- did the
#' title line get a final status note at all), \code{mins_since_update}
#' (numeric or NA), \code{status} one of \code{"no_checklist"},
#' \code{"done"}, \code{"stopped_gracefully"}, \code{"aborted"},
#' \code{"running"}, or \code{"possibly_hung"}.
#' @export
checklist_health <- function(config, index = 1, stale_after_mins = 30) {
  path <- if (missing(index) && !is.null(config$session_dir)) {
    file.path(config$session_dir, "checklist.txt")
  } else {
    # A project with no session_logs at all yet (never run) must read
    # as "no_checklist", not propagate session_checklist_path()'s own
    # "no existing sessions" error -- confirmed needed directly,
    # pipeline_status_all() hits exactly this for a fresh project dir.
    tryCatch(session_checklist_path(config, index), error = function(e) NA_character_)
  }
  if (is.na(path) || !file.exists(path))
    return(list(path = path, exists = FALSE, title = NA_character_, labeled = FALSE,
               mins_since_update = NA_real_, status = "no_checklist"))

  title <- readLines(path, n = 1, warn = FALSE)
  mins_since_update <- as.numeric(difftime(Sys.time(), file.info(path)$mtime, units = "mins"))

  status <- if (grepl("\\(done after", title)) "done"
  else if (grepl("\\(stopped gracefully after", title)) "stopped_gracefully"
  else if (grepl("\\(aborted after", title)) "aborted"
  else if (mins_since_update > stale_after_mins) "possibly_hung"
  else "running"

  list(path = path, exists = TRUE, title = title,
      labeled = status %in% c("done", "stopped_gracefully", "aborted"),
      mins_since_update = round(mins_since_update, 1), status = status)
}

#' Scan a set of project directories for current pipeline status, in one call
#'
#' A thin wrapper combining \code{\link{checklist_health}} and
#' \code{\link{studies_currently_processing}} across every project dir
#' given -- e.g. the full online/local x Ribo/RNA/disome/SSU set
#' documented in this server's own AGENTS.md, or any subset relevant to
#' the task at hand. Deliberately takes \code{project_dirs} as a plain
#' argument rather than guessing a default: which project dirs exist,
#' and which matter for a given task, varies by environment and by
#' what's actually being worked on -- hardcoding one server's own path
#' list here would silently stop generalizing the moment either
#' changes.
#' @param project_dirs character vector of project directories (each
#' one \code{pipeline_config()}'s own \code{project_dir}, e.g.
#' \code{"~/livemount/Bio_data/NGS_pipeline"}).
#' @return data.table(project_dir, status, title, mins_since_update,
#' n_active_processes) -- one row per project dir. A project dir with
#' no \code{session_logs} at all (never run) gets
#' \code{status = "no_checklist"}.
#' @export
pipeline_status_all <- function(project_dirs) {
  active <- studies_currently_processing()
  rows <- lapply(project_dirs, function(pd) {
    health <- tryCatch(checklist_health(list(project = pd, session_dir = NULL)),
                       error = function(e) list(title = NA_character_, mins_since_update = NA_real_,
                                                status = paste0("error: ", conditionMessage(e))))
    n_active <- sum(grepl(basename(pd), active$cmd_snippet, fixed = TRUE)) +
      sum(vapply(active$accession, function(a) !is.na(a) && grepl(a, pd, fixed = TRUE), logical(1)))
    data.table::data.table(project_dir = pd, status = health$status, title = health$title,
                           mins_since_update = health$mins_since_update, n_active_processes = n_active)
  })
  data.table::rbindlist(rows, fill = TRUE)
}

#' Regex matching massiveNGSpipe's own known install/lazy-load race error messages
#'
#' Shared by \code{\link{is_install_race_error}} (condition objects,
#' R/pipeline_async_exec.R) and \code{\link{summarize_last_errors}}
#' (plain strings already extracted from \code{error_logs}) -- same
#' pattern, two different input shapes.
#' @return character, a regex (for use with \code{perl = TRUE})
#' @noRd
install_race_error_regex <- function() {
  paste0("massiveNGSpipe\\.rdb['\"]? is corrupt|",
        "read failed on .*massiveNGSpipe\\.rdb|",
        "there is no package called .massiveNGSpipe.|",
        "cannot open file '[^']*massiveNGSpipe\\.rdb'")
}

#' Friendlier wrapper around \code{\link{last_session_errors}}: flags known-transient errors
#'
#' \code{last_session_errors()} returns raw error strings with no
#' indication of which ones are massiveNGSpipe's own known,
#' self-healing install/lazy-load race (already retried automatically
#' by \code{\link{run_experiment_subprocess}} where it occurs on that
#' path -- but older error records, or errors from a path that doesn't
#' go through that retry, still show up here looking identical to a
#' real failure). This separates the two so a real error doesn't get
#' lost in a sea of transient noise, without hiding either.
#' @inheritParams last_session_errors
#' @return data.table(error_key, likely_transient, message).
#' \code{error_key} is exactly \code{names(last_session_errors(...))}
#' -- NOT a clean experiment name on its own (it's
#' \code{"/<session-date>/<exp>.rds"} with \code{"error_logs"} and
#' \code{".rds"} stripped, see \code{\link{last_session_errors}}'s own
#' implementation), but still useful for matching against an
#' experiment name via \code{grepl(exp, error_key, fixed = TRUE)}.
#' Deliberately not named \code{key} -- that's a reserved parameter
#' name in \code{data.table()}'s own constructor (sets the table's
#' sort key rather than becoming a literal column), confirmed directly
#' when this collided during testing.
#' @export
summarize_last_errors <- function(config, index = 1, regex = NULL) {
  errors <- last_session_errors(config, index = index, regex = regex)
  if (length(errors) == 0)
    return(data.table::data.table(error_key = character(), likely_transient = logical(), message = character()))
  data.table::data.table(error_key = names(errors),
                         likely_transient = grepl(install_race_error_regex(), errors, perl = TRUE),
                         message = unname(errors))
}

#' One-call "what's going on" snapshot, for a human or an AI agent with no prior context
#'
#' Combines every other function in this file plus the live
#' \code{config} itself into one read-only snapshot -- built directly
#' from this session's own experience of a fresh agent (no memory of
#' prior sessions) needing to manually re-derive this exact picture,
#' piece by piece, via a dozen separate commands, every time it picked
#' the investigation back up. Never runs, mutates, or blocks anything.
#' @param config the mNGSp config object
#' @return a list: \code{active_processes} (see
#' \code{\link{studies_currently_processing}}), \code{checklist} (see
#' \code{\link{checklist_health}}), \code{recent_errors} (see
#' \code{\link{summarize_last_errors}}), \code{threads}
#' (\code{config$threads}), \code{thread_type} (character, the
#' configured BPPARAM constructor's name), \code{mode}, \code{project}.
#' @export
ai_context_snapshot <- function(config) {
  list(
    active_processes = tryCatch(studies_currently_processing(),
                                error = function(e) data.table::data.table()),
    checklist = tryCatch(checklist_health(config), error = function(e) NULL),
    recent_errors = tryCatch(summarize_last_errors(config), error = function(e) data.table::data.table()),
    threads = config$threads,
    thread_type = tryCatch(get_fun_name(config$thread_type, "BiocParallel"), error = function(e) NA_character_),
    mode = config$mode,
    project = config$project
  )
}

#' Generalized, read-only triage for one experiment/step -- the
#' \code{\link{pshift_triage}} philosophy, extended to any stage
#'
#' \code{pshift_triage()} (R/pipeline_qc_report.R) is cheap, read-only,
#' and built purely from precomputed files -- exactly the right pattern
#' for "what's wrong with this experiment", but scoped to the pshift
#' stage specifically. This is the same idea, generalized: whichever
#' step you ask about, it reads that step's own existing flags/markers/
#' errors/process-activity instead of needing a bespoke diagnostic
#' written for every stage individually.
#' @param exp character, experiment name
#' @param step character, a name in \code{config$flag} (e.g. \code{"trim"},
#' \code{"aligned"}, \code{"valid_pshift"})
#' @param config the mNGSp config object
#' @return a list: \code{step_done} (logical), \code{n_samples_done}/
#' \code{samples_done} (only meaningful for steps with per-sample
#' markers, see \code{\link{sample_flag_dir}}'s own roxygen for which
#' those are -- \code{NA}/\code{character()} otherwise), \code{currently_processing}
#' (logical, via \code{\link{studies_currently_processing}}),
#' \code{last_touch} (POSIXct, via \code{\link{last_real_touch}}),
#' \code{recent_errors} (named character, this exp's own entries from
#' \code{\link{last_session_errors}})
#' @export
stage_triage <- function(exp, step, config) {
  active <- tryCatch(studies_currently_processing(), error = function(e) data.table::data.table())
  errs <- tryCatch(last_session_errors(config), error = function(e) character())
  list(
    step_done = tryCatch(step_is_done(config, step, exp), error = function(e) NA),
    n_samples_done = tryCatch(n_samples_done(config, step, exp), error = function(e) NA_integer_),
    samples_done = tryCatch(samples_done(config, step, exp), error = function(e) character()),
    currently_processing = if (nrow(active) == 0) FALSE else
      any(grepl(exp, active$cmd_snippet, fixed = TRUE)) ||
        any(!is.na(active$accession) & mapply(grepl, active$accession, exp, fixed = TRUE)),
    last_touch = tryCatch(last_real_touch(exp, config), error = function(e) as.POSIXct(NA)),
    recent_errors = errs[grepl(exp, names(errs), fixed = TRUE)]
  )
}
