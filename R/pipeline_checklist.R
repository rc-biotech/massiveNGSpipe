#' Map a stage-group name to its sample-marker step id, if it has one
#'
#' Only stages with a plain in-process per-sample loop we control directly
#' get sample-level detail (see R/pipeline_sample_flags.R); everything else
#' (stages that hand off to ORFik's own internal BiocParallel dispatch, or
#' process a whole folder in one call) only ever gets study-level counts.
#' @param stage_name character, e.g. "pipe_trim_collapse"
#' @return character step id (e.g. "trim") or NA_character_
stage_marker_step <- function(stage_name) {
  switch(stage_name,
        pipe_fetch = "fetch",
        pipe_trim_collapse = "trim",
        pipe_align_clean = "aligned",
        pipe_exp_ofst = "ofst",
        # pipe_convert runs covRLE then bigwig; a single marker can only
        # track one, so bigwig is chosen (the later phase) -- during the
        # covRLE phase this stage just shows "queued" instead of live
        # per-sample progress, rather than showing stale/misleadingly-
        # complete covRLE counts once the bigwig phase has actually begun.
        pipe_convert = "bigwig",
        NA_character_)
}

#' One stage's live processing-speed label for the checklist, or NA
#'
#' Dispatches to whichever rate source (R/pipeline_checklist_rates.R)
#' applies to this stage: fetch's own growing download file (MB/s),
#' STAR's own Log.progress.out for align (M reads/hr, already computed
#' by STAR itself), a generic samples/hour for the other per-sample
#' marker stages, or studies/hour for stages with only study-level
#' counts. Only ever called for a "running" stage; every other state
#' (done/queued) always gets NA, so nothing is rendered for them.
#' @param stage_name character
#' @param marker_step character or NA, from \code{\link{stage_marker_step}}
#' @param state character, this stage's already-computed state
#' @param config the mNGSp config object
#' @param active_experiment character or NA
#' @param active_done integer or NA, samples done for active_experiment
#' @param n_done integer, studies done for this whole stage
#' @param pipelines the pipelines list
#' @return character label (e.g. "42.3 MB/s"), or NA_character_
#' @noRd
stage_progress_rate_label <- function(stage_name, marker_step, state, config,
                                      active_experiment, active_done, n_done,
                                      pipelines) {
  if (state != "running") return(NA_character_)

  if (stage_name == "pipe_fetch" && !is.na(active_experiment)) {
    run <- active_run_id(pipelines, config, "fetch", active_experiment)
    conf <- experiment_conf(pipelines, active_experiment)
    if (is.na(run) || is.null(conf)) return(NA_character_)
    rate <- fetch_progress_rate(conf["fastq"], run)
    if (is.na(rate)) return(NA_character_)
    return(sprintf("%.1f MB/s", rate))
  }

  if (stage_name == "pipe_align_clean" && !is.na(active_experiment)) {
    run <- active_run_id(pipelines, config, "aligned", active_experiment)
    conf <- experiment_conf(pipelines, active_experiment)
    if (is.na(run) || is.null(conf)) return(NA_character_)
    rate <- align_progress_rate(conf["bam"], run)
    if (is.na(rate)) return(NA_character_)
    return(sprintf("%.1f M reads/hr", rate))
  }

  if (!is.na(marker_step) && !is.na(active_experiment)) {
    rate <- generic_progress_rate(stage_name, active_experiment, active_done)
    if (is.na(rate)) return(NA_character_)
    return(sprintf("%.1f samples/hr", rate))
  }

  if (is.na(marker_step)) {
    rate <- generic_progress_rate(stage_name, NA_character_, n_done)
    if (is.na(rate)) return(NA_character_)
    return(sprintf("%.1f studies/hr", rate))
  }

  NA_character_
}

#' Format elapsed time since a start time as "N.N hours"
#' @param start POSIXct, e.g. config$init_time
#' @param end POSIXct, default Sys.time()
#' @return character, e.g. "5.2 hours"
#' @noRd
format_elapsed_hours <- function(start, end = Sys.time()) {
  paste0(round(as.numeric(difftime(end, start, units = "hours")), 1), " hours")
}

#' Count of experiments fetched but not yet past the fetch step
#'
#' Mirrors \code{pipe_fetch()}'s own gating computation exactly (factored
#' out here so \code{pipeline_checklist()} can display the same live
#' number \code{pipe_fetch()} uses to decide whether to pause -- see
#' \code{config$max_unprocessed_downloads}). \code{progress_report()}
#' itself is fairly verbose (prints a status table, queries system usage
#' again); suppressed here since this helper's only job is to return a
#' number, not to print -- in production this output is already silently
#' buffered away by BiocParallel's own log capture regardless (see
#' \code{\link{pipeline_checklist}}'s own docs), so suppressing it here
#' changes nothing observable for \code{pipe_fetch()}'s existing caller.
#' @param pipelines the pipelines list
#' @param config the mNGSp config object
#' @return integer, count of experiments whose progress is exactly at
#' the "fetch" step (fetched, but the next step -- normally "trim" --
#' isn't done yet). 0 if this config has no "fetch" step at all (e.g.
#' \code{mode = "local"}, where \code{pipe_fetch()} is never even
#' included in \code{config$pipeline_steps}).
#' @noRd
unprocessed_downloads_count <- function(pipelines, config) {
  if (!("fetch" %in% names(config$flag))) return(0L)
  flag_step <- which(names(config$flag) == "fetch")
  invisible(capture.output(suppressMessages({
    progress <- progress_report(pipelines, config, show_status_per_exp = FALSE,
                                show_done = FALSE, return_progress_vector = TRUE,
                                system_usage_stats = FALSE)
  })))
  sum(progress == flag_step)
}

#' Total sample count per experiment, from the live pipelines object
#'
#' Not derivable from flags alone -- needs the actual per-run metadata rows.
#' @param pipelines the pipelines list (see pipeline_init_all())
#' @return named integer vector, names are experiment ids
experiment_sample_counts <- function(pipelines) {
  counts <- list()
  for (p in pipelines) {
    for (organism in names(p$organisms)) {
      exp <- p$organisms[[organism]]$conf["exp"]
      counts[[exp]] <- nrow(p$study[ScientificName == organism])
    }
  }
  unlist(counts)
}

#' Base directory for this run's live checklist/console-log artifacts
#'
#' Session-scoped whenever \code{config$session_dir} is set (i.e. running
#' through \code{\link{run_pipeline}}, which sets it once per call via
#' \code{run_pipeline_set_up_session()} before dispatching any stage-group
#' work -- so every downstream function already receives a \code{config}
#' with it populated, no extra plumbing needed). This matters because two
#' concurrent \code{run_pipeline()} calls (e.g. processing two different
#' subsets of \code{pipelines} at once) would otherwise both write
#' \code{checklist.txt} and the same per-experiment console logs to the
#' exact same project-level path, silently clobbering each other -- the
#' whole reason \code{session_logs/<init_time>/} exists as a unique
#' per-call directory in the first place (see \code{\link{session_info_table}}
#' for browsing/swapping between past sessions, and \code{mNGSp_app()} for
#' a live Shiny view of one selected session at a time).
#'
#' Falls back to a fixed project-level location when called outside a real
#' session (e.g. calling a \code{pipeline_*()} function directly,
#' interactively, without going through \code{run_pipeline()}) -- that
#' fallback path does not change between calls, so it is not
#' concurrency-safe and should not be relied on for simultaneous runs.
#' @param config the mNGSp config object
#' @return character, a directory path (not guaranteed to exist yet)
pipeline_log_base <- function(config) {
  if (!is.null(config$session_dir)) config$session_dir
  else file.path(config$project, "log_pipeline")
}

#' Path to the live checklist snapshot file for one run_pipeline() session
#' @inheritParams pipeline_log_base
#' @return character, file.path(pipeline_log_base(config), "checklist.txt")
checklist_path <- function(config) file.path(pipeline_log_base(config), "checklist.txt")

#' List a project's session_logs session directories, newest first
#'
#' Same pattern as \code{\link{session_error_dirs}} (R/pipeline_logs.R),
#' applied to \code{session_logs} instead of \code{error_logs}.
#' @param config the mNGSp config object
#' @return character vector of directory paths
#' @noRd
session_log_dirs <- function(config) {
  sort(list.dirs(file.path(config$project, "session_logs"), recursive = FALSE), decreasing = TRUE)
}

#' Resolve one session's checklist.txt by index, newest first
#'
#' Deliberately derived from the \code{session_logs} directory listing,
#' not \code{config$session_dir} -- a fresh \code{config} object built
#' in a SEPARATE terminal/session (the whole point of
#' \code{\link{watch_pipeline_checklist}}, watching a run happening
#' elsewhere) never has \code{session_dir} set, so \code{index = 1}
#' (the newest \code{session_logs} entry) is what actually resolves to
#' a currently-running session's own live checklist -- not
#' \code{\link{checklist_path}}'s project-level fallback, which would
#' silently point at the wrong (and likely stale) file in that exact
#' "watch it from another terminal" scenario.
#' @param config the mNGSp config object
#' @param index integer, default 1 (newest). Same convention as
#' \code{\link{last_session_errors}}'s own \code{index} argument.
#' @return character, path to that session's checklist.txt
#' @noRd
session_checklist_path <- function(config, index = 1) {
  dirs <- session_log_dirs(config)
  if (length(dirs) < index)
    stop("You selected session ", index, ", but there ", if (length(dirs) == 1) "is" else "are",
        " only ", length(dirs), " existing session(s)")
  file.path(dirs[index], "checklist.txt")
}

#' Nextflow-style stage checklist
#'
#' Read-only: reports on experiment-level flags (existing) and sample-level
#' markers (new, only available for stages with a plain per-sample loop --
#' see stage_marker_step()). Never downloads, mutates config, or blocks.
#'
#' Always writes a plain-text snapshot to checklist_path(config) via a direct
#' file write (not message()/cat() to the R stdout/message streams). This
#' matters in production: run_pipeline()'s outer bplapply always runs with
#' config$parallel_conf's logdir set, and BiocParallel's own log=TRUE/logdir=
#' capture buffers *all* stdout/message output from a task and only flushes
#' it to disk once that whole task (one stage-group, potentially hours)
#' completes -- verified directly (SerialParam, MulticoreParam, both via
#' plain Rscript). message()/cat() calls made from inside parallel_wrap()
#' are therefore invisible, live, both on the console and in that per-task
#' log file, no matter which of the two is used. A direct file() write does
#' not go through that capture at all, so it is the only thing that actually
#' updates live; message() below is kept only as a best-effort extra for
#' contexts where parallel_wrap() is called outside of that bplapply
#' wrapping (e.g. directly, interactively).
#'
#' @param pipelines the pipelines list
#' @param config the mNGSp config object
#' @param print logical, default TRUE. If TRUE, also render via message()
#' (see Details above for why this alone is not enough in production).
#' @param run_status character, default NULL. Parenthetical shown after
#' the title line's timestamp. NULL means "still running": auto-computed
#' as \code{"(running for N.N hours)"} from \code{config$init_time}
#' (silently omitted if \code{init_time} isn't set, e.g. calling this
#' outside a real \code{run_pipeline()} session). A caller passes an
#' explicit string -- e.g. \code{"done after 5.2 hours"},
#' \code{"aborted after 5.2 hours"} -- for a final, one-off status
#' (see \code{\link{run_pipeline}}'s own on.exit handler).
#' @return invisible(data.table) with columns: stage, done, total, state
#' ("done"/"running"/"queued"), active_experiment, active_done, active_total
#' (the latter three NA when no marker evidence is available/applicable)
pipeline_checklist <- function(pipelines, config, print = TRUE, run_status = NULL) {
  exps <- pipelines_names(pipelines)
  sample_totals <- experiment_sample_counts(pipelines)

  rows <- lapply(names(config$flag_steps), function(stage_name) {
    flags <- config$flag_steps[[stage_name]]
    done_vec <- vapply(exps, function(e) all(vapply(flags, step_is_done, logical(1),
                                                     config = config, experiment = e)),
                       logical(1))
    n_done <- sum(done_vec)
    n_total <- length(exps)

    marker_step <- stage_marker_step(stage_name)
    active_experiment <- NA_character_
    active_done <- NA_integer_
    active_total <- NA_integer_
    state <- if (n_done == n_total) "done" else "queued"

    if (n_done < n_total && !is.na(marker_step)) {
      candidates <- exps[!done_vec]
      started <- vapply(candidates, function(e) n_samples_done(config, marker_step, e) > 0, logical(1))
      if (any(started)) {
        active_experiment <- candidates[started][1]
        active_done <- n_samples_done(config, marker_step, active_experiment)
        active_total <- unname(sample_totals[active_experiment])
        state <- "running"
      }
    } else if (is.na(marker_step) && n_done > 0 && n_done < n_total) {
      # No per-sample detail is possible for these (hand off to ORFik's
      # own internal BiocParallel dispatch or a whole-folder call), but
      # at least one study already finished this stage, so it HAS been
      # actively cycling through studies -- "queued" would otherwise
      # misleadingly suggest nothing has started yet.
      state <- "running"
    }

    rate_label <- stage_progress_rate_label(stage_name, marker_step, state, config,
                                            active_experiment, active_done, n_done,
                                            pipelines)

    data.table::data.table(stage = stage_name, done = n_done, total = n_total,
                           rate_label = rate_label,
                           state = state, active_experiment = active_experiment,
                           active_done = active_done, active_total = active_total)
  })
  tab <- data.table::rbindlist(rows)
  txt <- format_checklist(tab)
  usage <- format_system_usage_line(pipelines, config)

  # Fixed-height header, always exactly 6 lines before the stage table
  # (title / usage / drive-cap-note-or-blank / backlog-cap-note-or-blank /
  # blank / blank) whether or not either cap note has anything to say --
  # so the stage table always starts on the same line number and a
  # live-redrawing viewer (or someone just watching the plain file)
  # doesn't get thrown off by the line count shifting as either note
  # appears/disappears between polls.
  status_note <- if (!is.null(run_status)) {
    paste0(" (", run_status, ")")
  } else if (!is.null(config$init_time)) {
    paste0(" (running for ", format_elapsed_hours(config$init_time), ")")
  } else ""

  path <- checklist_path(config)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  header <- c(sprintf("Pipeline status as of %s%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), status_note),
             usage["line"], usage["drive_cap_note"], usage["backlog_cap_note"], "", "")
  cat(paste(c(header, txt), collapse = "\n"), "\n", sep = "", file = path)

  if (print) message(paste(c(usage["line"], usage["drive_cap_note"], usage["backlog_cap_note"],
                             "", "", txt), collapse = "\n"))
  invisible(tab)
}

#' One-line system usage summary for the checklist, plus two separate
#' fixed-slot cap-warning notes
#'
#' The usage line matches \code{ORFik::get_system_usage(one_liner = TRUE)}'s
#' own format (\code{"CPU (x%), Memory (y%), Drive <drive> (z%)"} -- built
#' directly from the list-mode return value here rather than capturing
#' that function's own \code{cat()} side effect, so the same call also
#' gives the numeric drive percentage needed for the cap check below
#' without a second \code{df}/\code{top} invocation).
#'
#' \code{pipe_fetch()} pauses new downloads for either of two independent
#' reasons -- drive usage at/above \code{stop_downloading_new_data_at_drive_usage},
#' or too many fetched-but-unprocessed experiments
#' (\code{max_unprocessed_downloads}, see \code{\link{unprocessed_downloads_count}})
#' -- and both are reported here as their own separate, always-reserved
#' line (not one combined note): appending one note's text length onto
#' the other's would make that line's length vary depending on which
#' cap(s) are active, which wraps unpredictably in a narrow terminal/UI
#' and pushes every line below it out of position. Same reasoning as
#' why the usage line and cap notes are already split from each other.
#' @param pipelines the pipelines list (needed for the backlog cap note;
#' see \code{\link{unprocessed_downloads_count}})
#' @param config the mNGSp config object
#' @return named character vector of length 3: \code{line} (the usage
#' summary), \code{drive_cap_note} (the drive-usage warning, or
#' \code{""} when below that cap), \code{backlog_cap_note} (the
#' unprocessed-downloads warning, or \code{""} when below that cap)
format_system_usage_line <- function(pipelines, config) {
  drive <- detect_drive(path.expand(config$config["ref"]))
  usage <- get_system_usage(drive)
  line <- paste0("CPU (", usage$CPU_Usage_Percent, "%),",
                " Memory (", usage$Memory_Usage_Percent, "%),",
                " Drive ", usage$Drive, " (", usage$Drive_Usage_Percent, "%)")

  cap <- config$stop_downloading_new_data_at_drive_usage
  drive_pct <- suppressWarnings(as.numeric(gsub("%", "", usage$Drive_Usage_Percent)))
  drive_cap_note <- ""
  if (!is.null(cap) && !is.na(drive_pct) && drive_pct >= cap) {
    drive_cap_note <- paste0("[DRIVE AT/ABOVE ", cap, "% CAP -- pipe_fetch() is pausing new downloads]")
  }

  backlog_cap <- config$max_unprocessed_downloads
  unprocessed <- unprocessed_downloads_count(pipelines, config)
  backlog_cap_note <- ""
  if (!is.null(backlog_cap) && !is.na(unprocessed) && unprocessed > backlog_cap) {
    backlog_cap_note <- paste0("[UNPROCESSED DOWNLOADS ", unprocessed, " > ", backlog_cap,
                               " CAP -- pipe_fetch() is pausing new downloads]")
  }
  c(line = line, drive_cap_note = drive_cap_note, backlog_cap_note = backlog_cap_note)
}

#' Live-watch the pipeline checklist in place, like a download progress bar
#'
#' Continuously redraws one session's checklist content in the terminal
#' in place -- like \code{curl}'s progress bar, or \code{docker compose
#' up}'s multi-line status block -- instead of printing a new block
#' underneath on every refresh. Purely a read-only viewer: it never
#' runs, mutates, or blocks the pipeline itself. Meant to be run in a
#' second terminal/session alongside a real \code{\link{run_pipeline}}
#' call happening elsewhere (in another terminal, or in the
#' background).
#'
#' Resolves which checklist to watch as follows: if \code{config$session_dir}
#' is set AND \code{index} was left at its default (i.e. not explicitly
#' requested), that session is used directly -- it's the authoritative
#' "this specific session" pointer when called from within an
#' already-running session, and strictly more correct than "newest by
#' directory listing" whenever a second, newer \code{run_pipeline()}
#' call has started concurrently elsewhere (which would otherwise make
#' \code{index = 1} resolve to that OTHER session). Otherwise -- the
#' normal case of a fresh \code{config} built in a separate
#' terminal/session, which never has \code{session_dir} set, the main
#' use case this function exists for -- the checklist is resolved via
#' \code{index} against the \code{session_logs} directory listing
#' (newest first), so \code{index = 1} means "whichever session is
#' newest", the currently-running one in the normal live-watch case.
#' An explicitly-passed \code{index} always wins over
#' \code{session_dir}, even when the latter is set, since that's a
#' clear request to inspect a specific past session.
#'
#' Uses ANSI cursor-movement escape codes (move up N lines, clear to end of
#' screen, redraw), so it needs a real terminal emulator -- an \code{Rscript}
#' or R session run from an actual shell. RStudio's own Console pane does
#' NOT support cursor repositioning (only plain text/color codes), so this
#' will not redraw in place there, only append -- use
#' \code{open_in_new_text_window = TRUE} in that case (or from a plain
#' shell with no RStudio at all, \code{watch -n 2 cat <path>} does the
#' same job directly, no R needed).
#'
#' @param config the mNGSp config object
#' @param interval numeric, seconds between redraws, default 2
#' @param index integer, default 1 (newest session). Same convention as
#' \code{\link{last_session_errors}}'s own \code{index} argument --
#' pass e.g. \code{2} to watch/inspect the previous session's final
#' checklist instead of the current/newest one.
#' @param open_in_new_text_window logical, default FALSE. When TRUE,
#' opens a real RStudio Terminal tab (via \code{rstudioapi::terminalCreate()})
#' running \code{watch -n <interval> cat <path>} there instead of
#' redrawing in THIS console -- the Terminal pane is a real terminal
#' emulator (unlike the Console), so \code{\r}/ANSI redraw-in-place
#' actually works. Requires an active RStudio session; errors
#' otherwise. Returns immediately (does not block this R session), and
#' does not use \code{interval} for its own redraw loop beyond passing
#' it through to \code{watch}.
#' @return invisible(NULL), or (when \code{open_in_new_text_window = TRUE})
#' invisibly the new terminal's id. The non-window form runs until
#' interrupted (Ctrl+C, or Esc in RStudio -- though see the terminal
#' note above).
#' @export
watch_pipeline_checklist <- function(config, interval = 2, index = 1, open_in_new_text_window = FALSE) {
  # config$session_dir -- when set AND index wasn't explicitly
  # requested -- takes priority over the index/session_logs listing:
  # it's the authoritative "this specific session" pointer when
  # called from within an already-running session, and strictly more
  # correct than "newest by directory listing" in the edge case where
  # a second, newer run_pipeline() call started concurrently elsewhere
  # (that would otherwise make index = 1 resolve to the OTHER session,
  # not the one this config actually belongs to). An explicit index
  # always wins, even when session_dir is set, since that's a clear
  # request to inspect a specific past session rather than "the
  # current one, however that's best determined".
  path <- if (missing(index) && !is.null(config$session_dir)) {
    file.path(config$session_dir, "checklist.txt")
  } else {
    session_checklist_path(config, index)
  }

  if (open_in_new_text_window) {
    if (!rstudioapi::isAvailable())
      stop("open_in_new_text_window = TRUE needs an active RStudio session.")
    term <- rstudioapi::terminalCreate(show = TRUE)
    rstudioapi::terminalSend(term, paste0("watch -n ", interval, " cat '", path, "'\n"))
    return(invisible(term))
  }

  n_prev_lines <- 0L
  repeat {
    content <- if (file.exists(path)) readLines(path) else
      "(waiting for checklist.txt to appear...)"
    if (n_prev_lines > 0) cat(sprintf("\033[%dA", n_prev_lines))
    cat("\033[J")
    cat(content, sep = "\n")
    cat("\n")
    n_prev_lines <- length(content) + 1L
    Sys.sleep(interval)
  }
}

#' Format a checklist data.table as a compact printable table
#' @param tab data.table, as returned by pipeline_checklist(print = FALSE)
#' @return character, one formatted line per stage
format_checklist <- function(tab) {
  lines <- vapply(seq_len(nrow(tab)), function(i) {
    row <- tab[i]
    mark <- switch(row$state, done = "\u2714", running = ">", "-")
    # "rate_label" %in% names(tab), not just row$rate_label -- a tab
    # built without this column (e.g. an existing caller/test
    # predating this field) would otherwise make row$rate_label
    # return NULL, and is.na(NULL) is logical(0), which errors inside
    # if(). Backward compatible: no column -> never rendered.
    rate_suffix <- if ("rate_label" %in% names(tab) && !is.na(row$rate_label)) {
      paste0(", ", row$rate_label)
    } else ""
    detail <- if (!is.na(row$active_experiment)) {
      sprintf(" -- active: %s (%d/%d samples%s)", row$active_experiment, row$active_done,
              row$active_total, rate_suffix)
    } else if (row$state == "running") {
      sprintf(" -- active (no per-sample detail available for this stage%s)", rate_suffix)
    } else ""
    sprintf("[%s] %-24s %d/%d studies done%s", mark, row$stage, row$done, row$total, detail)
  }, character(1))
  paste(lines, collapse = "\n")
}
