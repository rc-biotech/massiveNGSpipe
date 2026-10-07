#' Run one experiment's heavy step in a subprocess, polling for progress
#'
#' Wraps \code{callr::r_bg()} so that external-tool console output (fastp,
#' STAR, multiQC -- all invoked via \code{system()}/\code{system2()}, whose
#' output is NOT captured by BiocParallel's own \code{log=}/\code{logdir=}
#' mechanism, verified directly) lands in a log file instead of the calling
#' session's console, while the calling session itself stays free to poll
#' for progress (e.g. render the checklist) on a plain \code{Sys.sleep()}
#' loop. No extra thread or manually-managed background process is needed:
#' \code{callr::r_bg()} already gives an async handle, so the existing
#' loop that calls this function *is* the progress monitor.
#'
#' \code{func} must be fully self-contained: \code{package = "massiveNGSpipe"}
#' makes unexported massiveNGSpipe functions and \code{Depends}-only package
#' exports (data.table, ORFik) resolve unqualified inside it, but free
#' variables from the caller's enclosing scope are NOT available (verified)
#' -- every value \code{func} needs must be passed via \code{args}.
#'
#' @param func function, self-contained (see above)
#' @param args list, passed to \code{func}
#' @param logfile_out,logfile_err character, distinct file paths for this
#' call's stdout/stderr. These must NOT be the same path: verified that
#' \code{callr::r_bg()} corrupts output (silent truncation/interleaving)
#' when stdout and stderr are redirected to one shared file, because the
#' child's own R message stream and a nested \code{system()} grandchild's
#' writes race on the same file without coordination. Using two separate
#' files avoids this entirely.
#' @param on_poll function(), optional, called on every poll tick while the
#' subprocess is alive (e.g. to render a progress checklist). Takes no
#' arguments.
#' @param poll_interval numeric, seconds between polls, default 10.
#' @param install_race_retries integer, default 1. \code{package = "massiveNGSpipe"}
#' makes every single call here independently \code{loadNamespace()} this
#' same package from disk in a brand new process -- callr spawns a real
#' OS process, not a fork, so it shares no in-memory state with the
#' caller and cannot skip this. If a concurrent \code{R CMD INSTALL} of
#' massiveNGSpipe itself is in progress, that load can land on a
#' half-swapped install (R's own staged install does
#' \code{unlink(old); file.rename(new, old)} as two separate steps, not
#' one atomic op) and fail with a distinctive, transient error --
#' confirmed live, 2026-10-07, two unrelated studies
#' (PRJNA244941/PRJNA880902) both hit this within the same ~2 minute
#' window of an unrelated massiveNGSpipe reinstall mid-run. Retried
#' \code{install_race_retries} times (after \code{install_race_wait}
#' seconds each) ONLY when the error matches
#' \code{\link{is_install_race_error}}; any other error is raised
#' immediately on the first attempt, same as before this existed.
#' @param install_race_wait numeric, seconds to wait before a retry.
#' @return whatever \code{func} returns. Re-raises the child's error as a
#' normal R condition in the caller if \code{func} errored (matches
#' \code{callr::r()} semantics already relied on elsewhere), so existing
#' \code{try()}-based error handling one level up is unaffected. Before
#' returning (success, error, or interruption alike), both log files are
#' passed through \code{sanitize_progress_log()} to collapse any bare
#' carriage-return progress-bar spam a wrapped tool wrote into them.
run_experiment_subprocess <- function(func, args = list(),
                                      logfile_out, logfile_err,
                                      on_poll = NULL, poll_interval = 10,
                                      install_race_retries = 1,
                                      install_race_wait = 20) {
  stopifnot(logfile_out != logfile_err)
  dir.create(dirname(logfile_out), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(logfile_err), showWarnings = FALSE, recursive = TRUE)

  px <- NULL
  on.exit({
    if (!is.null(px) && px$is_alive()) px$kill_tree()
    sanitize_progress_log(logfile_out)
    sanitize_progress_log(logfile_err)
  }, add = TRUE)

  attempt <- 0
  repeat {
    attempt <- attempt + 1
    px <- callr::r_bg(func, args = args, package = "massiveNGSpipe",
                      stdout = logfile_out, stderr = logfile_err)
    while (px$is_alive()) {
      Sys.sleep(poll_interval)
      if (!is.null(on_poll)) on_poll()
    }
    outcome <- tryCatch(list(value = px$get_result()), error = function(e) list(error = e))
    if (is.null(outcome$error)) return(outcome$value)
    if (attempt <= install_race_retries && is_install_race_error(outcome$error)) {
      message("massiveNGSpipe lazy-load race detected in subprocess (attempt ", attempt,
              ") -- likely a concurrent massiveNGSpipe reinstall mid-run. Waiting ",
              install_race_wait, "s and retrying.")
      Sys.sleep(install_race_wait)
      next
    }
    stop(outcome$error)
  }
}

#' Does this condition look like the massiveNGSpipe install/lazy-load race?
#'
#' Matches the two concrete error messages confirmed live (see
#' \code{\link{run_experiment_subprocess}}'s own doc for the incident) plus
#' one anticipated variant (package briefly absent entirely, between R's
#' \code{unlink(old)} and \code{file.rename(new, old)}) that follows from
#' the same root cause but was not itself directly observed. Deliberately
#' scoped to massiveNGSpipe's own install path specifically, not a generic
#' "any lazy-load error" matcher, so an unrelated real corruption/package
#' problem elsewhere still surfaces immediately instead of being retried
#' and masked.
#' @param cnd a condition object (as caught by \code{tryCatch(..., error = )})
#' @return logical
#' @noRd
is_install_race_error <- function(cnd) {
  msg <- conditionMessage(cnd)
  grepl("massiveNGSpipe\\.rdb['\"]? is corrupt|read failed on .*massiveNGSpipe\\.rdb|there is no package called .massiveNGSpipe.",
        msg, perl = TRUE)
}

#' Collapse stray carriage-return progress spam in a captured log file
#'
#' Some external tools (a bare terminal progress meter, e.g. \code{aws s3
#' sync}'s default, or sratoolkit's \code{-p} flag) print live-updating
#' status with a bare carriage return, meant to overwrite one line in an
#' interactive terminal. Captured into a plain file by
#' \code{run_experiment_subprocess()}, that instead accumulates into one
#' enormous, unreadable line (verified directly against a real run: a
#' single real fetch produced one ~500KB line this way). Converting each
#' bare \code{\\r} into a real newline keeps every individual update as its
#' own line: turns an unreadable blob into a merely long, but normal, log
#' file.
#'
#' Best-effort only, and not a substitute for quieting a known-noisy call
#' at the source (see the \code{--no-progress}/\code{-q}/\code{progress_bar
#' = FALSE} settings in \code{download_sra_aws()}/\code{download_sra_ascp()}/
#' \code{sra_to_fastq()}) -- this exists as a backstop for tools we haven't
#' (or can't) silence that way.
#' @param path character, file to sanitize in place. No-op if it doesn't
#' exist, is empty, or contains no bare \code{\\r}.
#' @return invisible(NULL)
sanitize_progress_log <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  size <- file.info(path)$size
  if (is.na(size) || size == 0) return(invisible(NULL))
  txt <- readChar(path, size, useBytes = TRUE)
  if (!grepl("\r(?!\n)", txt, perl = TRUE, useBytes = TRUE)) return(invisible(NULL))
  txt <- gsub("\r\n", "\n", txt, fixed = TRUE, useBytes = TRUE)
  txt <- gsub("\r", "\n", txt, fixed = TRUE, useBytes = TRUE)
  cat(txt, file = path)
  invisible(NULL)
}
