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
#' @return whatever \code{func} returns. Re-raises the child's error as a
#' normal R condition in the caller if \code{func} errored (matches
#' \code{callr::r()} semantics already relied on elsewhere), so existing
#' \code{try()}-based error handling one level up is unaffected. Before
#' returning (success, error, or interruption alike), both log files are
#' passed through \code{sanitize_progress_log()} to collapse any bare
#' carriage-return progress-bar spam a wrapped tool wrote into them.
run_experiment_subprocess <- function(func, args = list(),
                                      logfile_out, logfile_err,
                                      on_poll = NULL, poll_interval = 10) {
  stopifnot(logfile_out != logfile_err)
  dir.create(dirname(logfile_out), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(logfile_err), showWarnings = FALSE, recursive = TRUE)
  px <- callr::r_bg(func, args = args, package = "massiveNGSpipe",
                    stdout = logfile_out, stderr = logfile_err)
  on.exit({
    if (px$is_alive()) px$kill_tree()
    sanitize_progress_log(logfile_out)
    sanitize_progress_log(logfile_err)
  }, add = TRUE)
  while (px$is_alive()) {
    Sys.sleep(poll_interval)
    if (!is.null(on_poll)) on_poll()
  }
  px$get_result()
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
