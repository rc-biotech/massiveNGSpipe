#' Reap zombie child processes of this R process
#'
#' processx (used by \code{callr::r_bg()}) replaces the SIGCHLD handler that
#' the parallel package installs once to reap forked children. After the first
#' \code{r_bg()} in a process, every later \code{parallel::mcparallel()} or
#' \code{BiocParallel::MulticoreParam} child stays a zombie: results still
#' arrive through pipes, but nobody calls \code{waitpid()}. Long-lived stage
#' workers mix both for hours, and on 2026-11-17 several thousand zombies hit the
#' container's process limit ("unable to fork"). Reproduced with plain
#' parallel + callr (R 4.5.3, processx 3.8.7); starting processx before the
#' first fork does not help, because every new \code{r_bg()} reinstalls its
#' handler.
#'
#' Reaps only children that are already zombies, and never the pids in
#' \code{keep} (live processx handles, whose exit status processx still needs).
#' KEEP THIS in the polling loops that mix callr and forks (see callers).
#' @param keep integer, pids not to reap
#' @return integer, number of reaped children (invisibly)
#' @noRd
reap_zombie_children <- function(keep = integer()) {
  if (.Platform$OS.type == "windows" || !dir.exists("/proc")) return(invisible(0L))
  me <- Sys.getpid(); n <- 0L
  for (st in Sys.glob("/proc/[0-9]*/stat")) {
    s <- tryCatch(suppressWarnings(readLines(st, warn = FALSE)), error = function(e) character())
    if (!length(s)) next
    f <- strsplit(sub("^.*\\) ", "", s[1]), " ", fixed = TRUE)[[1]]   # fields after "(comm) ": state, ppid, ...
    if (length(f) < 2 || f[1] != "Z" || f[2] != as.character(me)) next
    pid <- as.integer(sub(" .*", "", s[1]))
    if (pid %in% keep) next
    if (!inherits(tryCatch(suppressWarnings(sys::exec_status(pid, wait = FALSE)), error = function(e) e), "error")) n <- n + 1L
  }
  invisible(n)
}
