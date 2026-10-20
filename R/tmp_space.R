# Disk-space claims for large temporary files on a shared drive.
#
# The server's root filesystem (holding /tmp and every R tempdir()) is shared
# by the OS and the home folders of ~15 docker containers. On 2026-10-19 it hit
# 100 % because two pipelines ran fasterq-dump into their R tempdirs at the same
# time: each checked free space on its own (check_tempdir_has_space_sra_to_fastq())
# and neither knew about the other. This file adds a pre-emptive claim: before
# writing N GB to a directory a process registers "I will write N GB here", and
# a new claim is only granted while
#   free space - outstanding claims of live processes - N >= min_free_gb
# Outstanding = claimed GB minus what the claim's own paths already hold on disk
# (that part is already subtracted from df's free space).
#
# Claims live in <dir>/.mngsp_space_claims/<id>.rds (pid, host, gb, paths). A
# claim whose process is gone (same host) is stale: its claim file AND the paths
# it registered (only paths inside <dir>) are removed on the next claim/release
# call by anyone, so a crashed session cannot leave 70 GB behind forever.
# Style follows R/star_index_lock.R (marker files + pid liveness); the claim
# table itself is guarded by a short mkdir() mutex because grant decisions
# must see each other.
#
# Claims only coordinate processes in this container (pids are per container);
# the min_free_gb headroom is what protects the other containers on the disk.

#' Minimum free space (GB) to leave on a temp drive after a claim
#'
#' Option \code{mngsp.tmp_min_free_gb}, default 75.
#' @return numeric
#' @noRd
tmp_min_free_gb <- function() as.numeric(getOption("mngsp.tmp_min_free_gb", 75))

tmp_space_claim_dir <- function(dir) file.path(dir, ".mngsp_space_claims")

#' Free space in GB on the drive holding dir (df, not cached)
#' @noRd
tmp_free_gb <- function(dir) {
  out <- suppressWarnings(system2("df", c("--output=avail", "-B1", shQuote(dir)), stdout = TRUE, stderr = FALSE))
  as.numeric(utils::tail(out, 1)) / 1024^3
}

pid_alive <- function(pid) {
  !is.na(pid) && identical(suppressWarnings(system2("kill", c("-0", pid), stdout = FALSE, stderr = FALSE)), 0L)
}

host_name <- function() unname(Sys.info()[["nodename"]])

#' TRUE for paths strictly inside dir (never dir itself), resolving symlinks/..
#' @noRd
path_inside <- function(paths, dir) {
  if (!length(paths)) return(logical())
  d <- paste0(normalizePath(dir, mustWork = FALSE), "/")
  p <- normalizePath(paths, mustWork = FALSE)
  startsWith(p, d) & nchar(p) > nchar(d)
}

#' Total size in GB of files/directories (recursive), missing paths count 0
#' @noRd
paths_size_gb <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(0)
  files <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE) else p
  }))
  sum(file.size(files), na.rm = TRUE) / 1024^3
}

#' Run expr while holding the claim-table mutex of dir
#' @noRd
with_space_claim_lock <- function(dir, expr, stale_secs = 120, timeout_secs = 300) {
  cdir <- tmp_space_claim_dir(dir)
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  lock <- file.path(cdir, ".lock")
  t0 <- Sys.time()
  while (!dir.create(lock, showWarnings = FALSE)) {
    age <- as.numeric(difftime(Sys.time(), file.mtime(lock), units = "secs"))
    if (!is.na(age) && age > stale_secs) unlink(lock, recursive = TRUE)
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout_secs)
      stop("Could not acquire the space-claim lock: ", lock)
    Sys.sleep(stats::runif(1, 0.05, 0.3))
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  force(expr)
}

#' Live space claims on dir; purges stale ones (and their temp paths)
#'
#' Must be called while holding the lock (\code{with_space_claim_lock}).
#' @return list of claims (each a list with id, pid, host, gb, paths, time)
#' @noRd
tmp_space_claims_locked <- function(dir) {
  files <- list.files(tmp_space_claim_dir(dir), pattern = "\\.rds$", full.names = TRUE)
  live <- list()
  for (f in files) {
    cl <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(cl)) { file.remove(f); next }
    if (identical(cl$host, host_name()) && !pid_alive(cl$pid)) {
      stale <- cl$paths[path_inside(cl$paths, dir)]
      if (length(stale)) {
        message("- Removing temp files of dead process ", cl$pid, " (", cl$id, "): ",
                paste(basename(stale), collapse = ", "))
        unlink(stale, recursive = TRUE)
      }
      file.remove(f)
      next
    }
    live[[length(live) + 1]] <- cl
  }
  live
}

#' Outstanding (not yet written) GB of a set of claims
#' @noRd
claims_outstanding_gb <- function(claims) {
  if (!length(claims)) return(0)
  sum(vapply(claims, function(cl) max(0, cl$gb - paths_size_gb(cl$paths)), numeric(1)))
}

#' Claim space for temporary files before writing them
#'
#' Registers "this process will write \code{gb} GB under \code{paths} in
#' \code{dir}". Granted only if, after subtracting every live claim's still
#' outstanding GB and this request, at least \code{min_free_gb} remain free.
#' Always pair a granted claim with \code{\link{tmp_space_release}} in
#' \code{on.exit()}; if the process dies instead, the next claim/release on
#' the same dir removes the claim and deletes its paths.
#' @param dir character, directory the files go to (e.g. \code{tempdir()} or
#'  its parent); claims are shared per dir, so use the same dir for all
#'  processes writing to one drive (\code{tmp_claim_root()})
#' @param gb numeric, GB this process will write at most
#' @param paths character, files/directories (inside dir) that the claim
#'  covers: counted against the claim while they grow and deleted on release
#'  or when the owner dies
#' @param min_free_gb numeric, headroom to keep, default
#'  \code{getOption("mngsp.tmp_min_free_gb", 75)}
#' @param free_fun function(dir) returning free GB (for tests)
#' @return claim id (character), or NULL when the space cannot be granted
#' @export
tmp_space_claim <- function(dir, gb, paths = character(), min_free_gb = tmp_min_free_gb(),
                            free_fun = tmp_free_gb) {
  stopifnot(is.numeric(gb), length(gb) == 1, !is.na(gb), gb >= 0)
  stopifnot(all(path_inside(paths, dir)))
  with_space_claim_lock(dir, {
    live <- tmp_space_claims_locked(dir)
    free <- free_fun(dir)
    pending <- claims_outstanding_gb(live)
    left <- free - pending - gb
    if (is.na(left) || left < min_free_gb) {
      message(sprintf("- Temp space claim refused on %s: need %.1f GB, free %.1f GB, claimed by others %.1f GB, keep %.0f GB",
                      dir, gb, free, pending, min_free_gb))
      return(NULL)
    }
    id <- paste0(Sys.getpid(), "_", format(Sys.time(), "%Y%m%d%H%M%OS3"), "_", sample.int(1e6, 1))
    saveRDS(list(id = id, pid = Sys.getpid(), host = host_name(), gb = gb, paths = paths, time = Sys.time()),
            file.path(tmp_space_claim_dir(dir), paste0(id, ".rds")))
    id
  })
}

#' Release a space claim and delete its temporary paths
#' @param dir character, as given to \code{\link{tmp_space_claim}}
#' @param id character, the claim id (NULL is a no-op)
#' @param cleanup logical, default TRUE: delete the claim's paths (only
#'  those inside dir)
#' @return invisible(NULL)
#' @export
tmp_space_release <- function(dir, id, cleanup = TRUE) {
  if (is.null(id)) return(invisible(NULL))
  with_space_claim_lock(dir, {
    f <- file.path(tmp_space_claim_dir(dir), paste0(id, ".rds"))
    cl <- if (file.exists(f)) tryCatch(readRDS(f), error = function(e) NULL) else NULL
    if (cleanup && !is.null(cl)) unlink(cl$paths[path_inside(cl$paths, dir)], recursive = TRUE)
    if (file.exists(f)) file.remove(f)
    tmp_space_claims_locked(dir) # opportunistic stale purge
  })
  invisible(NULL)
}

#' Shared claim root for the drive of tempdir(): the parent of tempdir()
#'
#' Every R session has its own tempdir() (e.g. /tmp/RtmpXXXX), but all of
#' them share /tmp's drive, so claims are kept in the parent directory.
#' @noRd
tmp_claim_root <- function() dirname(tempdir())

#' Report (and optionally delete) stale R temp directories left by dead sessions
#'
#' R removes its tempdir() on a normal exit; a crashed or killed session leaves
#' it behind, including partial fasterq-dump output (tens of GB) and FastQC
#' adapter-check reports. A directory is a deletion candidate only if all of:
#' it is an R session tempdir (\code{Rtmp*}) other than this session's, nothing
#' inside changed for \code{older_than_days}, no process has a file open in it
#' or uses it as working directory, and every top-level entry is a known
#' pipeline artefact (SRA/ENA run files, fasterq-dump temp dirs, the
#' \code{fasterq_*} work dirs, FastQC report dirs, adapter candidate tables).
#' Anything else is reported but never deleted.
#' @param dir character, default the parent of tempdir() (i.e. /tmp)
#' @param older_than_days numeric, default 2
#' @param delete logical, default FALSE (report only)
#' @return data.table(path, size_gb, newest, open, only_pipeline_files, candidate, deleted)
#' @export
tmp_cleanup_stale <- function(dir = tmp_claim_root(), older_than_days = 2, delete = FALSE) {
  dirs <- list.files(dir, pattern = "^Rtmp", full.names = TRUE, include.dirs = TRUE)
  dirs <- dirs[dir.exists(dirs) & normalizePath(dirs) != normalizePath(tempdir())]
  if (!length(dirs)) return(data.table::data.table())
  open_paths <- tryCatch({
    fds <- list.files("/proc", pattern = "^[0-9]+$", full.names = TRUE)
    unlist(lapply(fds, function(p) c(suppressWarnings(Sys.readlink(file.path(p, "cwd"))),
                                     suppressWarnings(Sys.readlink(list.files(file.path(p, "fd"), full.names = TRUE))))))
  }, error = function(e) character())
  open_paths <- open_paths[!is.na(open_paths) & nzchar(open_paths)]
  pattern <- "^((SRR|ERR|DRR)[0-9]+.*|fasterq\\.tmp\\..*|fasterq_.*|file[0-9a-f]+|adapter_candidates\\.txt)$"
  rows <- lapply(dirs, function(d) {
    files <- list.files(d, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    newest <- if (length(files)) max(file.mtime(files), na.rm = TRUE) else file.mtime(d)
    top <- list.files(d, all.files = TRUE, no.. = TRUE)
    fastqc_dirs <- top[grepl("^file[0-9a-f]+$", top)]
    fastqc_only <- all(vapply(file.path(d, fastqc_dirs), function(x) {
      all(grepl("^stdin_fastqc", list.files(x, all.files = TRUE, no.. = TRUE)))
    }, logical(1)))
    only_pipe <- all(grepl(pattern, top)) && fastqc_only
    is_open <- any(startsWith(open_paths, paste0(normalizePath(d), "/")) | open_paths == normalizePath(d))
    old <- as.numeric(difftime(Sys.time(), newest, units = "days")) > older_than_days
    data.table::data.table(path = d, size_gb = sum(file.size(files), na.rm = TRUE) / 1024^3, newest = newest,
                           open = is_open, only_pipeline_files = only_pipe,
                           candidate = old && !is_open && only_pipe, deleted = FALSE)
  })
  res <- data.table::rbindlist(rows)
  if (delete && any(res$candidate)) {
    for (i in which(res$candidate)) res$deleted[i] <- unlink(res$path[i], recursive = TRUE) == 0
  }
  res[order(-size_gb)]
}
