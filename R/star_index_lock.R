# Cross-process claims on a STAR shared-memory genome index, to stop one
# process from removing (STAR --genomeLoad Remove) an organism's index
# while a DIFFERENT process (a separate run_pipeline() -- online vs
# local, Ribo-seq vs RNA-seq -- or an ad hoc manual script) is still
# actively aligning against it. pipeline_align_one_organism()
# (R/pipeline_preset_steps_sub.R) is the only caller; this file is pure
# file-system bookkeeping, no STAR/ORFik calls of its own.
#
# Confirmed live, 2026-10-08 (PRJNA1051101-homo_sapiens): a different
# process's "I'm done with human, remove the index" call racing against
# this process's still-running STAR alignment produced a ~70-minute
# silent stall (zero CPU, zero network, zero file progress, every
# worker idle on do_select) rather than a clean error -- removing the
# shared-memory segment out from under a live STAR process does not
# fail loudly.
#
# Style matches the existing flag-marker-file convention
# (R/pipeline_flags.R, R/pipeline_sample_flags.R): a file's existence
# IS the state, no atomic/flock primitive. The one addition needed here
# that those don't have is PID-liveness staleness recovery, since a
# claim (unlike a flag) must not survive its owning process's death.

#' Directory holding one organism's active STAR-index claims
#' @param index character, the organism's STAR_index dir (the same
#' value already threaded through pipeline$organisms[[organism]]$index
#' and used to build file.path(index, "genomeDir") everywhere else)
#' @return character path, not guaranteed to exist yet
#' @noRd
star_index_lock_dir <- function(index) file.path(index, ".align_claims")

#' Register this process as an active user of one organism's STAR index
#'
#' Call once, in the PARENT process, before dispatching the subprocess
#' that will actually run STAR with keep.index.in.memory = TRUE/"y" for
#' this organism. Always pair with a guaranteed release (on.exit(),
#' added right after this call) -- a claim must not outlive the logical
#' operation it protects, even on error; the PID-liveness check in
#' \code{\link{star_index_active_claims}} only protects against this
#' process dying entirely, not against it staying alive but forgetting
#' to release.
#' @param index character, the organism's STAR_index dir
#' @param owner_id character, default a fresh id combining this
#' process's PID and the current time -- unique enough to not collide
#' with a concurrent claim from a different process or an earlier claim
#' by this same process.
#' @return character, the owner_id (pass this to
#' \code{\link{star_index_release}} and \code{\link{star_index_safe_to_remove}})
#' @noRd
star_index_claim <- function(index, owner_id = paste0(Sys.getpid(), "_", as.integer(Sys.time()))) {
  d <- star_index_lock_dir(index)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  saveRDS(Sys.getpid(), file.path(d, paste0(owner_id, ".rds")))
  owner_id
}

#' Release this process's own claim on one organism's STAR index
#' @param index character, the organism's STAR_index dir
#' @param owner_id character, as returned by \code{\link{star_index_claim}}
#' @return invisible(NULL)
#' @noRd
star_index_release <- function(index, owner_id) {
  f <- file.path(star_index_lock_dir(index), paste0(owner_id, ".rds"))
  if (file.exists(f)) file.remove(f)
  invisible(NULL)
}

#' Which OTHER claims on one organism's STAR index are still genuinely
#' active, purging any stale (dead-owner) claim files found along the way
#'
#' A claim file surviving its owning process's death (crash, kill -9 --
#' anything that skips the owner's own on.exit()-guaranteed release)
#' would otherwise block every future removal forever; this check's PID
#' liveness test is what prevents that, at the cost of only catching the
#' "owner process no longer exists at all" case, not "owner process
#' still alive but got stuck/forgot to release" (guaranteed release via
#' on.exit() at the call site is what covers that case instead).
#' @param index character, the organism's STAR_index dir
#' @param exclude_owner character, this process's OWN owner_id (never
#' counted against itself, including a not-yet-released one still live
#' on disk)
#' @return character vector of still-active owner_ids (possibly empty)
#' @noRd
star_index_active_claims <- function(index, exclude_owner) {
  d <- star_index_lock_dir(index)
  if (!dir.exists(d)) return(character())
  files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  owners <- sub("\\.rds$", "", basename(files))
  keep <- owners != exclude_owner
  files <- files[keep]; owners <- owners[keep]
  if (length(files) == 0) return(character())

  alive <- vapply(files, function(f) {
    pid <- tryCatch(readRDS(f), error = function(e) NA_integer_)
    if (is.na(pid)) return(FALSE)
    # kill -0 signals nothing, just checks the PID exists and is ours
    # to signal; exit status 0 means alive, non-zero means gone.
    identical(suppressWarnings(system2("kill", c("-0", pid),
                                       stdout = FALSE, stderr = FALSE)), 0L)
  }, logical(1))

  stale <- files[!alive]
  if (length(stale) > 0) file.remove(stale)
  owners[alive]
}

#' Is it safe for THIS process to remove one organism's STAR index right now?
#' @inheritParams star_index_active_claims
#' @return logical, TRUE only when no other process holds an active claim
#' @noRd
star_index_safe_to_remove <- function(index, exclude_owner) {
  length(star_index_active_claims(index, exclude_owner)) == 0
}

#' Force-remove an organism's STAR index from shared memory if (and only
#' if) nothing currently holds a claim on it
#'
#' A maintenance utility, not wired into \code{\link{run_pipeline}}
#' automatically: the claim/release logic in
#' \code{\link{pipeline_align_one_organism}} deliberately errs toward
#' "keep loaded" whenever it cannot prove removal is safe, so an index
#' can end up staying resident in shared memory longer than strictly
#' necessary (never INCORRECTLY -- this only ever trades a little extra
#' memory for safety). Run this by hand between sessions, or from a
#' maintenance script, to reclaim shared memory for every organism
#' nothing is actively using. Uses the same manual remedy already
#' documented in ORFik's own STAR_Aligner/RNA_Align_pipeline.sh usage
#' text ("if STAR is stuck on load, run this line:
#' STAR --genomeDir <genomeDir> --genomeLoad Remove").
#' @param config the mNGSp config object
#' @param star.path character, default \code{ORFik::STAR.install()}
#' @return data.table(organism, index, removed, reason) -- one row per
#' organism folder found under \code{config$config["ref"]} that has a
#' STAR_index/genomeDir subfolder
#' @export
star_index_gc <- function(config, star.path = STAR.install()) {
  ref_base <- config$config["ref"]
  organism_dirs <- list.dirs(ref_base, recursive = FALSE)
  rows <- lapply(organism_dirs, function(org_dir) {
    index <- file.path(org_dir, "STAR_index")
    genome_dir <- file.path(index, "genomeDir")
    if (!dir.exists(genome_dir)) return(NULL)
    safe <- star_index_safe_to_remove(index, exclude_owner = "")
    removed <- FALSE
    if (safe) {
      status <- system2(star.path, c("--genomeDir", shQuote(genome_dir), "--genomeLoad", "Remove"),
                        stdout = FALSE, stderr = FALSE)
      removed <- identical(status, 0L)
    }
    data.table::data.table(organism = basename(org_dir), index = index, removed = removed,
                           reason = if (!safe) "other process holds an active claim" else if (!removed) "STAR remove call failed or nothing was loaded" else "ok")
  })
  data.table::rbindlist(rows[!vapply(rows, is.null, logical(1))], fill = TRUE)
}
