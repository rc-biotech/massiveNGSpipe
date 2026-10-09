# Dynamic BAM-sort RAM estimation + defer/crash logic for STAR alignment --
# built from a real production failure (PRJNA599943-homo_sapiens_RNA-seq,
# SRR11294057, 2026-10-09): STAR's --limitBAMsortRAM was a hardcoded 30GB
# (see ORFik commit 45fda30, which turned it into a real parameter), and a
# 19.84GB input needed ~62.88GB to sort -- never actually too large for
# this 502GB-RAM server, just for the fixed budget. See
# R/star_index_lock.R's own header comment for the matching style this
# follows (another resource-exhaustion incident this same session).

#' Rough proactive estimate of STAR's BAM-sort RAM need for one sample
#'
#' Calibrated from one real observed case (SRR11294057, 2026-10-09):
#' 19.84GB raw fastq needed ~62.88GB sort RAM, a ~3.17x ratio. Uses 4x as
#' a safety margin above that single data point, with a floor matching
#' the long-standing former hardcoded default so small files never get
#' an oddly tiny limit. A deliberately rough proactive guess, not the
#' authoritative number -- see \code{\link{estimate_bam_sort_ram_from_error}},
#' which reads STAR's own exact figure after a real failure instead.
#' Expect to recalibrate this multiplier as more real cases accumulate.
#' @param file1 character, the sample's (trimmed) input fastq path
#' @return numeric, bytes
#' @noRd
estimate_bam_sort_ram <- function(file1) {
  max(file.info(file1)$size * 4, 30e9)
}

#' STAR's own exact required value, parsed from its real failure text
#'
#' Far more reliable than any proactive guess -- STAR reports the EXACT
#' number it needed directly: "SOLUTION: re-run STAR with at least
#' --limitBAMsortRAM <N>". Used after a first attempt fails with this
#' specific error, to make the defer-vs-crash decision with real data
#' instead of a rough estimate.
#' @param err_message character, the caught error's conditionMessage()
#' @return numeric bytes, or NA_real_ if this isn't that specific error
#' @noRd
estimate_bam_sort_ram_from_error <- function(err_message) {
  m <- regmatches(err_message, regexpr("limitBAMsortRAM [0-9]+", err_message))
  if (length(m) == 0 || !nzchar(m)) return(NA_real_)
  as.numeric(sub("limitBAMsortRAM ", "", m))
}

#' Decide what to do about one sample's BAM-sort RAM need right now
#'
#' Three-way decision, directly from a real production incident's own
#' shape: if the file could never fit even with the whole machine to
#' itself, there's no point deferring or retrying -- crash now with a
#' clear message. If only CURRENTLY free memory is short but the total
#' would handle it, defer -- something else using memory right now may
#' free up by the time this sample is revisited. Otherwise proceed.
#' @param needed_bytes numeric
#' @param total_safety_fraction numeric, default 0.9 -- total memory
#' capped at this fraction before comparing, so a decision never relies
#' on using nearly 100% of system RAM for one sample's sort buffer (that
#' risks an OS-level OOM-kill affecting unrelated processes, not just
#' this one failing cleanly).
#' @return list(decision = one of "proceed"/"defer"/"crash", needed_gb,
#' total_gb, free_gb)
#' @noRd
bam_sort_ram_decision <- function(needed_bytes, total_safety_fraction = 0.9) {
  usage <- get_system_usage()
  total_gb <- usage$Memory_Total_GB * total_safety_fraction
  free_gb <- usage$Memory_Total_GB - usage$Memory_Usage_GB
  needed_gb <- needed_bytes / 1e9
  decision <- if (needed_gb > total_gb) "crash"
             else if (needed_gb > free_gb) "defer"
             else "proceed"
  list(decision = decision, needed_gb = needed_gb, total_gb = total_gb, free_gb = free_gb)
}

#' Path to one sample's "waiting for RAM" marker
#'
#' Same .rds-marker-file-existence-is-the-state convention as
#' R/pipeline_sample_flags.R, living alongside that step's own sample
#' flags.
#' @noRd
waiting_for_ram_path <- function(config, experiment, run)
  file.path(sample_flag_dir(config, "aligned", experiment), paste0(run, "_waiting_for_ram.rds"))

#' Mark one sample as waiting for RAM, if not already marked
#'
#' Deliberately does NOT reset an existing mark's timestamp -- the whole
#' point is tracking how long this sample has genuinely been waiting,
#' not how recently it was last checked.
#' @noRd
mark_waiting_for_ram <- function(config, experiment, run) {
  path <- waiting_for_ram_path(config, experiment, run)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(path)) saveRDS(Sys.time(), path)
  invisible(NULL)
}

#' Clear one sample's "waiting for RAM" mark, if any
#' @noRd
clear_waiting_for_ram <- function(config, experiment, run) {
  unlink(waiting_for_ram_path(config, experiment, run))
  invisible(NULL)
}

#' How long (minutes) one sample has been marked waiting for RAM
#' @return numeric minutes, or NA_real_ if it isn't currently marked
#' @noRd
waiting_for_ram_minutes <- function(config, experiment, run) {
  path <- waiting_for_ram_path(config, experiment, run)
  if (!file.exists(path)) return(NA_real_)
  as.numeric(difftime(Sys.time(), readRDS(path), units = "mins"))
}
