#' Detect colour-space (ABI SOLiD) FASTQ from the reads themselves
#'
#' Colour-space reads are a primer base followed by colour calls (0/1/2/3, '.' for no call), e.g. "T0120.3".
#' They must not enter base-space trimming/STAR, and a naive colour-to-base decoding propagates single colour errors
#' (codex-rna 2026-10-14, PRJNA603782). Instrument metadata alone does not decide this: SOLiD runs can also be
#' distributed in base space, so check the sequences.
#' @param file FASTQ(.gz) path
#' @param n_reads number of records inspected
#' @param min_frac fraction of inspected reads in colour-space alphabet needed to call it colour-space
#' @return TRUE if colour-space, FALSE otherwise (also for an empty/missing file)
#' @noRd
fastq_is_colorspace <- function(file, n_reads = 1000, min_frac = 0.9) {
  if (!file.exists(file)) return(FALSE)
  con <- if (grepl("\\.gz$", file)) gzfile(file, "r") else file(file, "r")
  on.exit(close(con))
  l <- readLines(con, n = 4 * n_reads, warn = FALSE)
  if (length(l) < 2) return(FALSE)
  seqs <- l[seq(2, length(l), by = 4)]
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
  if (!length(seqs)) return(FALSE)
  mean(grepl("^[ACGTN]?[0-3.]+$", seqs)) >= min_frac
}

#' Stop before base-space processing if any input FASTQ is colour-space
#' @param files FASTQ paths of one run
#' @param run Run id (for the message)
#' @return invisible(TRUE) if all are base space; otherwise an error
#' @noRd
assert_base_space_fastq <- function(files, run = "") {
  cs <- vapply(files, fastq_is_colorspace, logical(1))
  if (any(cs)) stop("Run ", run, ": colour-space (ABI SOLiD) FASTQ detected in ", paste(basename(files[cs]), collapse = ", "),
                    "; base-space trim/STAR refused (needs a colour-aware route)", call. = FALSE)
  invisible(TRUE)
}
