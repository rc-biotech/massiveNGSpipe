# Mitochondria-aware P-shifting for mitoribosome profiling.
#
# ORFik detects P-site offsets from the start-codon meta-profile of nuclear CDSs. Mitoribosomes translate only the
# mtDNA ORFs, most of which are leaderless (start codon at the mRNA 5' end), so the TIS profile is truncated and the
# nuclear signal is cytosolic contamination: automatic detection fails ("no periodic or accepted read-lengths",
# PRJNA542722, 2026-10-19) or, with a flat supplied table, is judged on the wrong genes. Analysis 2026-10-20
# (ribo_metadata_cleanup/reports/mitoribosome_shifting_2026-10-20.md): on mt CDSs human mitoribosome footprints
# show frame preference at the main lengths with offsets -12 (PRJNA632042) to -13 (PRJNA551033, 32-36 nt).
#
# For samples marked as mitoribosome profiling (FRACTION, see mito_fraction_terms()), offsets are estimated per
# read length on mitochondrial CDSs: a length is kept when it has enough reads and a clear frame preference in CDS
# bodies, and its offset is the value in 11..13 nt that puts the P-site on codon position 1 for the preferred
# frame (the 3-nt residue comes from the data, the range from the measured -12/-13 offsets). Samples with too few
# mt CDS reads fall back to -12 for 28-34 nt. The result is passed to ORFik as a known shift table, and every
# decision is written to <pshifted>/mito_shifts.csv.

#' FRACTION values that mark mitoribosome profiling (case-insensitive)
#' @noRd
mito_fraction_terms <- function() c("mito", "mitochondrial", "mitochondria", "mitoribosome", "mitoribosomes", "mitoriboseq")

#' Which rows of an experiment are mitoribosome profiling
#' @param df an ORFik experiment
#' @return logical, one per row
#' @noRd
is_mito_sample <- function(df) {
  fr <- if ("fraction" %in% colnames(df)) as.character(df$fraction) else rep("", nrow(df))
  fr[is.na(fr)] <- ""
  # ORFik joins varying columns into fraction (e.g. "Mitochondria__chl"): test the first token
  tolower(sub("__.*$", "", fr)) %in% mito_fraction_terms()
}

#' Mitochondrial CDSs of an experiment's annotation (>= 150 nt)
#' @noRd
mito_cds <- function(df) {
  cds <- ORFik::loadRegion(df, "cds")
  mt <- intersect(GenomeInfoDb::seqlevels(cds), c("MT", "chrM", "M", "Mt", "mitochondrion", "Mito", "chrMT"))
  if (!length(mt)) return(cds[0])
  chr <- as.character(GenomicRanges::seqnames(unlist(range(cds))))
  cds <- cds[chr %in% mt]
  cds <- cds[ORFik::widthPerGroup(cds, FALSE) >= 150]
  cds[!duplicated(names(cds))]
}

#' Per-length frame summary of 5' ends on mitochondrial CDS bodies
#' @param reads GRanges of reads (ofst import, score = read multiplicity)
#' @param cds GRangesList of mitochondrial CDSs
#' @param lengths integer, read lengths to evaluate
#' @return data.table(len, reads, best_frame, frame_frac)
#' @noRd
mito_frame_table <- function(reads, cds, lengths = 25:40) {
  empty <- data.table::data.table(len = integer(), reads = numeric(), best_frame = integer(), frame_frac = numeric())
  if (!length(cds) || !length(reads)) return(empty)
  reads <- reads[as.character(GenomicRanges::seqnames(reads)) %in% GenomeInfoDb::seqlevels(cds)]
  w <- S4Vectors::mcols(reads)$score  # ofst: collapsed reads, multiplicity in score
  if (is.null(w)) w <- rep(1, length(reads))
  len <- ORFik::readWidths(reads)
  keep <- len %in% lengths
  reads <- reads[keep]; w <- w[keep]; len <- len[keep]
  if (!length(reads)) return(empty)
  five <- GenomicRanges::resize(GenomicRanges::granges(reads), 1L, fix = "start")
  pos <- GenomicFeatures::mapToTranscripts(five, cds)
  if (!length(pos)) return(empty)
  xi <- S4Vectors::mcols(pos)$xHits
  cw <- ORFik::widthPerGroup(cds, TRUE)[as.character(GenomicRanges::seqnames(pos))]
  dt <- data.table::data.table(len = len[xi], w = w[xi], p = GenomicRanges::start(pos), cdsw = cw)
  dt <- dt[p > 30 & p < cdsw - 60]
  if (!nrow(dt)) return(empty)
  fr <- dt[, .(w = sum(w)), by = .(len, frame = (p - 1L) %% 3L)]
  fr[, `:=`(reads = sum(w), frac = w / sum(w)), by = len]
  best <- fr[order(-frac), .SD[1], by = len]
  best[, .(len, reads, best_frame = as.integer(frame), frame_frac = frac)][order(len)]
}

#' Offsets (negative, ORFik convention) from a mito frame table
#'
#' P-site = 5' end + offset; the offset in 11..13 whose residue puts the P-site of the preferred frame on codon
#' position 1 (frame f: 5' ends at CDS positions p with (p - 1) %% 3 == f, so offset o needs (f + o) %% 3 == 0).
#' @noRd
mito_offsets_from_frames <- function(ft, min_reads = 1000, min_frame_frac = 0.42) {
  ok <- ft[reads >= min_reads & frame_frac >= min_frame_frac]
  if (!nrow(ok)) return(data.table::data.table(fraction = integer(), offsets_start = integer()))
  off <- vapply(ok$best_frame, function(f) (11:13)[((f + 11:13) %% 3L) == 0L], integer(1))
  data.table::data.table(fraction = as.integer(ok$len), offsets_start = -as.integer(off))
}

#' Shift tables for the mitoribosome samples of an experiment
#' @param df an ORFik experiment (only rows with is_mito_sample() are used)
#' @param fallback_lengths,fallback_offset used when a sample has no length passing the frame test
#' @return named list (names = ofst paths) of data.tables(fraction, offsets_start); writes mito_shifts.csv
#' @noRd
mito_shift_tables <- function(df, fallback_lengths = 28:34, fallback_offset = -12L, min_reads = 1000,
                              min_frame_frac = 0.42) {
  rows <- which(is_mito_sample(df))
  if (!length(rows)) return(list())
  cds <- mito_cds(df)
  paths <- ofst_filepath(df[rows, ])
  log <- list(); tables <- list()
  for (k in seq_along(rows)) {
    reads <- ORFik::fimport(paths[k])
    ft <- mito_frame_table(reads, cds)
    tab <- mito_offsets_from_frames(ft, min_reads, min_frame_frac)
    method <- "mt_cds_frame"
    if (!nrow(tab)) {
      tab <- data.table::data.table(fraction = as.integer(fallback_lengths),
                                    offsets_start = rep(as.integer(fallback_offset), length(fallback_lengths)))
      method <- "fallback_default"
    }
    tables[[paths[k]]] <- tab
    log[[k]] <- merge(data.table::data.table(Run = ORFik::runIDs(df[rows[k], ]), len = tab$fraction,
                                             offset = tab$offsets_start, method = method),
                      ft, by = "len", all.x = TRUE)
    message("- Mito P-shift ", ORFik::runIDs(df[rows[k], ]), ": ", method, ", lengths ",
            paste0(tab$fraction, ":", tab$offsets_start, collapse = " "))
  }
  out_dir <- dirname(pshifted_filepath(df[rows[1], ]))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # keep the decisions of runs that were not reshifted this time (subset reruns), replace those that were
  new_log <- data.table::rbindlist(log, fill = TRUE)
  f <- file.path(out_dir, "mito_shifts.csv")
  if (file.exists(f)) {
    old <- data.table::fread(f, colClasses = c(Run = "character"))
    new_log <- data.table::rbindlist(list(old[!(Run %in% new_log$Run)], new_log), fill = TRUE)
  }
  data.table::fwrite(new_log[order(Run, len)], f)
  tables
}
