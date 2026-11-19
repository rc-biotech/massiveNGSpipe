# Single-read post-shift QC (t79).
#
# Profile of the legacy block on GSE73553-mus_musculus (884 libraries, t79/PROFILE.md): almost all time was spent in
# ORFik::fimport() (GRanges build, ~10 s per 1.2 M-row library) and GRanges coverage code, repeated 4 times per library
# (explicit shiftPlots on ALL libraries, shiftPlots on 39, orfFrameDistributions, save_pshifted_length_distributions).
# The helpers below read each pshifted .ofst ONCE as a plain data.table (fst, no GRanges) and derive from that single
# read: the frame-per-read-length table (same numbers as ORFik::regionPerReadLength(withFrames, "frameSumPerL")),
# the read-length distribution (same as read_length_distribution()), and, only for libraries that are plotted, the GRanges
# needed by ORFik::windowPerReadLength() (the plot data are unchanged).

#' Read a pshifted/ofst file once as data.table (no GRanges)
#' @param path character, one .ofst file
#' @return data.table with columns seqnames (character), start, strand (character), size, score; 0 rows when empty
#' @noRd
qc_read_ofst <- function(path) {
  cols <- fst::metadata_fst(path)$colNames
  want <- intersect(c("seqnames", "start", "width", "strand", "size", "score"), cols)
  dt <- fst::read_fst(path, columns = want, as.data.table = TRUE)
  if (is.null(dt$score)) dt[, score := rep(1L, .N)]
  # ORFik::readWidths(): width-1 (p-shifted) reads use the size column (else score); wider reads their width
  if (nrow(dt) == 0) return(dt)
  one_based <- is.null(dt$width) || all(dt$width == 1L)
  dt[, size := if (one_based) { if (is.null(dt$size)) score else dt$size } else width]
  dt[, seqnames := as.character(seqnames)][, strand := as.character(strand)]
  dt[]
}

#' Read-length distribution table from an already-read library (same output as read_length_distribution())
#' @noRd
qc_length_distribution <- function(dt) {
  if (nrow(dt) == 0) {
    return(data.table::data.table(read_length = integer(), count = numeric(), percent = numeric()))
  }
  agg <- dt[, .(count = sum(score)), by = .(read_length = size)][order(read_length)]
  agg[, percent := round(100 * count / sum(count), 2)]
  agg[]
}

#' Flat CDS exon table with the 0-based position of each exon start inside its CDS, built ONCE per study
#' @param cds GRangesList of CDS (exons in transcription order, as loadRegion(df, "cds"))
#' @return data.table keyed on seqnames, start, end
#' @noRd
qc_cds_exon_table <- function(cds) {
  ul <- unlist(cds, use.names = FALSE)
  nm <- rep.int(seq_along(cds), lengths(cds))
  e <- data.table::data.table(seqnames = as.character(GenomicRanges::seqnames(ul)),
                              start = BiocGenerics::start(ul), end = BiocGenerics::end(ul),
                              strand = as.character(BiocGenerics::strand(ul)), tx = nm)
  e[, w := end - start + 1L]
  e[, before := cumsum(w) - w, by = tx]
  e[, ord := seq_len(.N)]
  data.table::setkey(e, seqnames, start, end)
  e[]
}

#' Frame-sum per read length for one library, from the already-read table
#'
#' Reproduces \code{ORFik::regionPerReadLength(cds, lib, withFrames = TRUE, scoring = "frameSumPerL", weight = "score",
#' drop.zero.dt = TRUE)} for width-1 (p-shifted) reads: a read is counted at its position in every CDS it hits on the
#' same strand, frame = (position in CDS - 1) mod 3, rows ordered as ORFik returns them (read length ascending, frames in
#' order of first coverage along CDS order).
#' @param dt data.table from \code{qc_read_ofst()}
#' @param exons data.table from \code{qc_cds_exon_table()}
#' @return data.table(fraction = read length, frame, score) (the caller renames like shift_qc_one_sample())
#' @noRd
qc_frame_sums <- function(dt, exons) {
  empty <- data.table::data.table(fraction = integer(), frame = numeric(), score = numeric())
  if (nrow(dt) == 0) return(empty)
  r <- dt[, .(score = sum(score)), by = .(seqnames, start, strand, size)]
  r[, end := start]
  data.table::setkey(r, seqnames, start, end)
  h <- data.table::foverlaps(r, exons, by.x = c("seqnames", "start", "end"), by.y = c("seqnames", "start", "end"),
                             type = "within", nomatch = NULL)
  # foverlaps: columns of x (r) are prefixed i. when names clash
  h <- h[strand == i.strand]
  if (nrow(h) == 0) return(empty)
  h[, pos0 := before + data.table::fifelse(strand == "-", end - i.start, i.start - start)]
  h[, frame := as.numeric(pos0 %% 3L)]
  h[, firstkey := tx * 4294967296 + pos0]
  res <- h[, .(score = sum(score), firstkey = min(firstkey)), by = .(fraction = size, frame)]
  data.table::setorder(res, fraction, firstkey)
  res[, firstkey := NULL]
  data.table::setcolorder(res, c("fraction", "frame", "score"))   # ORFik column order
  res[]
}

#' TIS-window exon table (positions -upstream..downstream around the CDS start), built ONCE per study
#'
#' The windows are what \code{ORFik::windowPerReadLength()} builds itself with \code{startRegion(cds, mrna, TRUE, upstream,
#' downstream)}; here they are flattened like \code{qc_cds_exon_table()} so every library can be counted with one join.
#' @return list(exons = data.table, W = window width, n = number of windows) or NULL when windows differ in width
#' (the caller then falls back to ORFik)
#' @noRd
qc_window_table <- function(plot_cds, plot_mrna, upstream, downstream) {
  windows <- ORFik::startRegion(plot_cds, plot_mrna, TRUE, upstream, downstream)
  W <- unique(ORFik::widthPerGroup(windows))
  if (length(W) != 1L) return(NULL)
  list(exons = qc_cds_exon_table(windows), W = as.integer(W), n = length(windows))
}

#' Shift-barplot data (hitMap) for one library from the already-read table
#'
#' Reproduces \code{ORFik::windowPerReadLength(plot_cds, plot_mrna, lib, upstream, downstream)} (scoring
#' "transcriptNormalized", weight "score", drop.zero.dt FALSE): the per-window, per-position read counts are made with one
#' join (as \code{qc_frame_sums()}) instead of ORFik's GRanges coverage; the SAME \code{ORFik::coverageScorings()} turns them
#' into the score, so the numbers are computed by the same summation.
#' @param dt data.table from \code{qc_read_ofst()}
#' @param wt result of \code{qc_window_table()}
#' @return data.table(position, score, fraction), as windowPerReadLength()
#' @noRd
qc_hitmap <- function(dt, wt, upstream) {
  if (nrow(dt) == 0) return(data.table::data.table(position = integer(), score = integer(), fraction = integer()))
  r <- dt[, .(seqnames, start, strand, size, score)]
  r[, end := start]
  data.table::setkey(r, seqnames, start, end)
  h <- data.table::foverlaps(r, wt$exons, by.x = c("seqnames", "start", "end"), by.y = c("seqnames", "start", "end"),
                             type = "within", nomatch = NULL)
  h <- h[strand == i.strand]
  h[, pos0 := before + data.table::fifelse(strand == "-", end - i.start, i.start - start)]
  agg <- h[, .(count = sum(score)), by = .(size, tx, pos0)]
  W <- wt$W; n <- wt$n
  genes <- rep.int(seq_len(n), rep.int(W, n)); position <- rep.int(seq_len(W), n) - (as.integer(upstream) + 1L)
  out <- lapply(sort(unique(dt$size)), function(l) {
    a <- agg[size == l]
    full <- data.table::data.table(count = numeric(n * W), genes = genes, position = position)
    data.table::set(full, i = (a$tx - 1L) * W + a$pos0 + 1L, j = "count", value = as.numeric(a$count))
    sc <- ORFik::coverageScorings(full, "transcriptNormalized", copy.dt = FALSE)
    sc[, fraction := l]
    sc
  })
  data.table::rbindlist(out)
}
