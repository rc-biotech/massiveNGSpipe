# Per-sample-cached replacement for shift_qc()'s whole-experiment
# reload, confirmed live on GSE151959-homo_sapiens, 2026-10-07: fixing
# ONE sample out of 25 and rerunning cost ~33 real minutes once the
# ofst/covrle/bigwig per-sample resume was working correctly, and
# nearly all of that remaining time was this one stage. Root cause:
# shift_qc() reloads EVERY sample's pshifted library from disk TWICE
# per call -- once inside shiftPlots() for the shift-barplot hitMap,
# once again inside orfFrameDistributions()->outputLibs() for the
# frame/periodicity table -- with no per-sample cache, so an unchanged
# sibling's data is reloaded and recomputed exactly as much as the one
# sample that actually changed. shiftPlots() additionally accepts a
# BPPARAM argument that is never used anywhere in its body (confirmed
# by reading ORFik's source directly) -- its per-sample loop is
# hard-serial.
#
# shift_qc_cached() replaces that call: each sample is loaded AT MOST
# ONCE per call (the barplot hitMap and the frame table are both
# derived from the same in-memory library), results are cached
# per-sample on disk under QCfolder(df)/per_sample_qc/, and a sample
# whose cache is already newer than its own pshifted file is never
# reloaded at all. The combined outputs (Ribo_frames_all.csv,
# Ribo_frames_badzero.csv, pshifts_barplots.png, the good/warning/
# no_data verdict) are therefore just a concatenation/reassembly of
# per-sample pieces -- cached or freshly computed -- not a fresh
# whole-experiment pass every time.
#
# A sample whose single-load computation errors (shift_qc_one_sample()
# returning a "shift_qc_sample_error") gets no cache written and is
# excluded from this call's aggregate, rather than failing the whole
# experiment's QC or being silently counted as "zero reads" -- it will
# be retried (not skipped) on the next call, since no cache means
# shift_qc_cache_valid() stays FALSE for it.

#' Directory holding one experiment's cached per-sample pshift-QC pieces
#'
#' Namespaced by mapper mode (\code{uniqueMappers(df)}): \code{QCfolder(df)}
#' itself does NOT vary with mapper mode (confirmed in ORFik's own
#' \code{QCfolder()} method -- both passes share one \code{QC_STATS/}
#' dir), but the two passes' underlying pshifted files are at different
#' paths and must never validate each other's cache. Without this, the
#' unique-mappers pass could see a just-written all-mappers cache whose
#' mtime happens to postdate ITS OWN (different-path) pshifted file and
#' wrongly treat it as already fresh -- silently skipping real
#' computation and serving the wrong mapper mode's data.
#' @param df an ORFik experiment
#' @return character, directory path (not guaranteed to exist yet)
#' @noRd
shift_qc_cache_dir <- function(df) {
  file.path(QCfolder(df), "per_sample_qc", if (uniqueMappers(df)) "unique" else "all")
}

#' The cache file paths for one sample (one row of df)
#' @param df_one_row an ORFik experiment subset to exactly one row
#' @return list(hitmap = character, frames = character)
#' @noRd
shift_qc_cache_paths <- function(df_one_row) {
  run_id <- ORFik::runIDs(df_one_row)
  cache_dir <- shift_qc_cache_dir(df_one_row)
  list(hitmap = file.path(cache_dir, paste0(run_id, "_hitmap.rds")),
       frames = file.path(cache_dir, paste0(run_id, "_frames.csv")))
}

#' Thin wrapper around \code{ORFik::filepath(df_one_row, "pshifted")}
#'
#' Exists purely so tests can mock ONE plain massiveNGSpipe function
#' instead of needing \code{filepath()} itself to dispatch on a fake
#' experiment class -- \code{filepath()} is a plain function in ORFik
#' (not an S4 generic) with an internal \code{stopifnot(is(df, "experiment"))},
#' so it can never be given a method for a test stub class.
#' @param df_one_row an ORFik experiment subset to exactly one row
#' @return character, path to the pshifted file
#' @noRd
pshifted_filepath <- function(df_one_row) filepath(df_one_row, "pshifted")

#' Is this sample's QC cache present and at least as new as its pshifted file?
#'
#' The freshness check (not just existence) is what makes a redone
#' sample (e.g. a barcode fix followed by a reshift) correctly trigger
#' recomputation -- its pshifted file gets a newer mtime than the old
#' cache, which was written against the PREVIOUS version of that file.
#' @param df_one_row an ORFik experiment subset to exactly one row
#' @return logical
#' @noRd
shift_qc_cache_valid <- function(df_one_row) {
  pshifted_path <- pshifted_filepath(df_one_row)
  if (!file.exists(pshifted_path)) return(FALSE)
  paths <- unlist(shift_qc_cache_paths(df_one_row))
  if (!all(file.exists(paths))) return(FALSE)
  cache_mtime <- min(file.info(paths)$mtime)
  cache_mtime >= file.info(pshifted_path)$mtime
}

#' Compute one sample's shift-QC pieces from exactly ONE library load
#'
#' Loads this sample's pshifted library once and derives both the
#' shift-barplot \code{hitMap} (\code{ORFik::shiftPlots()}'s own
#' per-sample computation) and the frame/periodicity table
#' (\code{ORFik::regionPerReadLength()}'s own per-library computation,
#' as used inside \code{orfFrameDistributions()}) from that single
#' in-memory object -- see this file's header for why that matters.
#' @param df_one_row an ORFik experiment subset to exactly one row
#' @param cds,mrna GRangesList, loaded ONCE for the whole experiment by
#' the caller (never reloaded per sample)
#' @param upstream,downstream numeric, passed to \code{windowPerReadLength()}
#' @param weight character, default "score"
#' @return list(hitmap = data.table, frames = data.table, name_short =
#' character) on success. On failure, an object of class
#' \code{c("shift_qc_sample_error", "error", "condition")} with a
#' \code{message} element -- the caller must skip caching for this
#' sample and exclude it from aggregation (see \code{\link{shift_qc_cached}}),
#' never treat a failure as if the sample had zero reads.
#' @noRd
shift_qc_one_sample <- function(df_one_row, cds, mrna, upstream, downstream, weight = "score") {
  tryCatch({
    style <- seqinfo(df_one_row)
    lib <- ORFik::fimport(pshifted_filepath(df_one_row), style)
    name_short <- ORFik::bamVarName(df_one_row, skip.experiment = TRUE)

    hitmap <- ORFik::windowPerReadLength(cds, mrna, lib, upstream = upstream, downstream = downstream)
    hitmap[, frame := position %% 3]

    frames <- ORFik::regionPerReadLength(cds, lib, withFrames = TRUE, scoring = "frameSumPerL",
                                         weight = weight, drop.zero.dt = TRUE,
                                         exclude.zero.cov.grl = length(lib) == 0,
                                         BPPARAM = BiocParallel::SerialParam())
    if (nrow(frames) > 0) {
      frames[, length := fraction]
      frames[, fraction := rep(name_short, nrow(frames))]
    }
    list(hitmap = hitmap, frames = frames, name_short = name_short)
  }, error = function(e) {
    structure(list(message = conditionMessage(e)),
             class = c("shift_qc_sample_error", "error", "condition"))
  })
}

#' Write one sample's shift-QC cache
#' @param df_one_row an ORFik experiment subset to exactly one row
#' @param result a successful \code{\link{shift_qc_one_sample}} result
#' @return invisible(NULL)
#' @noRd
write_shift_qc_cache <- function(df_one_row, result) {
  paths <- shift_qc_cache_paths(df_one_row)
  dir.create(dirname(paths$hitmap), showWarnings = FALSE, recursive = TRUE)
  saveRDS(result$hitmap, paths$hitmap)
  data.table::fwrite(result$frames, paths$frames)
  invisible(NULL)
}

#' Read one sample's shift-QC cache back
#' @param df_one_row an ORFik experiment subset to exactly one row
#' @return list(hitmap = data.table, frames = data.table)
#' @noRd
read_shift_qc_cache <- function(df_one_row) {
  paths <- shift_qc_cache_paths(df_one_row)
  list(hitmap = readRDS(paths$hitmap), frames = data.table::fread(paths$frames))
}

#' Cached, per-sample replacement for \code{shift_qc()}
#'
#' Same external contract as \code{shift_qc()} (writes
#' \code{Ribo_frames_all.csv}, \code{Ribo_frames_badzero.csv},
#' \code{pshifts_barplots.png}, the good/warning/no_data verdict flag,
#' and the adapter/barcode QC check) but loads each sample's pshifted
#' library AT MOST once per call, and only for samples whose cache is
#' missing or stale -- see this file's own header.
#' @param df an ORFik experiment
#' @param BPPARAM BiocParallel param, used to parallelize across
#' samples that need (re)computation
#' @param max_no_adapter_removed_pct passed to
#' \code{check_adapter_barcode_quality()}
#' @return invisible(NULL)
#' @noRd
#' Transcript annotation context \code{shift_qc_one_sample()} needs:
#' upstream/downstream window sizes plus cds/mrna regions
#'
#' Factored out of \code{shift_qc_cached()} itself so it's one plain,
#' directly mockable massiveNGSpipe function in tests, the same way
#' \code{shift_qc_one_sample()}/\code{shift_qc_build_combined_plot()}
#' already are -- \code{ORFik::filterTranscripts()}/\code{loadRegion()}
#' need a real transcript annotation to run at all, which a lightweight
#' test fixture (no real genome) can't provide.
#' @param df an ORFik experiment
#' @return list(upstream, downstream, cds, mrna)
#' @noRd
shift_qc_annotation_context <- function(df) {
  has_leaders <- length(filterTranscripts(df, 5, 0, 0, stopOnEmpty = FALSE))
  list(upstream = ifelse(has_leaders, 5, 0), downstream = 20,
      cds = loadRegion(df, part = "cds"), mrna = loadRegion(df, part = "mrna"))
}

shift_qc_cached <- function(df, BPPARAM = bpparam(), max_no_adapter_removed_pct = 80) {
  subset <- if (nrow(df) >= 40) seq(39) else seq(nrow(df))

  needs_compute <- !vapply(seq_len(nrow(df)), function(i) shift_qc_cache_valid(df[i, ]), logical(1))

  if (any(needs_compute)) {
    # Annotation (transcript filtering, cds/mrna) is only ever needed
    # to compute a sample's QC from scratch -- skip it entirely when
    # every sample's cache is already valid, instead of always paying
    # for a genome/annotation load regardless of whether anything
    # changed.
    ctx <- shift_qc_annotation_context(df)
    upstream <- ctx$upstream; downstream <- ctx$downstream
    cds <- ctx$cds; mrna <- ctx$mrna
    idx <- which(needs_compute)
    results <- BiocParallel::bplapply(idx, function(i, df, cds, mrna, upstream, downstream) {
      res <- shift_qc_one_sample(df[i, ], cds, mrna, upstream, downstream)
      # Write AS SOON AS this one sample finishes, inside the worker,
      # instead of returning it to the master for a deferred write
      # after the WHOLE batch completes. Each sample's cache files are
      # its own, distinct paths (shift_qc_cache_paths()), so concurrent
      # workers writing their own sample's files is safe. This is what
      # makes an interrupted run actually resumable: a kill/crash
      # partway through the batch now keeps every already-finished
      # sample cached, instead of losing the whole batch's progress
      # because nothing had been persisted yet. Confirmed live,
      # 2026-10-08: the OLD deferred-write design meant a multi-hour
      # PRJNA637713 valid_pshift run, if ever killed, would have lost
      # every sample's work regardless of how many had already
      # finished computing well before the kill.
      if (!inherits(res, "shift_qc_sample_error")) write_shift_qc_cache(df[i, ], res)
      res
    }, df = df, cds = cds, mrna = mrna, upstream = upstream, downstream = downstream,
       BPPARAM = BPPARAM)
    for (j in seq_along(idx)) {
      res <- results[[j]]
      if (inherits(res, "shift_qc_sample_error")) {
        warning("shift_qc_cached(): sample ", ORFik::runIDs(df[idx[j], ]),
               " failed, skipping cache for this sample: ", res$message)
      }
    }
  }

  ok <- vapply(seq_len(nrow(df)), function(i) shift_qc_cache_valid(df[i, ]), logical(1))
  QCFolder <- QCfolder(df)
  if (!any(ok)) {
    warning("shift_qc_cached(): no sample succeeded, nothing to aggregate for ", name(df))
    return(invisible(NULL))
  }
  ok_idx <- which(ok)
  cached <- lapply(ok_idx, function(i) read_shift_qc_cache(df[i, ]))

  # fill = TRUE: a sample whose frames table came back with ZERO rows
  # (e.g. no periodicity-covered regions at all) never goes through
  # shift_qc_one_sample()'s own `if (nrow(frames) > 0)` column-
  # augmentation (length/fraction renaming), so its cached frames.csv
  # keeps a shorter column set (3 columns: fraction/frame/score) than
  # a non-empty sample's (4: +length). rbindlist() errors on that
  # mismatch without fill=TRUE -- safe to add here specifically
  # because a 0-row table contributes zero ROWS to the result
  # regardless of its own column set, so this can never introduce an
  # NA-padded row into frameQC; it only lets a genuinely-empty sample's
  # table merge in as a no-op instead of crashing the whole
  # experiment's aggregation. Confirmed live, 2026-10-08,
  # PRJNA637713-zea_mays/SRR13808095 (21-byte, header-only frames.csv).
  frameQC <- data.table::rbindlist(lapply(cached, `[[`, "frames"), fill = TRUE)
  if (nrow(frameQC) > 0) {
    frameQC[, frame := as.factor(frame)]
    frameQC[, fraction := as.factor(fraction)]
    frameQC[, percent := (score / sum(score)) * 100, by = fraction]
    frameQC[, percent_length := (score / sum(score)) * 100, by = .(fraction, length)]
    frameQC[, best_frame := (percent_length / max(percent_length)) == 1, by = .(fraction, length)]
    frameQC[, fraction := factor(fraction, levels = unique(fraction), ordered = TRUE)]
  } else {
    frameQC <- data.table::data.table(frame = factor(), fraction = factor(), length = numeric(),
                                      score = numeric(), percent = numeric(),
                                      percent_length = numeric(), best_frame = logical())
  }
  zero_frame <- frameQC[frame == 0 & best_frame == FALSE, ]

  data.table::fwrite(frameQC, file = file.path(QCFolder, "Ribo_frames_all.csv"))
  data.table::fwrite(zero_frame, file = file.path(QCFolder, "Ribo_frames_badzero.csv"))

  plot_idx <- intersect(ok_idx, subset)
  if (length(plot_idx) > 0) shift_qc_build_combined_plot(df, plot_idx)

  trimmed_dir <- file.path(bam_dir_from_df(df), "trim")
  check_adapter_barcode_quality(bam_dir_from_df(df), trimmed_dir, max_no_adapter_removed_pct)
  periodicity_check_flag(zero_frame, QCFolder, has_cds_data = nrow(frameQC) > 0)
  invisible(NULL)
}

#' Reassemble the combined shift-barplot from cached per-sample hitMaps
#'
#' Mirrors \code{ORFik::shiftPlots()}'s own plot-assembly tail exactly
#' (same \code{pSitePlot()} call, same \code{arrangeGrob()}/
#' \code{ggsave()} sizing formulas, same output paths/filenames), but
#' reads each sample's \code{hitMap} from its cache file instead of
#' calling \code{fimport()} -- no library reload needed for this part
#' at all.
#' @param df an ORFik experiment
#' @param plot_idx integer, row indices (into \code{df}) to include,
#' already capped to \code{shift_qc_cached()}'s own max-39 subset
#' @return invisible(NULL)
#' @noRd
shift_qc_build_combined_plot <- function(df, plot_idx) {
  lib_names <- ORFik::bamVarName(df[plot_idx, ], skip.experiment = TRUE)
  plots <- lapply(seq_along(plot_idx), function(x) {
    hitmap <- read_shift_qc_cache(df[plot_idx[x], ])$hitmap
    miniTitle <- gsub("_", " ", lib_names[x])
    p <- ORFik::pSitePlot(hitmap, scoring = "transcriptNormalized",
                          facet = TRUE, frameSum = TRUE, title = miniTitle)
    attr(p, "data") <- hitmap
    p
  })
  data <- data.table::rbindlist(lapply(plots, function(p) attr(p, "data")), idcol = TRUE)
  res <- do.call(gridExtra::arrangeGrob, c(plots, ncol = 1, top = "Ribo-seq"))

  output <- file.path(QCfolder(df), "pshifts_barplots.png")
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  dpi <- ifelse(length(plot_idx) < 22, 300, 200)
  height_scaler <- 95
  ggplot2::ggsave(output, res, width = 225, height = (length(res) - 1) * height_scaler,
                  units = "mm", dpi = dpi, limitsize = FALSE)
  fst::write_fst(data, paste0(output, ".fst"))
  message("Saved pshift plots to location: ", output)
  invisible(NULL)
}
