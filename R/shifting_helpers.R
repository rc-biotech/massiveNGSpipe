#' Try to shift 3 times, 1 strict shifting (full FFT), 2. less strict (wavelet),
#' 3. Hard 12 for allowed species
#' @noRd
shiftFootprintsByExperimentSafe <- function(df, shifting_table, accepted_lengths,
                                            allowed_hard12_species, BPPARAM,
                                            max_no_adapter_removed_pct = 80) {
  fft_files <- fft_strength_files(df)
  fft_strengths <- names(fft_files)

  fft_strength <- fft_strengths[1]
  res <- tryCatch(shiftFootprintsByExperiment(df, output_format = "ofst",
                                              accepted.lengths = accepted_lengths,
                                              shift.list = shifting_table,
                                              BPPARAM = BPPARAM),
                  error = function(e) {
                    message(conditionMessage(e))
                    message(name(df))
                    message("P-shifting failed with strict FFT, trying weak!")
                    return(e)
                  })

  if(inherits(res, "error")) {
    fft_strength <- fft_strengths[2]
    res <- tryCatch(shiftFootprintsByExperiment(df, output_format = "ofst",
                                                accepted.lengths = accepted_lengths,
                                                shift.list = shifting_table,
                                                strict.fft = FALSE,
                                                BPPARAM = BPPARAM),
                    error = function(e) {
                      message(conditionMessage(e))
                      message(name(df))
                      message("P-shifting failed also failed with weak FFT")
                      return(e)
                    })
    if(inherits(res, "error")) {

      if (organism(df) %in% allowed_hard12_species) {
        fft_strength <- fft_strengths[3]
        shifting_table <- template_shift_table(df, accepted_lengths)
        res <- tryCatch(shiftFootprintsByExperiment(df, output_format = "ofst",
                                                    accepted.lengths = accepted_lengths,
                                                    shift.list = shifting_table,
                                                    strict.fft = FALSE,
                                                    BPPARAM = BPPARAM),
                        error = function(e) {
                          message(conditionMessage(e))
                          message(name(df))
                          message("P-shifting failed also failed for hard 12nt,
                                    Fix manually (skipping to next project!)")
                          return(e)
                        })
      } else {
        message("Fix manually (skipping to next project!)")
        bad_pshifting_report(df, max_no_adapter_removed_pct)
      }
    }
  }
  if(!inherits(res, "error")) {
    suppressWarnings(file.remove(fft_files[!(names(fft_files) %in% fft_strength)]))
    dir.create(QCfolder(df), showWarnings = FALSE)
    suppressWarnings(saveRDS(fft_strength, fft_files[fft_strength]))
  }

  return(res)
}

shifts_load_safe <- function(df, reuse_shifts_if_existing) {
  shifting_table <- NULL
  if (reuse_shifts_if_existing) {
    shifts <- suppressWarnings(try(shifts_load(df), silent = TRUE))
    if (!is(shifts, "try-error")) {
      if (length(shifts) > 0) {
        shifting_table <- shifts
        length_original <- length(shifts)
        rfp_files <- filepath(df, "ofst")
        if (!all(names(shifting_table) %in% rfp_files)) {
          hits_order <- sapply(runIDs(df), function(x) grep(x, names(shifting_table)))
          names(shifting_table) <- rfp_files[unlist(hits_order)]
        }
        stopifnot(length(shifting_table) == length_original)
        if (length(shifting_table) != nrow(df)) {
          shifting_table <- shifting_table[rfp_files]
          names(shifting_table) <- rfp_files
        }
        stopifnot(length(shifting_table) == nrow(df))
      }
    }
  }
  return(shifting_table)
}

fft_strength_files <- function(df) {
  fft_strengths <- c("strong", "weak", "manual_12")
  fft_files <- file.path(QCfolder(df), paste0(fft_strengths, "_periodicity.rds"))
  names(fft_files) <- fft_strengths
  return(fft_files)
}

shift_qc <- function(df, BPPARAM = bpparam(), max_no_adapter_removed_pct = 80) {
  # Plot max 39 libraries!
  subset <- if (nrow(df) >= 40) {seq(39)} else {seq(nrow(df))}
  has_leaders <- length(filterTranscripts(df, 5, 0, 0, stopOnEmpty = FALSE))
  upstream <- ifelse(has_leaders, 5, 0)
  message("- Shift barplots")
  invisible(shiftPlots(df[subset,], output = "auto", plot.ext = ".png",
                       upstream = upstream, BPPARAM = BPPARAM))
  # Check frame usage
  message("- Frame distributions")
  frameQC <- orfFrameDistributions(df, BPPARAM = BPPARAM)
  remove.experiments(df)
  zero_frame <- frameQC[frame == 0 & best_frame == FALSE,]

  QCFolder <- QCfolder(df)
  data.table::fwrite(frameQC, file = file.path(QCFolder, "Ribo_frames_all.csv"))
  data.table::fwrite(zero_frame, file = file.path(QCFolder, "Ribo_frames_badzero.csv"))

  # Adapter/barcode-quality diagnostic signal (see R/pshift_diagnostics.R):
  # recorded here, alongside the periodicity verdict, since this is the
  # first point after trim where a P-shift-failure classification is
  # actually attempted.
  trimmed_dir <- file.path(bam_dir_from_df(df), "trim")
  check_adapter_barcode_quality(bam_dir_from_df(df), trimmed_dir,
                                max_no_adapter_removed_pct)

  # Store a flag that says good / bad / no-data shifting.
  # nrow(frameQC) == 0 (zero CDS frame data at all, e.g. because almost
  # nothing aligned to real ORFs) is a different situation from
  # nrow(zero_frame) == 0 (perfect periodicity for every read length) --
  # the corrected has_cds_data tells these two apart.
  periodicity_check_flag(zero_frame, QCFolder, has_cds_data = nrow(frameQC) > 0)
}

orfFrameDistributions <- function(df, type = "pshifted", weight = "score",
                                  orfs = loadRegion(df, part = "cds"),
                                  libraries = outputLibs(df, type = type, output.mode = "envirlist", BPPARAM = BPPARAM),
                                  BPPARAM = BiocParallel::bpparam()) {
  frame_sum_per <- regionPerReadLengthPerLib(orfs, libraries, scoring = "frameSumPerL",
                                             weight, BPPARAM)
  if (nrow(frame_sum_per) == 0) frame_sum_per <- data.table(frame = numeric(),
                                                            fraction = numeric(),
                                                            length = numeric(),
                                                            score = numeric())
  frame_sum_per[, frame := as.factor(frame)]
  frame_sum_per[, fraction := as.factor(fraction)]
  frame_sum_per[, percent := (score / sum(score))*100, by = fraction]
  frame_sum_per[, percent_length := (score / sum(score))*100, by = .(fraction, length)]
  frame_sum_per[, best_frame := (percent_length / max(percent_length)) == 1, by = .(fraction, length)]
  frame_sum_per[, fraction := factor(fraction, levels = names(libraries),
                                     labels = gsub("^RFP_", "", names(libraries)), ordered = TRUE)]

  frame_sum_per[, fraction := factor(fraction, levels = unique(fraction), ordered = TRUE)]
  frame_sum_per[]
  return(frame_sum_per)
}

regionPerReadLengthPerLib <- function(grl, libraries, scoring = "frameSumPerL",
                         weight = "score", BPPARAM = BiocParallel::bpparam()) {
  stopifnot(is(libraries, "list"))
  # Frame distribution over all
  frame_sum_per1 <- bplapply(libraries, FUN = function(lib, grl, weight) {
    name <- attr(lib, "name_short")
    message("- ", name)
    total <- regionPerReadLength(grl, lib,
                                 withFrames = TRUE, scoring = "frameSumPerL",
                                 weight = weight, drop.zero.dt = TRUE,
                                 exclude.zero.cov.grl = length(lib) == 0)
    if (nrow(total) > 0) {
      total[, length := fraction]
      total[, fraction := rep(name, nrow(total))]
    }
    return(total)
  }, grl = grl, weight = weight, BPPARAM = BPPARAM)
  return(rbindlist(frame_sum_per1))
}

template_shift_table_exps <- function(exps, accepted.lengths = c(20, 21, 25:33)) {
  lapply(exps, function(exp) {
    message(exp)
    df <- read.experiment(exp)
    l <- template_shift_table(df, accepted.lengths = accepted.lengths)
    shifts_save(l, file.path(libFolder(df), "pshifted"))
  })
  return(invisible(TRUE))
}

#' Classify a P-shifted experiment's periodicity as good/warning/no_data
#'
#' \code{nrow(zero_frame) == 0} on its own is ambiguous: it's true both
#' when periodicity is perfect for every read length (no bad rows to
#' report) AND when there is zero CDS frame data at all (e.g. almost
#' nothing aligned to real ORFs, so there was never anything to compute
#' frame usage from in the first place) -- two very different situations.
#' \code{has_cds_data} disambiguates them into a third \code{"no_data"}
#' status, distinct from a genuine \code{"good"} verdict.
#' @param zero_frame data.table, the frame==0 & !best_frame rows from
#' \code{orfFrameDistributions()}'s output (see \code{shift_qc()})
#' @param QCFolder character, directory the status flag files are
#' written to (see \code{QCfolder()})
#' @param has_cds_data logical, default \code{nrow(zero_frame) > 0} --
#' this default deliberately reproduces the historical (buggy) boolean
#' for any caller other than \code{shift_qc()}, which passes the
#' corrected \code{nrow(frameQC) > 0} instead.
#' @return character, one of \code{"good"}, \code{"warning"}, \code{"no_data"}
#' @noRd
periodicity_check_flag <- function(zero_frame, QCFolder,
                                   has_cds_data = nrow(zero_frame) > 0) {
  status_flag_files <- file.path(QCFolder, paste0(c("warning", "good", "no_data"), ".rds"))
  names(status_flag_files) <- c("warning", "good", "no_data")
  if (!has_cds_data) {
    warning("No CDS frame data at all -- nothing aligned to real ORFs?")
    status <- "no_data"
  } else if (any(zero_frame$percent_length < 25)) {
    warning("Some libraries contain shift that is < 25% of CDS coverage")
    status <- "warning"
  } else {
    status <- "good"
  }
  saveRDS(TRUE, status_flag_files[status])
  suppressWarnings(file.remove(status_flag_files[names(status_flag_files) != status]))
  update_qc_diagnostics(file.path(QCFolder, "qc_diagnostics.rds"), periodicity_status = status)
  return(status)
}
