#' Try to shift 3 times, 1 strict shifting (full FFT), 2. less strict (wavelet),
#' 3. Hard 12 for allowed species
#' @noRd
shiftFootprintsByExperimentSafe <- function(df, shifting_table, accepted_lengths,
                                            allowed_hard12_species, BPPARAM) {
  fft_files <- fft_strength_files(df)
  fft_strengths <- names(fft_files)

  fft_strength <- fft_strengths[1]
  res <- tryCatch(shiftFootprintsByExperiment(df, output_format = "ofst",
                                              accepted.lengths = accepted_lengths,
                                              shift.list = shifting_table,
                                              BPPARAM = BPPARAM),
                  error = function(e) {
                    message(e)
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
                      message(e)
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
                          message(e)
                          message(name(df))
                          message("P-shifting failed also failed for hard 12nt,
                                    Fix manually (skipping to next project!)")
                          return(e)
                        })
      } else {
        message("Fix manually (skipping to next project!)")
        bad_pshifting_report(df)
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

shift_qc <- function(df) {
  # Plot max 39 libraries!
  subset <- if (nrow(df) >= 40) {seq(39)} else {seq(nrow(df))}
  invisible(shiftPlots(df[subset,], output = "auto", plot.ext = ".png"))
  # Check frame usage
  frameQC <- orfFrameDistributions(df)
  remove.experiments(df)
  zero_frame <- frameQC[frame == 0 & best_frame == FALSE,]
  # Store a flag that says good / bad shifting
  QCFolder <- QCfolder(df)
  data.table::fwrite(frameQC, file = file.path(QCFolder, "Ribo_frames_all.csv"))
  data.table::fwrite(zero_frame, file = file.path(QCFolder, "Ribo_frames_badzero.csv"))
  status_flag_files <- file.path(QCFolder, paste0(c("warning", "good"), ".rds"))
  names(status_flag_files) <- c("warning", "good")
  status <- "good"
  any_wrong_frame <- any(zero_frame$percent_length < 25) | nrow(zero_frame) == 0
  if (any_wrong_frame) {
    warning("Some libraries contain shift that is < 25% of CDS coverage")
    status <- "warning"
  }
  saveRDS(TRUE, status_flag_files[status])
  suppressWarnings(file.remove(status_flag_files[names(status_flag_files) != status]))
}

orfFrameDistributions <- function(df, type = "pshifted", weight = "score",
                                  orfs = loadRegion(df, part = "cds"),
                                  libraries = outputLibs(df, type = type, output.mode = "envirlist"),
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
