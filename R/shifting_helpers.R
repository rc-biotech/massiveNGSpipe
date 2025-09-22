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
  any_wrong_frame <- any(zero_frame$percent_length < 25)
  if (any_wrong_frame) {
    warning("Some libraries contain shift that is < 25% of CDS coverage")
    status <- "warning"
  }
  saveRDS(TRUE, status_flag_files[status])
  suppressWarnings(file.remove(status_flag_files[names(status_flag_files) != status]))
}
