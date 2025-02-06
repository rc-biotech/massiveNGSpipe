#' Validate pshifting
#'
#' It will do these two things: Plot all TIS regions and
#' make a frame distribution table for all libs.
#' It then save a rds file called either 'warning.rds' or 'good.rds'.
#' Depending on the results was accepted or not:
#' \code{any(zero_frame$percent_length < 25)}
#' @param df_list a list of ORFik experiments
#' @param config a pipeline config
#' @return invisible(NULL)
pipeline_validate_shifts <- function(df_list, config) {
  #
  # TODO: This function can be improved, especially how we detect bad shifts!
  for (df in df_list) {
    if (!step_is_next_not_done(config, "valid_pshift", name(df))) next
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
    if (any(zero_frame$percent_length < 25)) {
      warning("Some libraries contain shift that is < 25% of CDS coverage")
      saveRDS(FALSE, file.path(QCFolder, "warning.rds"))
    } else saveRDS(TRUE, file.path(QCFolder, "good.rds"))
    set_flag(config, "valid_pshift", name(df))
  }
}
