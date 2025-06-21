#' If pshifting fails, Check qc and shift plots
#'
#' @param df an ORFik experiment object
#' @return invisible(NULL)
bad_pshifting_report <- function(df) {
  message("Running failed pshift report")
  report_flag_path <- file.path(QCfolder(df), "bad_shift_report_done.rds")
  report_exists <- file.exists(report_flag_path)
  if (report_exists) return(invisible(NULL))

  lapply()
  envExp(df) <- new.env()
  out_dir <- file.path(QCfolder(df), "before_pshifting")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  barplot_path <- file.path(out_dir, "before_pshift_barplot")
  try(QCreport(df, out_dir, complex.correlation.plots = FALSE, create.ofst = FALSE))
  try(invisible(shiftPlots(df[seq(min(nrow(df), 39)),], output = barplot_path,
                       plot.ext = ".png")))
  saveRDS(TRUE, report_flag_path)
  return(invisible(NULL))
}
