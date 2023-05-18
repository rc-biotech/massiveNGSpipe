#' If pshifting fails, Check qc and shift plots
bad_pshifting_report <- function(df) {
  report_flag_path <- file.path(QCfolder(df), "bad_shift_report_done.rds")
  report_exists <- file.exists(report_flag_path)
  if (report_exists) return(invisible(NULL))

  envExp(df) <- new.env()
  try(QCreport(df, complex.correlation.plots = FALSE, create.ofst = FALSE))
  try(invisible(shiftPlots(df[seq(min(nrow(df), 39)),], output = "auto",
                       plot.ext = ".png")))
  saveRDS(TRUE, report_flag_path)
  return(invisible(NULL))
}
