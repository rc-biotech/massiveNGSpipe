#' If pshifting fails, Check qc and shift plots
#'
#' Also classifies the likely failure cause (see
#' \code{\link{classify_bad_shift_cause}}) into \code{qc_diagnostics.rds}
#' every call, independent of the \code{bad_shift_report_done.rds} guard
#' below -- that guard exists only to skip the expensive
#' \code{QCreport()}/\code{shiftPlots()} calls on repeat runs, so
#' classification (cheap: just reads already-written files) stays current
#' as upstream trim/align diagnostics evolve even for a study stuck here
#' across many pipeline runs.
#' @param df an ORFik experiment object
#' @param max_no_adapter_removed_pct numeric, percent, see
#' \code{\link{check_adapter_barcode_quality}} (default 80, matching
#' \code{pipeline_config()}'s own default)
#' @return invisible(NULL)
bad_pshifting_report <- function(df, max_no_adapter_removed_pct = 80) {
  message("Running failed pshift report")

  bam_dir <- bam_dir_from_df(df)
  trimmed_dir <- file.path(bam_dir, "trim")
  check_adapter_barcode_quality(bam_dir, trimmed_dir, max_no_adapter_removed_pct)
  cause <- classify_bad_shift_cause(bam_dir)
  message("- Likely cause: ", cause)

  report_flag_path <- file.path(QCfolder(df), "bad_shift_report_done.rds")
  report_exists <- file.exists(report_flag_path)
  if (report_exists) return(invisible(NULL))

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
