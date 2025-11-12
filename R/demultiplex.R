
#' Given ..., demultiplex to fastq
#' @param path path to directory of
#' @param csv path to csv file of
#' @return invisible(NULL), files saved to disc
demultiplex <- function(path, csv, demultiplex_tool_bin = "~/bin/demultiplexer") {
  call <- paste(demultiplex_tool_bin, "-i", path, "-c", csv)
  ret <- system(call)

  if (ret != 0) stop("Demultiplexing failed")
  return(invisible(NULL))
}
