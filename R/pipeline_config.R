#' Initial path config
#'
#' Set up all paths and flag directories
#' @param project_dir where will specific pipeline outputs be put
#' @param config path, default \code{ORFik::config()}, where will
#' fastq, bam, references and ORFik experiments go
#' @return a list with a defined config
#' @export
pipeline_config <- function(project_dir = file.path(dirname(config)[1], "NGS_pipeline"),
                            config = ORFik::config()) {
  metadata <- paste0(project_dir, "/metadata")
  flags <- pipeline_flags(project_dir)
  return(list(project = project_dir, config = config, flag = flags, metadata = metadata))
}
