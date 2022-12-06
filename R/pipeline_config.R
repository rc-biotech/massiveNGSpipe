pipeline_config <- function(project_dir, config = ORFik::config()) {
  metadata <- paste0(project_dir, "/metadata")
  flags <- pipeline_flags(project_dir)
  return(list(project = project_dir, config = config, flag = flags, metadata = metadata))
}