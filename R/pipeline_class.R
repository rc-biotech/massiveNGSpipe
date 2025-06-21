#' Get study data.table
#' @export
study_pipeline <- function(pipelines) {
  rbindlist(lapply(pipelines, function(pipeline) pipeline$study))
}

#' Get study directories
#' @export
list_dirs_of_pipeline <- function(pipeline) {
  stopifnot(length(pipeline) == 3)
  conf <- lapply(pipeline$organisms, function(exp) {
    dirs <- exp$conf
    dirs <- c(dirs,
              trim = file.path(dirs["bam"], "trim"),
              aligned = file.path(dirs["bam"], "aligned"))
    dirs <- c(dirs,
              ofst = file.path(dirs["aligned"], "ofst"),
              pshifted = file.path(dirs["aligned"], "pshifted"),
              cov = file.path(dirs["aligned"], "cov_RLE"),
              bigwig = file.path(dirs["aligned"], "bigwig"))
    return(dirs)
  })
  return(conf)
}

#' Get files of study in dir type
#' @export
list_files_of_type <- function(pipeline, type = "fastq", full_names = FALSE) {
  all_dirs <- list_dirs_of_pipeline(pipeline)
  stopifnot(type %in% names(all_dirs[[1]]))
  list_files <- lapply(all_dirs, function(dirs) list.files(dirs[type], full.names = full_names))
  return(list_files)
}

#' Get missing files of study in dir type
#' @export
list_files_of_type_missing <- function(pipelines, type = "fastq") {
  format_regex <- "\\.fastq$|\\.fastq.gz$|\\.fasta\\.fasta.gz\\.bam|\\.ofst|_pshifted.*"
  study <- study_pipeline(pipelines)
  lapply(pipelines, function(pipeline)
    study[!(Run %in% gsub(format_regex, "", unlist(list_files_of_type(pipeline, type), use.names = FALSE)))])
}
