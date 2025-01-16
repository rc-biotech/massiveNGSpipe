get_experiment_names <- function(pipelines) {
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
}

organism_merged_exp_name <- function(organisms) {
  paste0("all_merged-", gsub(" ", "_", organisms))
}
organism_merged_modalities_exp_name <- function(organisms) {
  paste0(organism_merged_exp_name(organisms), "_modalities")
}

organism_collection_exp_name <- function(organisms) {
  paste0("all_samples-", gsub(" ", "_", organisms))
}

reference_folder_name <- function(organism, full = FALSE,
                                  full_path = ORFik::config()["ref"]) {
  res <- gsub(" ", "_", trimws(tolower(organism)))
  if (full) res <- file.path(full_path, res)
  return(res)
}
