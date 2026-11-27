#' Extract every organism's experiment name from a pipelines list
#' @param pipelines the pipelines list (see `pipeline_init_all()`)
#' @return a flat (unlisted) list/vector of experiment id strings
#' @noRd
get_experiment_names <- function(pipelines) {
  exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
  exp <- unlist(exp, recursive = FALSE)
}

#' Experiment name for an organism's pooled-merge track
#' @param organisms character vector of scientific names
#' @return character, `"all_merged-<organism, spaces to underscores>"`
#' @noRd
organism_merged_exp_name <- function(organisms) {
  paste0("all_merged-", gsub(" ", "_", organisms))
}
#' Experiment name for an organism's pooled multi-modality track
#' @inheritParams organism_merged_exp_name
#' @return character, `organism_merged_exp_name()` with `"_modalities"` appended
#' @noRd
organism_merged_modalities_exp_name <- function(organisms) {
  paste0(organism_merged_exp_name(organisms), "_modalities")
}

#' Experiment name for an organism's sample collection
#' @inheritParams organism_merged_exp_name
#' @return character, `"all_samples-<organism, spaces to underscores>"`
#' @noRd
organism_collection_exp_name <- function(organisms) {
  paste0("all_samples-", gsub(" ", "_", organisms))
}

#' Organism names that always resolve to a specific strain key
#'
#' The default Ensembl Bacteria lookup takes the first strain of a species, which for "Zymomonas mobilis" is ATCC 10988,
#' while all our Zymomonas data is ZM4 (ATCC 31821). The ATCC 10988 reference was retired 2026-11-27 (Hakon), so a plain
#' "Zymomonas mobilis" label (any case) is mapped to "Zymomonas mobilis ZM4" before a reference is looked up or built.
#' Names that are not listed are returned unchanged.
#' @param organism character, scientific name(s)
#' @return character, same length
#' @noRd
organism_alias <- function(organism) {
  aliases <- c("zymomonas mobilis" = "Zymomonas mobilis ZM4")
  hit <- match(gsub("\\s+", " ", trimws(tolower(organism))), names(aliases))
  organism[!is.na(hit)] <- unname(aliases[hit[!is.na(hit)]])
  organism
}

#' Reference-folder-safe organism name
#' @param organism character, scientific name
#' @param full logical, default FALSE. If TRUE, join onto `full_path`.
#' @param full_path character, base reference directory, only evaluated
#' (and only needed) when `full = TRUE`
#' @return character: `organism` lowercased/trimmed/space-to-underscore,
#' or that joined onto `full_path` if `full = TRUE`
#' @noRd
reference_folder_name <- function(organism, full = FALSE,
                                  full_path = ORFik::config()["ref"]) {
  organism <- organism_alias(organism)
  res <- gsub(" ", "_", trimws(tolower(organism)))
  if (full) res <- file.path(full_path, res)
  return(res)
}
