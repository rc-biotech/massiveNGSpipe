#' For docker usage copy
#'
#' If you have download docker and release docker.
#' This is useful to copy files from download to release
#' @param config a config object
#' @param pipelines pipeline objects, default:
#' \code{pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE)}
#' @param new_exp_dir character, default "~/livemount/Bio_data/ORFik_experiments/"
#' @param docker_conversion the sed conversion string: 's/livemount\///g'
#' @param individual_studies logical, default TRUE. Should raw ORFik experiments
#' be exported per study.
#' @param merged_by_organism logical, default TRUE.
#' Create a new study, merge the superset of all per organism into single file.
#' This requires that you have run: pipeline_merge_org
#' @param all_by_organism logical, default TRUE.
#' Create a new study, the superset of all studies of all per organism.
#' This requires that you have run: pipeline_merge_org
#' @param all_modalities_by_org logical, default TRUE.
#' Should the all modalities tracks be exported per organism.
#' @return logical, TRUE if sucessful for all.
#' @export
docker_copy_done_experiments <- function(config,
                                         pipelines = pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE),
                                         new_exp_dir = "~/livemount/Bio_data/ORFik_experiments/",
                                         docker_conversion = "'s/livemount\\///g'",
                                         individual_studies = TRUE,
                                         merged_by_organism = TRUE,
                                         all_by_organism = TRUE,
                                         all_modalities_by_org = TRUE) {
  if (length(pipelines) == 0) stop("Empty pipelines object given!")
  message("- Individual")
  done_stats <- progress_report(pipelines, config, return_progress_vector = TRUE)
  max_step <- length(unlist(config$flag_steps))
  done_all_steps <- done_stats == max_step
  res <- c()
  csv_names <- paste0(pipelines_names(pipelines), ".csv")
  csv_names <- csv_names[done_all_steps]
  old_exp_dir <- config$config["exp"]

  if (length(csv_names) == 0) {
    warning("No experiments are done, returning directly!")
    return(FALSE)
  }
  if (individual_studies) {
    res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
  }

  if (merged_by_organism | all_by_organism | all_modalities_by_org) {
    names(pipelines) <- NULL
    exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
    exp <- unlist(exp, recursive = FALSE)
    # Merge all per organism
    done_exp_list <- exp[done_all_steps]
    done_organisms <- unique(names(done_exp_list))

    if (merged_by_organism) {
      message("- Merged by organism")
      csv_names <- paste0(organism_merged_exp_name(done_organisms), ".csv")
      existing_exps <- list.files(config$config["exp"])
      merged_is_made <- csv_names %in% existing_exps
      csv_names <- csv_names[merged_is_made]
      res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
    }
    if (all_by_organism) {
      message("- All by organism")
      csv_names <- paste0(organism_collection_exp_name(done_organisms), ".csv")
      existing_exps <- list.files(config$config["exp"])
      merged_is_made <- csv_names %in% existing_exps
      csv_names <- csv_names[merged_is_made]

      res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
    }

    if (all_modalities_by_org) {
      message("- All modalities by organism")
      csv_names <- paste0(organism_merged_modalities_exp_name(done_organisms), ".csv")
      existing_exps <- list.files(config$config["exp"])
      merged_is_made <- csv_names %in% existing_exps
      csv_names <- csv_names[merged_is_made]

      res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
    }
  }
  return(all(res == 0))
}

#' Copy done experiment csvs to new directory
#' @examples
#' csvs <- list.files(ORFik::config()["exp"], pattern = "modalities\\.csv$")
#' copy_experiments_to(csvs)
copy_experiments_to <- function(csv_names, old_exp_dir = ORFik::config()["exp"],
                                new_exp_dir = sub("_local$|_local/$", "", old_exp_dir),
                                docker_conversion = "'s/livemount\\///g'") {
  done_exp_old_path <- file.path(old_exp_dir, csv_names)
  done_exp_new_path <- file.path(new_exp_dir, csv_names)
  res <- c()
  for (i in seq_along(csv_names)) {
    file.copy(done_exp_old_path[i], done_exp_new_path[i], overwrite = TRUE)
    cmd <- system(paste0("sed -i ", docker_conversion, " ",done_exp_new_path[i]))
    res <- c(res, cmd)
  }
  return(res)
}

#' Get name of function
#' @importFrom rlang ensyms
#' @importFrom purrr map
#' @param ... functions
#' @return character, name of function
#' @export
name_of_function <- function(...) unlist(purrr:::map(rlang::ensyms(...) , as.character), use.names = FALSE)

#' @inherit rstudioapi::filesPaneNavigate
#' @importFrom rstudioapi filesPaneNavigate
navigate <- rstudioapi::filesPaneNavigate


