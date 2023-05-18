#' For docker usage copy
#'
#' If you have download docker and release docker.
#' This is useful to copy files from download to release
#' @param config a config object
#' @param pipelines pipeline objects, default pipeline_init_all(config)
#' @param new_exp_dir character, default "~/livemount/Bio_data/ORFik_experiments/"
#' @param docker_conversion the sed conversion string: 's/livemount\///g'
#' @return logical, TRUE if sucessful for all.
#' @export
docker_copy_done_experiments <- function(config, pipelines = pipeline_init_all(config),
                                         new_exp_dir = "~/livemount/Bio_data/ORFik_experiments/",
                                         docker_conversion = "'s/livemount\\///g'",
                                         merged_by_organism = TRUE) {
  if (length(pipelines) == 0) stop("Empty pipelines object given!")
  done_stats <- progress_report(pipelines, config, return_progress_vector = TRUE)
  max_step <- length(unlist(config$flag_steps))
  done_all_steps <- done_stats == max_step
  csv_names <- paste0(pipelines_names(pipelines), ".csv")
  csv_names <- csv_names[done_all_steps]

  if (length(csv_names) == 0) {
    warning("No experiments are done, returning directly!")
    return(FALSE)
  }

  done_exp_old_path <- paste0(config$config["exp"], csv_names)
  done_exp_new_path <- paste0(new_exp_dir, csv_names)
  res <- c()
  for (i in seq_along(done_exp_old_path)) {
    file.copy(done_exp_old_path[i], done_exp_new_path[i], overwrite = TRUE)
    cmd <- system(paste0("sed -i ", docker_conversion, " ",done_exp_new_path[i]))
    res <- c(res, cmd)
  }


  if (merged_by_organism) {
    names(pipelines) <- NULL
    exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
    exp <- unlist(exp, recursive = FALSE)
    # Merge all per organism
    done_exp_list <- exp[done_all_steps]
    done_organisms <- unique(names(done_exp_list))
    csv_names <- paste0(organism_merged_exp_name(done_organisms), ".csv")

    existing_exps <- list.files(config$config["exp"])
    merged_is_made <- csv_names %in% existing_exps
    csv_names <- csv_names[merged_is_made]

    done_exp_old_path <- paste0(config$config["exp"], csv_names)
    done_exp_new_path <- paste0(new_exp_dir, csv_names)
    for (i in seq_along(csv_names)) {
      file.copy(done_exp_old_path[i], done_exp_new_path[i], overwrite = TRUE)
      cmd <- system(paste0("sed -i ", docker_conversion, " ",done_exp_new_path[i]))
      res <- c(res, cmd)
    }
  }
  return(all(res == 0))
}
