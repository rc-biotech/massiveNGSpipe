#' Validate pshifting
#'
#' It will do these two things: Plot all TIS regions and
#' make a frame distribution table for all libs.
#' It then save a rds file called either 'warning.rds' or 'good.rds'.
#' Depending on the results was accepted or not:
#' \code{any(zero_frame$percent_length < 25)}
#' @param df_list a list of ORFik experiments
#' @param config a pipeline config
#' @return invisible(NULL)
pipeline_validate_shifts <- function(df_list, config) {
  #
  # TODO: This function can be improved, especially how we detect bad shifts!
  for (df in df_list) {
    if (!step_is_next_not_done(config, "valid_pshift", name(df))) next
    # Plot max 39 libraries!
    subset <- if (nrow(df) >= 40) {seq(39)} else {seq(nrow(df))}
    invisible(shiftPlots(df[subset,], output = "auto", plot.ext = ".png"))
    # Check frame usage
    frameQC <- orfFrameDistributions(df)
    remove.experiments(df)
    zero_frame <- frameQC[frame == 0 & best_frame == FALSE,]
    # Store a flag that says good / bad shifting
    QCFolder <- QCfolder(df)
    data.table::fwrite(frameQC, file = file.path(QCFolder, "Ribo_frames_all.csv"))
    data.table::fwrite(zero_frame, file = file.path(QCFolder, "Ribo_frames_badzero.csv"))
    if (any(zero_frame$percent_length < 25)) {
      warning("Some libraries contain shift that is < 25% of CDS coverage")
      saveRDS(FALSE, file.path(QCFolder, "warning.rds"))
    } else saveRDS(TRUE, file.path(QCFolder, "good.rds"))
    set_flag(config, "valid_pshift", name(df))
  }
}

#' Report summary of current progress
#'
#' How many studies have finished etc
#' @inheritParams run_pipeline
#' @param show_stats logical, default TRUE, output trim/alignment stats
#' + plots, set to FALSE if you only want progress report.
#' @param show_done logical, default TRUE. If FALSE, display only status
#' of projects that are not done. Stats will still show for all.
#' @param status_plot plot an  interactive plot of total status
#' @param return_progress_vector logical, default FALSE. If true,
#' return progress vector.
#' @return invisible(NULL) / or progress vector
#' @importFrom plotly plot_ly layout
#' @export
progress_report <- function(pipelines, config, show_stats = FALSE,
                            show_done = TRUE, status_plot = FALSE,
                            return_progress_vector = FALSE) {
  n_bioprojects <- sum(unlist(lapply(pipelines, function(p) length(p$organisms))))
  steps <- names(config[["flag"]])
  negative_message <- steps; names(negative_message) <- steps
  negative_message[c("start","fetch","trim")] <- c("started","downloaded", "trimmed")
  negative_message <- paste(" - Not", negative_message);
  names(negative_message) <- steps
  index <- 1; done <- 0
  progress_index <- projects <- c()
  dt <- dt.trim <- data.table()
  for (pipeline in pipelines) {
    for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    project <- conf["exp"]
    projects <- c(projects, project)
    bio.index <- paste0("(", index, "/",n_bioprojects,") ")
    index <- index + 1
    go_to_next <- FALSE
    progress_index_this <- 0
    for (step in steps) {
      if (!step_is_done(config, step, project)) {
        message(bio.index, project, negative_message[step])
        go_to_next <- TRUE
        break
      }
      progress_index_this <- progress_index_this + 1
    }
    progress_index <- c(progress_index, progress_index_this)
    if (go_to_next) next
    if (show_done) message(bio.index, project, " - Done")
    if (show_stats) {
      out.aligned <- conf["bam"] #TODO, fix when it works to subset by name!
      trimmed.out <- file.path(out.aligned, "trim")
      alignment.stats <- file.path(out.aligned, "full_process_SINGLE.csv")
      dt <- rbindlist(list(dt, fread(alignment.stats, header = TRUE)))
      dt.trim <- rbindlist(list(dt.trim, ORFik:::trimming.table(trimmed.out)))
    }

    done <- done + 1
    }
  }
  cat("Number of studies completed\n")
  cat(done, " / ", n_bioprojects, "\n")
  ret <- invisible(NULL)
  if (return_progress_vector) ret <- progress_index
  if (status_plot) {
    ret <- status_plot(steps, progress_index, projects, n_bioprojects, done)
  }
  if (show_stats) {
    save_report(dt, dt.trim, done, total = n_bioprojects, config$project)
  } else return(ret)
}

save_report <- function(dt, dt.trim, done, total, report_dir) {
  if (done == 0) {
    cat("Nothing done, Returning without creating summary")
    return(invisible(NULL))
  }
  cat("Number of samples (total)", "\n")
  cat(nrow(dt), "\n")
  cat("-- Total Raw reads (in Billions):", "\n")
  cat(round(sum(dt.trim$raw_reads) / 1e9, 2), "\n")
  cat("-- Total Trimmed reads (in Billions):", "\n")
  cat(round(sum(dt.trim$trim_reads) / 1e9, 2), "\n")
  cat("-- Total mapped reads (in Billions):", "\n")
  cat(round(sum(dt$`total mapped reads #-genome`) / 1e9, 2), "\n")
  #print(dt.trim)
  #print(nrow(dt.trim))
  read.dist <- table(dt.trim$trim_mean_length)
  #read.dist <- read.dist[read.dist > 50]
  cat("-- Total read length usage (in percentages (bottom))", "\n")
  suppressWarnings(print(round((read.dist / sum(read.dist)) * 100, 1)))
  #print(dt)
  # Save
  fwrite(dt.trim, file.path(report_dir, "raw_trimmed_reads_stats.csv"))
  library(ggplot2); library(scales)
  plot <- ggplot(dt, aes(x = dt$`total mapped reads %-genome`)) +
    geom_histogram(bins = 20,
                   aes(y = after_stat(width*density))) +
    geom_boxplot(alpha = 0.1, color = "darkblue") +
    ylab("% of libraries") + xlab("% mapped reads") +
    scale_y_continuous(labels = percent_format()) +
    scale_x_continuous(n.breaks = 10) +
    ggtitle("Ribo-seq mapped read % (all libraries)",
            "With flipped box-plot for quantiles") +
    coord_cartesian(ylim = c(0, 0.5)) + theme_classic()
  plot(plot)
  ggsave(filename = file.path(report_dir, "genome_alignment_rate.png"),
         plot, width = 5, height = 5)
  return(invisible(NULL))
}

#' Status plot of pipeline
#'
#' Given all studies, report how far they have come as a heatmap (y-axis
#' is steps, x-axis is studies)
status_plot <- function(steps, progress_index, projects, n_bioprojects, done) {
  status_matrix <- lapply(seq(length(steps)) - 1, function(x) as.data.table(matrix(as.integer(x <= progress_index), nrow = 1, byrow = TRUE)))
  status_matrix <- as.matrix(rbindlist(status_matrix))
  fig1 <- plot_ly(z = status_matrix, x = projects, y = steps, type = "heatmap", colors = c("red", "green"), showscale=FALSE) %>%
    layout(title = list(text = paste0('NGS Processing Status',
                                      '<br>',
                                      '<sup>',
                                      paste0("Completed: ", done, " / ", n_bioprojects),
                                      '</sup>')), margin=list(t = 75))
  return(fig1)
}

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
