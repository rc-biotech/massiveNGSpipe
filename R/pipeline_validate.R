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
                            return_progress_vector = FALSE,
                            check_merged_org = FALSE) {
  n_bioprojects <- sum(unlist(lapply(pipelines, function(p) length(p$organisms))))
  steps <- names(config[["flag"]])
  negative_message <- steps; names(negative_message) <- steps
  negative_message[c("start","fetch","trim")] <- c("started","downloaded", "trimmed")
  negative_message <- paste(" - Not", negative_message);
  names(negative_message) <- steps
  index <- 1; done <- 0
  progress_index <- projects <- all_organism <- c()
  dt <- dt.trim <- data.table()
  for (pipeline in pipelines) {
    for (organism in names(pipeline$organisms)) {
      all_organism <- c(all_organism, organism)
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
      if (!file.exists(alignment.stats))
        alignment.stats <- file.path(out.aligned, "full_process.csv")

      dt <- rbindlist(list(dt, fread(alignment.stats, header = TRUE)))
      dt.trim <- rbindlist(list(dt.trim, ORFik:::trimming.table(trimmed.out)))
    }

    done <- done + 1
    }
  }
  if (check_merged_org) {
    organism_report(all_organism, config, progress_index)
  }
  cat("Number of studies completed\n")
  cat(done, " / ", n_bioprojects, "\n")
  cat("Last update:")
  cat("(", round(last_update_diff(config, units = "hours"), 1), " hours ago): ",
      format(last_update(config),usetz=TRUE), "\n", sep = "")
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

organism_report <- function(all_organism, config, progress_index) {
  all_organism_unique <- sort(unique(all_organism))
  cat("--------------------\n")
  cat("Processed organisms (Number of studies done/total, merged files done or not)\n")
  file_name <- sapply(all_organism_unique, organism_merged_exp_name)
  progress_index_done <- progress_index == length(config$flag)
  studies_per_org <- table(all_organism)
  studies_per_org_done <- table(all_organism[progress_index_done])
  studies_per_org_done <- studies_per_org_done[all_organism_unique]
  studies_per_org_done[is.na(studies_per_org_done)] <- 0
  names(studies_per_org_done) <- all_organism_unique
  merged_org_exists <-
    file.exists(file.path(config$config["exp"], paste0(file_name, ".csv")))
  names(merged_org_exists) <- all_organism_unique
  print(paste0(names(merged_org_exists), " (",studies_per_org_done,
               "/", studies_per_org,
               ",", merged_org_exists, ")"))
  return(invisible(NULL))
}

