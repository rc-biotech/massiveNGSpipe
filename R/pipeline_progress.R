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
#' return integer of current flag step per experiment, set
#' named_progress_vector to TRUE, to get experiment names as names.
#' @param check_merged_org logical, default FALSE. Check status for merged organisms file
#' relative to pipeline.
#' @param named_progress_vector logical, default FALSE. Add experiment names as names
#' to progress_vector if returned, ignored otherwise.
#' @param system_usage_stats logical, default TRUE. Show system usage stats.
#' @param show_status_per_exp logical, default TRUE.
#'  Show a message per experiment with progression status.
#' @return invisible(NULL) / or progress vector
#' @importFrom plotly plot_ly layout
#' @export
progress_report <- function(pipelines, config,
                            show_status_per_exp = TRUE,
                            show_done = TRUE, status_plot = FALSE,
                            show_stats = FALSE,
                            return_progress_vector = FALSE,
                            check_merged_org = FALSE,
                            named_progress_vector = FALSE,
                            system_usage_stats = TRUE) {
  message("- Progress report:")

  status_per_study_list <- status_per_study(pipelines, config,
                                            show_status_per_exp, show_done,
                                            show_stats)
  progress_index <- status_per_study_list$progress_index
  n_bioprojects <- status_per_study_list$n_bioprojects
  if (check_merged_org) {
    organism_report(all_organism, config, progress_index)
  }
  if (system_usage_stats) {
    usage <- get_system_usage(one_liner = TRUE)
  }
  cat(status_per_study_list$message, "\n")

  ret <- invisible(NULL)
  if (return_progress_vector) {
    ret <- progress_index
    if (named_progress_vector) {
      names(ret) <- pipelines_names(pipelines)
    }
  }
  if (status_plot) {
    ret <- status_plot(status_per_study_list)
  }
  if (show_stats) {
    save_report(status_per_study_list)
  } else return(ret)
}

#' Get running processes
#'
#' Nice to kill a processes group by pgid
#' @param fields "pid,pgid,user,pcpu,pmem,nice,priority,time,command"
#' @param intern logical, default FALSE. If TRUE, capture as character vector.
#' @return either NULL or captured output
#' @export
get_running_processes <- function(fields = "pid,pgid,user,pcpu,pmem,nice,priority,time,command",
                                  intern = FALSE) {
  system(paste("ps xo ", fields, "--sort=pcpu,size"), intern = intern)
}

get_pgid <- function() as.integer(system("ps -o pgid= -p $$", intern = TRUE))

get_pgid_all <- function() {
  as.integer(get_running_processes("pgid", TRUE)[-1])
}

session_info_save <- function(config, pipelines) {
  if (is.null(config$session_dir)) stop("Can only save session info for active session!")
  info <- list(pipeline_names = names(pipelines),
               pgid = get_pgid(),
               init_time = config$init_time,
               status = "started")

  saveRDS(info, file.path(config$session_dir, "session_info.rds"))
  return(invisible(NULL))
}

session_info_read <- function(config) {
  if (is.null(config$session_dir)) stop("Can only save session info for active session!")
  return(readRDS(file.path(config$session_dir, "session_info.rds")))
}

#' Is the last call to run_pipeline still running / zombie process
#' @param config a pipeline config
#' @return logical, TRUE / FALSE
#' @export
last_session_stil_active <- function(config) {
  session_info_path <- list.files(sort(list.dirs(file.path(config$project, "session_logs")), decreasing = TRUE)[1],
                                  full.names = TRUE)
  if (length(session_info_path) == 0) stop("Last session has no valid info file")
  session_info <- readRDS(session_info_path)

  return(session_info$pgid %in% get_pgid_all())
}

#' Status per study
#' @inheritParams progress_report
#' @return a list structured as: list(done, progress_index, projects,
#'  alignment.stats.all, trimmed.out.all, n_bioprojects,
#'  report_dir = config$project, message)
#' @export
status_per_study <- function(pipelines, config, show_status_per_exp, show_done,
                             show_stats) {
  n_bioprojects <- sum(unlist(lapply(pipelines, function(p) length(p$organisms))))
  steps <- names(config[["flag"]])
  negative_message <- steps; names(negative_message) <- steps
  negative_message[c("start","fetch","trim")] <- c("started","downloaded", "trimmed")
  negative_message <- paste(" - Not", negative_message);
  names(negative_message) <- steps
  index <- 1; done <- 0
  progress_index <- projects <- all_organism <- c()
  alignment.stats.all <- trimmed.out.all <- c()
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
          if (show_status_per_exp) message(bio.index, project, negative_message[step])
          go_to_next <- TRUE
          break
        }
        progress_index_this <- progress_index_this + 1
      }
      progress_index <- c(progress_index, progress_index_this)
      if (go_to_next) next
      if (show_done & show_status_per_exp) message(bio.index, project, " - Done")
      if (show_stats) {
        out.aligned <- conf["bam"] #TODO, fix when it works to subset by name!
        trimmed.out <- file.path(out.aligned, "trim")
        alignment.stats <- file.path(out.aligned, "full_process_SINGLE.csv")

        if (!file.exists(alignment.stats))
          alignment.stats <- file.path(out.aligned, "full_process.csv")
        alignment.stats.all <- c(alignment.stats.all, alignment.stats)
        trimmed.out.all <- c(trimmed.out.all, trimmed.out)
      }
      done <- done + 1
    }
  }
  message <- paste(capture.output({
    cat("Number of studies completed\n")
    percentage_done <- round(100*(done / n_bioprojects), 2)
    percentage_done <- paste0(" (", percentage_done, "%)")
    cat(done, " / ", n_bioprojects, percentage_done, "\n")
    print(current_step_table(progress_index, steps))
    last_update_message(config)
  }), collapse = "\n")

  return(list(done = done, progress_index = progress_index,
              projects = projects, alignment.stats.all = alignment.stats.all,
              trimmed.out.all = trimmed.out.all, n_bioprojects = n_bioprojects,
              report_dir = config$project, message= message))
}

#' Store total trimming and alignment stats for pipeline
#' @param status_per_study_list a list, the output of \link{status_per_study}
#' @return invisible(NULL)
#' @noRd
save_report <- function(status_per_study_list) {

  alignment.stats.all <- status_per_study_list$alignment.stats.all
  trimmed.out.all <- status_per_study_list$trimmed.out.all
  done <- status_per_study_list$done
  total <- status_per_study_list$n_bioprojects
  report_dir <- status_per_study_list$report_dir

  message("- Loading alignment stats for all studies..")
  dt <- rbindlist(lapply(alignment.stats.all, function(f) {
    alignment.stats.dt <- try(fread(f, header = TRUE), silent = TRUE)
    if (is(alignment.stats.dt, "try-error")) alignment.stats.dt <- data.table()
    return(alignment.stats.dt)
  }), fill = TRUE)
  message("- Loading trimming stats for all studies..")
  dt.trim <- rbindlist(lapply(trimmed.out.all, function(f) {
    dt.trim.single <- try(ORFik:::trimming.table(f), silent = TRUE)
    if (is(dt.trim.single, "try-error")) dt.trim.single <- data.table()
    return(dt.trim.single)
  }), fill = TRUE)

  if (file.exists(file.path(report_dir, "FINAL_LIST.csv"))) {
    valid_runs <- fread(file.path(report_dir, "FINAL_LIST.csv"))
    dt_temp <- dt[gsub("collapsed_trimmed_|trimmed_|trimmed2_", "", sample) %in% valid_runs$Run]
    if (nrow(dt.trim) > 0) {
      if (nrow(dt_temp) > 0) dt <- dt_temp
      dt_temp <- dt.trim[raw_library %in% valid_runs$Run]
      if (nrow(dt_temp) > 0) dt.trim <- dt_temp
    }
  }
  report_dir <- file.path(report_dir, "summary_statistics")
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

  no_tables_made <- (nrow(dt) == 0 & nrow(dt.trim) == 0)
  if (done == 0) {
    cat("Nothing done, Returning without creating summary")
    return(invisible(NULL))
  }
  if (no_tables_made) {
    cat("No tables to plot, Returning without creating summary")
    return(invisible(NULL))
  }
  cat("Number of samples (total)", "\n")
  cat(nrow(dt), "\n")
  cat("-- Total Raw reads (in Billions):", "\n")
  cat(round(sum(dt.trim$raw_reads) / 1e9, 2), "\n")
  cat("-- Total Trimmed reads (in Billions):", "\n")
  cat(round(sum(dt.trim$trim_reads) / 1e9, 2), "\n")
  cat("-- Total mapped reads (in Billions):", "\n")
  cat("Read length distribution: ")
  print(summary(dt.trim$trim_mean_length))

  total_mapped_reads <- if (!is.null(dt$`total mapped reads #-genome`)) {
    dt$`total mapped reads #-genome`
  } else {dt$`total mapped reads #`}
  cat(round(sum(total_mapped_reads) / 1e9, 2), "\n")
  read.dist <- table(dt.trim$trim_mean_length)
  cat("-- Total read length usage (in percentages (bottom))", "\n")
  suppressWarnings(print(round((read.dist / sum(read.dist)) * 100, 1)))

  # Save
  fwrite(dt.trim, file.path(report_dir, "raw_trimmed_reads_stats.csv"))
  fwrite(dt, file.path(report_dir, "aligned_reads_stats.csv"))
  mapping_rate_plot(dt, dt.trim, report_dir)

  return(invisible(NULL))
}

current_step_table <- function(progress_index, steps) {
  all_flag_steps <- c(steps, "complete")

  tab <- table(progress_index)
  old_index <- names(tab)
  cat("Current steps active: name(flag_index):\n")
  names(tab) <- paste0(all_flag_steps[as.numeric(old_index) + 1], "(", old_index, ")")
  return(tab)
}

#' Mapping rate plot
#'
#' @param dt a data.table of alignment stats
#' @param dt.trim a data.table of trim stats
#' @param report_dir character path to directory
#' @return invisible(NULL)
#' @noRd
#' @import ggplot2 scales
mapping_rate_plot <- function(dt, dt.trim, report_dir) {
  mapping_rates_all <- dt$`total mapped reads %-genome`
  if (is.null(mapping_rates_all)) mapping_rates_all <- dt$`total mapped reads %`
  if (!is.null(mapping_rates_all) & nrow(dt) > 0) {
    dt[, mapping_rates := mapping_rates_all]

    stats_all_libraries_plot(ggplot(dt, aes(x = mapping_rates)),
                             file.path(report_dir, "genome_alignment_rates.png"),
                             xlab = "% mapped reads",
                             title = "Genome alignment rate % (all libraries)")
  }
  if (nrow(dt.trim) > 0) {
    stats_all_libraries_plot(ggplot(dt.trim, aes(x = `% trimmed`)),
                             file.path(report_dir, "trim_rates.png"))
    stats_all_libraries_plot(ggplot(dt.trim, aes(x = `trim_mean_length`)),
                             file.path(report_dir, "read_lengths_processed.png"),
                             xlab = "Read length",
                             title = "Processed read lengths (all libraries)")
  }



  return(invisible(NULL))
}

stats_all_libraries_plot <- function(ggplot, filename, xlab = "% reads removed by trimming",
                                                     n.breaks = 10, title = "Raw reads trimmed out % (all libraries)",
                                                     subtitle = "With flipped box-plot for quantiles",
                                                     ylim = c(0, 0.5), width = 7, height = 5) {
  plot <- ggplot +
    geom_histogram(bins = 20,
                   aes(y = after_stat(width*density))) +
    geom_boxplot(alpha = 0.1, color = "darkblue", outliers = FALSE) +
    ylab("% of libraries") + xlab(xlab) +
    scale_y_continuous(labels = percent_format()) +
    ggtitle(title, subtitle) +
    coord_cartesian(ylim = ylim) + theme_classic()
  if (!is.null(n.breaks)) {
    plot <- plot + scale_x_continuous(n.breaks = n.breaks)
  }
  plot(plot)
  ggsave(filename = filename, plot, width = width, height = height)
}

#' Status plot of pipeline
#'
#' Given all studies, report how far they have come as a heatmap (y-axis
#' is steps, x-axis is studies)
status_plot <- function(status_per_study_list) {

  steps <- status_per_study_list$steps
  progress_index <- status_per_study_list$progress_index
  projects <- status_per_study_list$projects
  n_bioprojects <- status_per_study_list$n_bioprojects
  done <- status_per_study_list$done

  status_matrix <- lapply(seq(length(steps)) - 1, function(x) as.data.table(matrix(as.integer(x <= progress_index), nrow = 1, byrow = TRUE)))
  status_matrix <- as.matrix(rbindlist(status_matrix))
  fig1 <- plot_ly(z = status_matrix, x = projects, y = steps, type = "heatmap", colors = c("red", "green")[unique(as.vector(status_matrix)) + 1], showscale=FALSE) %>%
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

