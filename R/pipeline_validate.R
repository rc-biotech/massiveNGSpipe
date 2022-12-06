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

progress_report <- function(pipelines, config) {
  n_bioprojects <- sum(unlist(lapply(pipelines, function(p) length(p$organisms))))
  steps <- names(config[["flag"]])
  negative_message <- steps; names(negative_message) <- steps
  negative_message[c("fetch","trim")] <- c("started", "trimmed")
  negative_message <- paste(" - Not", negative_message); 
  names(negative_message) <- steps
  index <- 1; done <- 0
  dt <- dt.trim <- data.table()
  for (pipeline in pipelines) {
    for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    project <- conf["exp"]
    bio.index <- paste0("(", index, "/",n_bioprojects,") ")
    index <- index + 1
    go_to_next <- FALSE
    for (step in steps) {
      if (!step_is_done(config, step, project)) {
        message(bio.index, project, negative_message[step])
        go_to_next <- TRUE
        break
      }
    }
    if (go_to_next) next
    message(bio.index, project, " - Done")
    out.aligned <- conf["bam"] #TODO, fix when it works to subset by name!
    trimmed.out <- file.path(out.aligned, "trim")
    alignment.stats <- file.path(out.aligned, "full_process_SINGLE.csv")
    dt <- rbindlist(list(dt, fread(alignment.stats, header = TRUE)))
    dt.trim <- rbindlist(list(dt.trim, ORFik:::trimming.table(trimmed.out)))
    done <- done + 1
    }
  }
  save_report(dt, dt.trim, done, total = n_bioprojects, config$project)
}

save_report <- function(dt, dt.trim, done, total, report_dir) {
  cat("Number of studies completed initial: download, alignment, pshifting and merge!\n")
  cat(done, " / ", total, "\n")
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