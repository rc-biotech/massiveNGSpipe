#' Analysis of Ribo-seq and RNA-seq
#' @param df.rfp an ORFik experiment
#' @param df.rna an ORFik experiment
#' @param canonical_isoforms a character vector with tx names,
#'  default: canonical_isoforms(df.rfp)
#' @return invisible(NULL), results saved to disc
#' @export
analysis_pipeline_ribo_rna <- function(df.rfp = read.experiment("Eleonora-homo_sapiens"),
                                       df.rna = read.experiment("Eleonora-homo_sapiens_RNA-seq"),
                                       canonical_isoforms = canonical_isoforms(df.rfp)) {
  df.rfp <- read.experiment("vaccine_mod_bnt_111024-homo_sapiens_x_sars_cov2")
  df.rna <- read.experiment("vaccine_mod_bnt_111024-homo_sapiens_x_sars_cov2_RNA-seq")
  message("- Ribo-seq & RNA-seq analysis_pipeline starting")

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Convert files and run Annotation vs alignment QC
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- QC..")
  # General QC
  ORFikQC(df.rfp, create.ofst = FALSE)
  ORFikQC(df.rna, create.ofst = FALSE)
  ggsave(filename = file.path(QCfolder(df.rfp), paste0("PCAplot_", df.rfp@experiment, ".png")),
         ORFik:::pcaExperiment(df.rfp), height = 7, width = 7)
  ggsave(filename = file.path(QCfolder(df.rna), paste0("PCAplot_", df.rna@experiment, ".png")),
         ORFik:::pcaExperiment(df.rna), height = 7, width = 7)
  remove.experiments(df.rfp) # Remove loaded data (it is not pshifted)
  remove.experiments(df.rna)
  RiboQC.plot(df.rfp, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
  remove.experiments(df.rfp)

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Create heatmaps (Ribo-seq)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- Ribo heatmaps")
  # Pre-pshifting
  heatMapRegion(df.rfp, region = c("TIS", "TTS"), shifting = "5prime", type = "ofst",
                outdir = file.path(QCfolder(df.rfp), "heatmaps/pre-pshift/"))
  remove.experiments(df.rfp)
  # After pshifting
  heatMapRegion(df.rfp, region = c("TIS", "TTS"), shifting = "5prime", type = "pshifted",
                outdir = file.path(QCfolder(df.rfp), "heatmaps/pshifted/"))


  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Differential translation analysis (condition: WT vs CO)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- Differential analysis (DESeq2 model)")
  # The design is by default chosen by this factor: The condition column in this case
  design(df.rfp, multi.factor = FALSE)
  # We now run, and here get 11 unique DTEG genes
  RFP_counts <- countTable(df.rfp, region = "cds", type = "summarized")
  RNA_counts <- countTable(df.rna, region = "cds", type = "summarized")
  RFP_counts_canonical <- RFP_counts[rownames(RFP_counts) %in% canonical_isoforms,]
  RNA_counts <- countTable(df.rna, region = "cds", type = "summarized")
  RNA_counts_canonical <- RNA_counts[rownames(RNA_counts) %in% canonical_isoforms,]

  res <- DTEG.analysis(df.rfp, df.rna, RFP_counts = RFP_counts_canonical, RNA_counts = RNA_counts_canonical)
  res_extended <- append_gene_symbols(res, symbols(df.rfp), extend_id = FALSE)
  fwrite(res_extended, file.path(QCfolder(df.rfp), "DTEG_analysis.csv"))
  res_extended <- fread(file.path(QCfolder(df.rfp), "DTEG_analysis.csv"))
  gg <- DTEG.plot(res_extended, QCfolder(df.rfp), width = 8, height = 6,
                  plot.ext = ".pdf", plot_to_console = FALSE)
  gg <- DTEG.plot(res_extended, QCfolder(df.rfp), width = 8, height = 6,
            plot.ext = ".pdf", plot_to_console = FALSE)
  gg <- DTEG.plot(res_extended, QCfolder(df.rfp), width = 8, height = 6,
                  plot.ext = ".png", plot_to_console = FALSE)
  gg <- DTEG.plot(res_extended, QCfolder(df.rfp), width = 8, height = 6,
                  plot.ext = ".jpg", plot_to_console = FALSE)
  res_extended[, external_gene_name := NULL]; res_extended[, uniprot_id := NULL]; res_extended[, ensembl_gene_id := NULL]
  res_extended_id <- append_gene_symbols(res_extended, symbols(df.rfp), extend_id = TRUE)

  canonical_isoforms <- fread(file.path("/home/rstudio/livemount/Bio_data/references/homo_sapiens/", "canonical_isoforms.txt"), header = FALSE)[[1]]
  res_extended_id <- res_extended_id[id_original %in% canonical_isoforms,]

  res_deg <- DEG.analysis(df.rfp, counts = RFP_counts_canonical)
  res_deg[, external_gene_name := NULL]; res_deg[, uniprot_id := NULL]; res_deg[, ensembl_gene_id := NULL]
  res_deg <- append_gene_symbols(res_deg, symbols(df.rfp), extend_id = TRUE)
  gg_deg <- DEG.plot.static(res_deg, NULL, width = 8, height = 6,
                            plot.ext = ".pdf")

  res <- res_extended_id
  res <- res_deg
  DEG_plot(res_deg)
  DEG_plot(res_extended_id)

  DEG_plot <- function(dt, draw_non_regulated = TRUE,
                       xlim = ifelse(two_dimensions, "bidir.max", "auto"),
                       ylim = "bidir.max",
                       xlab = ifelse(two_dimensions, "RNA fold change (log2)", "Mean counts (log2)"),
                       ylab = ifelse(two_dimensions, "RFP fold change (log2)",  "Fold change (log2)"),
                       two_dimensions = ifelse("LFC" %in% colnames(dt), FALSE, TRUE),
                       color.values = c("No change" = "black", "Significant" = "red",
                                        "Buffering" = "purple", "mRNA abundance" = "darkgreen",
                                        "Expression" = "blue", "Forwarded" = "yellow",
                                        "Inverse" = "aquamarine", "Translation" = "orange4"),
                       format = "png") {
    DEG_plot_input_validation()
    if (nrow(dt) == 0) return(invisible(NULL))

    dt[, y_axis := if("rfp.lfc" %in% colnames(dt)) {rfp.lfc} else LFC]
    dt[, x_axis := if("rna.lfc" %in% colnames(dt)) {rna.lfc} else log2(meanCounts)]
    # dt[, size := c(1, 0.8)[1 + (Regulation == "No change")]]
    dt[, trace := paste0("id: ", id, "\nReg: ", Regulation)]

    contrasts <- unique(dt$contrast)
    x_titles <- xlab
    y_titles <- ylab
    zerolines <- two_dimensions

    shared_df <- SharedData$new(dt, key = ~id, group = "A")
    filter_box <- filter_select("id_search", "Search by ID:", shared_df, ~id)

    showlegends <- c(TRUE, rep(FALSE, length(contrasts) - 1))
    names(showlegends) <- contrasts
    if (two_dimensions) {
      values <- abs(c(dt$y_axis, dt$x_axis))
      max <- max(values[is.finite(values)], na.rm = TRUE) + 0.5
      max <- c(-max, max)
    } else {
      values <- dt$x_axis
      values <- values[is.finite(values)]
      max <- c(min(c(-0.5, values), na.rm = TRUE), max(values + 0.5, na.rm = TRUE))
    }
    print(max)

    gg2 <- lapply(contrasts, function(cont) {
      message(cont)
      shared_df_sub <- SharedData$new(dt[contrast == cont,], key = ~id, group = "A")

      plot_ly(shared_df_sub, x = ~x_axis, y = ~y_axis, color = ~ Regulation, colors = color.values,
              text = ~trace, hoverinfo = "text", legendgroup = ~Regulation,
              type = "scattergl", showlegend = showlegends[cont],
              mode = "markers", marker = list(opacity = 0.3), name = ~ Regulation) %>%
        layout(xaxis = list(zerolinecolor = "red", zeroline = zerolines),
               yaxis = list(title = list(text = y_titles,
                                         font = list(color = "black", weight = "bold", size = 16)),
                            zerolinecolor = "red"),
               shapes = if (two_dimensions) {
                 list(
                   list(type = "line", xref = "x", yref = "y",
                        x0 = max[1], x1 = max[2],
                        y0 = max[1], y1 = max[2],
                        line = list(color = "rgba(128, 128, 128, 0.5)", dash = "dash", width = 1)
                   )
                 )}
        )
    })
    combined <- subplot(gg2, nrows = 1, shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = TRUE,
                        which_layout = "merge") %>%
      layout(
        annotations = list(
          list(x = 0.5, y = -0.15, text = x_titles[1], showarrow = FALSE, xref = "paper", yref = "paper",
               font = list(color = "black", weight = "bold", size = 16))),
        margin = list(b = 55),
        xaxis = list(range = c(max[1], max[2]))
      ) %>% plotly::config(toImageButtonOptions = list(format = format), displaylogo = FALSE)

    widths <- cumsum(rep(1 / length(gg2), length(gg2)))
    center <- widths[1] / 2
    combined <- layout(combined, annotations = lapply(seq_along(contrasts), function(i) {
      list(x = widths[i] - center, y = 1.05, text = sub("Comparison: ", "", contrasts[i]), showarrow = FALSE, xref = 'paper', yref = 'paper', xanchor = 'center')
    }))
    combined_with_box <- bscols(list(filter_box, tags$br(), combined    # your plotly subplot or plot
    ), widths = 12
    )
    return(combined_with_box)
  }

  DEG_plot_input_validation <- function() {
    with(rlang::caller_env(), {
      # colnames(dt)[colnames(dt) == "rna.lfc"] <- "rna"
      # colnames(dt)[colnames(dt) == "rfp.lfc"] <- "rfp"
      # colnames(dt)[colnames(dt) == "te.lfc"] <- "te"

      stopifnot(is(dt, "data.table"))
      if (is.character(xlim)) stopifnot(xlim %in% c("bidir.max", "auto"))
      if (is.character(ylim)) stopifnot(ylim %in% c("bidir.max", "auto"))
      # Remove this line in next bioc
      if("variable" %in% colnames(dt))
        colnames(dt) <- gsub("variable", "contrast", colnames(dt))
      columns_must_exists <- if (two_dimensions) {
        c("contrast", "Regulation", "id", "rna.lfc", "rfp.lfc", "te.lfc")
      } else c("contrast", "Regulation", "id", "LFC", "meanCounts")

      dt[, Regulation := factor(Regulation, levels = names(color.values))]

      if (!(all(columns_must_exists %in% colnames(dt))))
        stop("dt must minimally contain the columns: ",
             paste(columns_must_exists, collapse = ", "))
      if (!draw_non_regulated) dt <- dt[Regulation != "No change",]
      if (nrow(dt) == 0) {
        warning("dt input had no valid rows to plot")
      }
    })
  }



  bs_inlined <- attachDependencies(combined_with_box,
                                   htmlDependencies(combined_with_box),
                                   append = TRUE)
  htmltools::save_html(htmltools::browsable(bs_inlined), file = file.path(QCfolder(df.rfp), "DTEG_analysis.html"))
  browseURL(file.path(QCfolder(df.rfp), "DTEG_analysis.html"))

  gorilla_output_dir <- file.path(QCfolder(df.rfp), "DTEG_analysis_subsets")
  DEG_gorilla(res_extended, gorilla_output_dir)


  #sapply(all_gorilla_ids$url, function(x) browseURL(x))
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Codon analysis (From WT rep 1 & HSR rep 1)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- Codon dwell times..")

  cds <- loadRegion(df.rfp, "cds", filterTranscripts(df.rfp, minThreeUTR = NULL))
  cds <- cds[names(cds) %in% canonical_isoforms]
  stopifnot(length(cds) > 0)

  pshifted_libs <- outputLibs(df.rfp, type = "cov", output.mode = "envirlist")
  codon_table <- codon_usage_exp(df.rfp, pshifted_libs, cds = cds)
  fwrite(codon_table, file.path(QCfolder(df.rfp), "codon_dwell_times_table.csv"))
  ggsave(file = file.path(QCfolder(df.rfp), "codon_dwell_time_with_start_stop.png"),
         codon_usage_plot(codon_table), width = 5, height = 10) # There is an increased dwell time on (R:CGC) of A-sites in both conditions
  ggsave(file = file.path(QCfolder(df.rfp), "codon_dwell_time_without_start_stop.png"),
         codon_usage_plot(codon_table, ignore_start_stop_codons = TRUE), width = 5, height = 10)
  # There is an increased dwell time on (R:CGG) of A-sites of HSP condition, why ?

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Peak detection (strong peaks in CDS)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- Ribo peaks")
  dt_peaks <- rbindlist(lapply(pshifted_libs, function(lib) findPeaksPerGene(cds, reads = lib, type = "max")), idcol = "file")
  fwrite(dt_peaks, file.path(QCfolder(df.rfp), "extreme_cds_peaks.csv"))
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # All ORF type predictions
  # Prediction using peridicity (Similar to RiboCode, ORFScore, minimum coverage, and comparison
  # to upstream and downstream window)
  # Will create 3 files in format (.rds), GRangesList of candidate ORFs, of predicted ORFs and a table
  # of all scores per ORF used for prediction
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Run on 2 first libraries
  message("-- Translon predictions")
  ORFik::detect_ribo_orfs(df.rfp, prediction_output_folder,
                          c("uORF", "uoORF", "annotated", "NTE", "NTT", "doORF", "dORF"),
                          startCodon = "ATG|CTG|TTG|GTG",
                          mrna = loadRegion(df.rfp, "mrna", names(cds)),
                          cds = cds) # Human also has a lot of ACG uORFs btw
  message("- Ribo-seq & RNA-seq analysis pipeline done")
}
