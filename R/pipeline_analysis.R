#' Analysis of Ribo-seq and RNA-seq
#'
#' Result will be saved in respective QC_STATS folders:
#' i.e. ribo-seq goes to QCfolder(df.rfp) and rna-seq goes to
#' QCfolder(df.rna). Update to your home directory if
#' you do not have livemount access, by doing
#' \code{df.rna@resultFolder <- "~/my_dir/"}.
#' @param df.rfp an ORFik experiment
#' @param df.rna an ORFik experiment
#' @param selected_isoforms a character vector with tx names,
#'  default: canonical_isoforms(df.rfp). For subset of analysis tha can
#'  subset transcripts, which to use.
#' @param exp_design character vector, default design(df.rfp, multi.factor = FALSE).
#' Which factors to run DTEG analysis on.
#' @return invisible(NULL), results saved to disc
#' @export
#' @examples
#' df.rfp <- read.experiment("vaccine_mod_bnt_111024-homo_sapiens_x_sars_cov2")
#' df.rna <- read.experiment("vaccine_mod_bnt_111024-homo_sapiens_x_sars_cov2_RNA-seq")
#' analysis_pipeline_ribo_rna(df.rfp, df.rna)
analysis_pipeline_ribo_rna <- function(df.rfp = read.experiment("Eleonora-homo_sapiens"),
                                       df.rna = read.experiment("Eleonora-homo_sapiens_RNA-seq"),
                                       output_dir = "auto",
                                       selected_isoforms = canonical_isoforms(df.rfp),
                                       exp_design = design(df.rfp, multi.factor = FALSE)) {
  message("- Ribo-seq & RNA-seq analysis_pipeline starting")
  if (identical(QCfolder(df.rfp), QCfolder(df.rna))) {
    stop("QCfolder of ribo-seq and rna-seq must be different!")
  }
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Run Annotation and alignment QC
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  output_dir <- QCfolder(df.rfp)
  qc_pipeline_generic(df.rna)
  analysis_pipeline_ribo(df.rfp, output_dir, selected_isoforms)
  analysis_pipeline_DTEG(df.rfp, df.rna, output_dir, selected_isoforms, exp_design)
  message("- Ribo-seq & RNA-seq analysis pipeline done")
}

qc_pipeline_generic <- function(df, output_dir = QCfolder(df)) {

  message("-- QC..")
  # General QC
  ORFikQC(df, out.dir = dirname(output_dir), create.ofst = FALSE,
          use_simplified_reads = FALSE)

  # PCA
  ggsave(filename = file.path(output_dir, paste0("PCAplot_", name(df), ".png")),
         ORFik:::pcaExperiment(df), height = 7, width = 7)

  is_ribo <- identical("RFP", libraryTypes(df))
  if (is_ribo) {
    RiboQC.plot(df, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
    message("-- Ribo heatmaps")
    remove.experiments(df)
    # before pshifting
    heatMapRegion(df, region = c("TIS", "TTS"), shifting = "5prime", type = "ofst",
                  outdir = file.path(output_dir, "heatmaps/pre-pshift/"))
    heatMapRegion(df, region = c("TIS", "TTS"), shifting = "5prime", type = "pshifted",
                  outdir = file.path(output_dir, "heatmaps/pshifted/"))
  }
  remove.experiments(df) # Remove loaded data (it is not pshifted)
}

analysis_pipeline_DTEG <- function(df.rfp = read.experiment("Eleonora-homo_sapiens"),
                                   df.rna = read.experiment("Eleonora-homo_sapiens_RNA-seq"),
                                   output_dir = QCfolder(df.rfp),
                                   selected_isoforms = canonical_isoforms(df.rfp),
                                   exp_design = design(df.rfp, multi.factor = FALSE)) {
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Differential translation analysis (e.g. condition: WT vs CO)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("-- Differential analysis (DESeq2 model)")
  symbols <- symbols(df.rfp)
  # We now run, and here get 11 unique DTEG genes
  RFP_counts <- countTable(df.rfp, region = "cds", type = "summarized")
  RNA_counts <- countTable(df.rna, region = "cds", type = "summarized")

  RFP_counts_canonical <- RFP_counts[rownames(RFP_counts) %in% selected_isoforms,]
  RNA_counts_canonical <- RNA_counts[rownames(RNA_counts) %in% selected_isoforms,]

  custom_for_sars_cov2 <- grep("x_sars_cov2$", name(df.rfp))
  if (custom_for_sars_cov2) {

    mat <- t(as.matrix(colSums(tail(assay(RFP_counts_canonical), 2))))
    rownames(mat) <- "vaccine"
    range <- rowRanges(tail(RFP_counts_canonical, 1))
    names(range) <- "vaccine"
    RFP_counts_canonical <- SummarizedExperiment(rbind(assay(RFP_counts_canonical), mat), colData = colData(RFP_counts_canonical),
                                                 rowRanges = c(rowRanges(RFP_counts_canonical), range))

    mat <- t(as.matrix(colSums(tail(assay(RNA_counts_canonical), 2))))
    rownames(mat) <- "vaccine"
    RNA_counts_canonical <- SummarizedExperiment(rbind(assay(RNA_counts_canonical), mat), colData = colData(RNA_counts_canonical),
                                                 rowRanges = c(rowRanges(RNA_counts_canonical), range))
    symbols <- rbindlist(list(symbols, data.table(ensembl_gene_id = "vaccine", external_gene_name = "vaccine", ensembl_tx_name = "vaccine")), fill = TRUE)
  }

  # If design is not full rank, paste columns together and empty the others.
  res <- DTEG.analysis(df.rfp, df.rna, output_dir, design = exp_design,
                       RFP_counts = RFP_counts_canonical, RNA_counts = RNA_counts_canonical)
  res_extended <- append_gene_symbols(res, symbols, extend_id = FALSE)
  fwrite(res_extended, file.path(output_dir, "DTEG_analysis.csv"))
  fwrite(res_extended[, .N, by = .(contrast, Regulation)],
         file.path(output_dir, "DTEG_analysis_summary.csv"))
  res_extended <- fread(file.path(output_dir, "DTEG_analysis.csv"))
  sink <- DTEG.plot(res_extended, output_dir, width = 8, height = 6,
                    plot.ext = c(".pdf", ".png",".jpg"), plot_to_console = FALSE)
  res_extended[, external_gene_name := NULL]; res_extended[, uniprot_id := NULL]; res_extended[, ensembl_gene_id := NULL]
  res_extended_id <- append_gene_symbols(res_extended, symbols, extend_id = TRUE)
  fwrite(res_extended, file.path(output_dir, "DTEG_analysis_with_gene_symbols.csv"))


  combined_with_box <- RiboCrypt:::DEG_plot(res_extended_id)
  bs_inlined <- htmltools::attachDependencies(combined_with_box,
                                              htmltools::htmlDependencies(combined_with_box),
                                              append = TRUE)
  htmltools::save_html(htmltools::browsable(bs_inlined), file = file.path(output_dir, "DTEG_analysis.html"))
  browseURL(file.path(output_dir, "DTEG_analysis.html"))

  gorilla_output_dir <- file.path(output_dir, "DTEG_analysis_subsets")
  res_extended_id_gorilla <- copy(res_extended_id)
  res_extended_id_gorilla[, total_genes_by_category := .N, by = .(contrast, Regulation)]
  res_extended_id_gorilla <- res_extended_id_gorilla[total_genes_by_category > 10]

  DEG_gorilla(res_extended_id_gorilla, gorilla_output_dir, organism(df.rfp))
  copy_results <- try(DEG_gorilla_copy_to_local(gorilla_output_dir))
  if (!is(copy_results, "try-error")) {
    dt <- DEG_gorilla_local_load_data(gorilla_output_dir)
    g <- DEG_gorilla_plot(dt)
    go_plot_file <- file.path(gorilla_output_dir, "GO_Enrichment_plot.jpg")
    ggsave(g, file = go_plot_file, width = 30, height = 20, dpi = 300)
  } else {
    saveRDS(FALSE, file.path(output_dir, "FAILED_TO_FETCH_GORILLA"))
    warning("Could not fetch htmls from gorilla, rerun later, will skip for now!")
  }
  #sapply(all_gorilla_ids$url, function(x) browseURL(x))
  return(invisible(NULL))
}

analysis_pipeline_ribo <- function(df.rfp, output_dir = QCfolder(df.rfp),
                                   selected_isoforms = canonical_isoforms(df.rfp)) {
  # qc_pipeline_generic(df.rfp, output_dir)

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Codon analysis (From WT rep 1 & HSR rep 1)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- Codon dwell times..")

  cds <- loadRegion(df.rfp, "cds", filterTranscripts(df.rfp, minThreeUTR = NULL))
  cds <- cds[names(cds) %in% selected_isoforms]
  stopifnot(length(cds) > 0)

  pshifted_libs <- outputLibs(df.rfp, type = "cov", output.mode = "envirlist")
  codon_table <- codon_usage_exp(df.rfp, pshifted_libs, cds = cds)
  fwrite(codon_table, file.path(output_dir, "codon_dwell_times_table.csv"))
  ggsave(file = file.path(output_dir, "codon_dwell_time_with_start_stop.png"),
         codon_usage_plot(codon_table), width = 5, height = 10) # There is an increased dwell time on (R:CGC) of A-sites in both conditions
  ggsave(file = file.path(output_dir, "codon_dwell_time_without_start_stop.png"),
         codon_usage_plot(codon_table, ignore_start_stop_codons = TRUE), width = 5, height = 10)
  # There is an increased dwell time on (R:CGG) of A-sites of HSP condition, why ?

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Peak detection (strong peaks in CDS)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- Ribo peaks")
  dt_peaks <- rbindlist(lapply(pshifted_libs, function(lib) findPeaksPerGene(cds, reads = lib, type = "max")), idcol = "file")
  fwrite(dt_peaks, file.path(output_dir, "extreme_cds_peaks.csv"))
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # All ORF type predictions
  # Prediction using peridicity (Similar to RiboCode, ORFScore, minimum coverage, and comparison
  # to upstream and downstream window)
  # Will create 3 files in format (.rds), GRangesList of candidate ORFs, of predicted ORFs and a table
  # of all scores per ORF used for prediction
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Run on 2 first libraries
  message("-- Translon predictions")
  prediction_output_folder <- file.path(output_dir, "translon_predictions")
  ORFik::detect_ribo_orfs(df.rfp, prediction_output_folder,
                          c("uORF", "uoORF", "annotated", "NTE", "NTT", "doORF", "dORF", "ncORF"),
                          startCodon = "ATG|CTG|TTG|GTG",
                          mrna = loadRegion(df.rfp, "mrna", names(cds)),
                          cds = cds) # Human also has a lot of ACG uORFs btw
}

analysis_pipeline_DEG <- function(df.rfp, output_dir = QCfolder(df.rfp),
                                  selected_isoforms = canonical_isoforms(df.rfp)) {
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Differential translation analysis (e.g. condition: WT vs CO)
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

  message("-- Differential analysis (DESeq2 model)")
  # The design is by default chosen by this factor: The condition column in this case
  design <- design(df.rfp, multi.factor = FALSE)
  symbols <- symbols(df.rfp)
  # We now run, and here get 11 unique DTEG genes
  RFP_counts <- countTable(df.rfp, region = "cds", type = "summarized")
  RFP_counts_canonical <- RFP_counts[rownames(RFP_counts) %in% selected_isoforms,]
  df.rfp$fraction <- NULL
  RFP_counts_canonical$fraction <- NULL

  pairs <- combn.pairs(unlist(df.rfp[, design]))
  pairs <- lapply(pairs, function(pair) {
    should_be_last <- c("WT", "ctrl", "control", "other")
    if (pair[1] %in% should_be_last) pair <- rev(pair)
    return(pair)
  })
  res <- DEG.analysis(df.rfp, counts = RFP_counts_canonical, pairs = pairs)
  res_extended <- append_gene_symbols(res, symbols, extend_id = FALSE)
  fwrite(res_extended, file.path(output_dir, "DTEG_analysis.csv"))
  res_extended <- fread(file.path(output_dir, "DTEG_analysis.csv"))
  sink <- DEG.plot.static(res_extended, output_dir, width = 8, height = 6,
                    plot.ext = c(".pdf", ".png",".jpg"))
  res_extended[, external_gene_name := NULL]; res_extended[, uniprot_id := NULL]; res_extended[, ensembl_gene_id := NULL]
  res_extended_id <- append_gene_symbols(res_extended, symbols, extend_id = TRUE)


  combined_with_box <- RiboCrypt:::DEG_plot(res_extended_id)
  bs_inlined <- htmltools::attachDependencies(combined_with_box,
                                              htmltools::htmlDependencies(combined_with_box),
                                              append = TRUE)
  htmltools::save_html(htmltools::browsable(bs_inlined), file = file.path(output_dir, "DEG_analysis.html"))
  browseURL(file.path(output_dir, "DEG_analysis.html"))

  gorilla_output_dir <- file.path(output_dir, "DEG_analysis_subsets")
  DEG_gorilla(res_extended_id, gorilla_output_dir, organism(df.rfp))
  DEG_gorilla_copy_to_local(gorilla_output_dir)

  #sapply(all_gorilla_ids$url, function(x) browseURL(x))
  return(invisible(NULL))
}


