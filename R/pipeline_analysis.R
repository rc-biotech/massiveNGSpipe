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

  canonical_isoforms <- canonical_isoforms(df.rfp)
  res_extended_id <- res_extended_id[id_original %in% canonical_isoforms,]

  combined_with_box <- DEG_plot(res_extended_id)
  bs_inlined <- attachDependencies(combined_with_box,
                                   htmlDependencies(combined_with_box),
                                   append = TRUE)
  htmltools::save_html(htmltools::browsable(bs_inlined), file = file.path(QCfolder(df.rfp), "DTEG_analysis.html"))
  browseURL(file.path(QCfolder(df.rfp), "DTEG_analysis.html"))
  # # Ribo-seq only analysis if no RNA-seq
  # res_deg <- DEG.analysis(df.rfp, counts = RFP_counts_canonical)
  # res_deg[, external_gene_name := NULL]; res_deg[, uniprot_id := NULL]; res_deg[, ensembl_gene_id := NULL]
  # res_deg <- append_gene_symbols(res_deg, symbols(df.rfp), extend_id = TRUE)
  # gg_deg <- DEG.plot.static(res_deg, NULL, width = 8, height = 6,
  #                           plot.ext = ".pdf")

  gorilla_output_dir <- file.path(QCfolder(df.rfp), "DTEG_analysis_subsets")
  DEG_gorilla(res_extended_id, gorilla_output_dir, organism(df.rfp))
  DEG_gorilla_copy_to_local(gorilla_output_dir)

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
