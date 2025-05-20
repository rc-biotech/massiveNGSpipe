analysis_pipeline_ribo_rna <- function(df.rfp = read.experiment("Eleonora-homo_sapiens"),
                                       df.rna = read.experiment("Eleonora-homo_sapiens_RNA-seq"),
                                       canonical_isoforms = fread(file.path(refFolder(df.rfp), "canonical_isoforms.txt"), header = FALSE)[[1]]) {
  message("- Ribo-seq & RNA-seq analysis_pipeline starting")

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Convert files and run Annotation vs alignment QC
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  message("-- QC..")
  # General QC
  ORFikQC(df.rfp, create.ofst = FALSE)
  ORFikQC(df.rna, create.ofst = FALSE)
  ggsave(filename = file.path(dirname(df.rfp$filepath[1]), "QC_STATS",
                              paste0("PCAplot_", df.rfp@experiment, ".png")), ORFik:::pcaExperiment(df.rfp), height = 7, width = 7)
  ggsave(filename = file.path(dirname(df.rna$filepath[1]), "QC_STATS",
                              paste0("PCAplot_", df.rna@experiment, ".png")), ORFik:::pcaExperiment(df.rna), height = 7, width = 7)
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
  res <- DTEG.analysis(df.rfp, df.rna,
                       RFP_counts = countTable(df.rfp, region = "cds", type = "summarized",
                                               count.folder = "pshifted"))
  fwrite(res, file.path(QCfolder(df.rfp), "DTEG_analysis.csv"))
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
