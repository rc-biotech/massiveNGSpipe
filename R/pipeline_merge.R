# Pipeline, Merge libraries: per study -> all merged per organism
pipe_merge_study <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "merged_lib", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipeline_merge_study_single(df_list, config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, merge_exp, study: ", experiments[1])
  }
}

pipeline_merge_study_single <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "merged_lib", name(df))) next
    ORFik::mergeLibs(df, file.path(libFolder(df), "pshifted_merged"),
                     "lib", "pshifted", FALSE)
    set_flag(config, "merged_lib", name(df))
  }
}

#' Merge genome for mixed genome to align against
#'
#' @param species character pre-existing genomes in the ref_dir: use basename (e.g. homo_sapiens)
#' @param primary_organism character, scientific name (e.g. "Homo sapiens").
#' The name of the primary organism, so that biomart can search for symbols etc.
#' @param custom_sequences character, path to custom fasta sequences
#' @param custom_gtfs character, default NULL, path to custom fasta sequences if existing
#' @param ref_dir location of existing references, default: ORFik::config()["ref"]
#' @param contaminants_sequences path to STAR index to copy for contaminants, default "primary_only". Copies
#' the contaminants for species[1]
#' @param gene_symbols path to gene symbol fst to copy, default "primary_only". Copies
#' the fst for species[1]
#' @return invisible(NULL), results saved to disc
#' @export
make_mixed_genome <- function(species, custom_sequences, custom_gtfs = NULL,
                              primary_organism = stringr::str_to_title(gsub("_", " ", species[1])),
                              ref_dir = ORFik::config()["ref"],
                              contaminants_sequences = "primary_only",
                              gene_symbols = "primary_only") {
  stopifnot(length(species) > 0 & is.character(species))
  message("Making mix of ", length(species), " species")
  message("- ", paste(species, collapse = ", "))
  if (!is.null(custom_gtfs)) stopifnot(all(file.exists(mixed_gtf_path)))
  if (!is.null(custom_sequences)) {
    stopifnot(all(file.exists(custom_sequences)))
    message("With ", length(custom_sequences), " custom sequences")
    message("- ", paste(custom_sequences, collapse = ", "))
  }
  species_dirs <- file.path(ref_dir, species)
  primary_dir <- species_dirs[1]
  stopifnot(all(dir.exists(species_dirs)))
  genome_config <- lapply(file.path(species_dirs, "outputs.rds") , function(x) readRDS(x))



  contaminants_index <- NULL
  if (contaminants_sequences == "primary_only") {
    contaminants_index <- file.path(primary_dir, "STAR_index", "contaminants_genomeDir")
    if (!file.exists(contaminants_index)) {
      stop("STAR index contamint does not exists")
    }
    # annotation_primary <- getGenomeAndAnnotation("temp", primary_dir, GTF = output["gtf"], genome = output["genome"])
    # contaminants <- annotation_primary["contaminants"]
    message("Contamination used: \n", contaminants_index)
  }


  all_files <- unlist(genome_config)
  gtfs <- all_files[names(all_files) == "gtf"]
  genomes <- all_files[names(all_files) == "genome"]
  genomes <- c(genomes, custom_sequences)


  mix_species_name <- paste(species, collapse = "_x_")
  mix_species_folder <- file.path(ref_dir, mix_species_name)
  dir.create(mix_species_folder)
  mixed_gtf_path <- file.path(mix_species_folder, paste0(mix_species_name, ".gtf"))
  mixed_genome_path <- file.path(mix_species_folder, paste0(mix_species_name, ".fasta"))
  genome_cat <- paste("cat", paste(genomes, collapse = " "), ">", mixed_genome_path)
  system(genome_cat)
  indexFa(mixed_genome_path)


  gtf_cat <- paste("cat", paste(gtfs[1], collapse = " "), ">", mixed_gtf_path)
  system(gtf_cat)

  gtf_cat <- paste("tail -n +5", paste(gtfs[-1], collapse = " "), ">>", mixed_gtf_path)
  system(gtf_cat)

  output <- c(gtf = mixed_gtf_path, genome = mixed_genome_path)
  saveRDS(output, file.path(mix_species_folder, "outputs.rds"))

  annotation <- getGenomeAndAnnotation("temp", mix_species_folder, GTF = output["gtf"], genome = output["genome"])
  annotation <- c(annotation) # Contaminants from main species

  makeTxdbFromGenome(output["gtf"], output["genome"], organism = primary_organism, optimize = TRUE)
  if (gene_symbols == "primary_only") {
    gene_symbols <- file.path(primary_dir, "gene_symbol_tx_table.fst")
    if (file.exists(gene_symbols)) {
      file.copy(gene_symbols, file.path(mix_species_folder, "gene_symbol_tx_table.fst"))
    }
  }
  index <- STAR.index(annotation, remake = TRUE) # Only primary, contam is copied
  contaminants_index_files <- list.files(contaminants_index, recursive = FALSE, include.dirs = FALSE, full.names = TRUE)
  new_genome_dir <-  file.path(mix_species_folder, "STAR_index", "contaminants_genomeDir")
  dir.create(new_genome_dir)
  file.copy(contaminants_index_files, new_genome_dir)
  message("Done (files saved to):")
  message("- ", mix_species_folder)
  return(invisible(NULL))
}

