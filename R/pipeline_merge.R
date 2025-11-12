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
    if (config$all_mappers)
      ORFik::mergeLibs(df, file.path(libFolder(df), "pshifted_merged"),
                       "lib", "pshifted", FALSE)

    if (config$split_unique_mappers) {
      uniqueMappers(df) <- TRUE
      ORFik::mergeLibs(df, file.path(libFolder(df), "pshifted_merged"),
                       "lib", "pshifted", FALSE)
    }

    set_flag(config, "merged_lib", name(df))
  }
}

#' Merge genome for mixed genome to align against
#'
#' All main species must have valid outputs.rds with gtf and fasta genome defined in
#' reference folder. Custom sequences (plasmids etc) requires fasta only, gtf can
#' be created for you.
#' TODO: Currently it fails if fasta or gtf files does not end with new line
#' @param species character pre-existing genomes in the ref_dir: use basename (e.g. homo_sapiens)
#' @param primary_organism character, scientific name (e.g. "Homo sapiens").
#' Default: \code{stringr::str_to_sentence(gsub("_", " ", species[1]))}
#' The name of the primary organism, so that biomart can search for symbols etc.
#' @param custom_sequences character, path to custom fasta sequences
#' @param custom_gtfs character, default make_gtf_from_sequences(custom_sequences),
#' path to custom fasta sequences gtfs if existing, NULL if custom_sequences is length 0.
#' @param ref_dir location of existing references, default: ORFik::config()["ref"]
#' @param contaminants_sequences path to STAR index to copy for contaminants, default "primary_only". Copies
#' the contaminants for species[1]
#' @param gene_symbols path to gene symbol fst to copy, default "primary_only". Copies
#' the fst for species[1]
#' @return invisible(NULL), results saved to disc
#' @export
#' @examples
#' # Make outputs.rds object
#' annotation <- c(genome = "~/livemount/Bio_data/references/influenza/PR8_JWY_RefGenome_v4.fasta",
#' gtf = "~/livemount/Bio_data/references/influenza/PR8_JWY_RefGenome_v4.gtf")
#' #saveRDS(annotation, file.path(dirname(annotation[1]), "outputs.rds"))
#' #indexFa(annotation["genome"])
#' # Custom plasmid
#' #custom_sequences <- "~/livemount/Bio_data/references/influenza/9SNLucRev_Seg9.fasta"
#' #make_mixed_genome(c("homo_sapiens", "influenza"), custom_sequences)
make_mixed_genome <- function(species, custom_sequences, custom_gtfs = make_gtf_from_sequences(custom_sequences),
                              primary_organism = stringr::str_to_sentence(gsub("_", " ", species[1])),
                              ref_dir = ORFik::config()["ref"],
                              contaminants_sequences = "primary_only",
                              gene_symbols = "primary_only") {
  stopifnot(length(species) > 0 & is.character(species))
  message("Making mix of ", length(species), " species")
  message("- ", paste(species, collapse = ", "))
  if (!is.null(custom_gtfs)) stopifnot(all(file.exists(custom_gtfs)))
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
  gtfs <- c(gtfs, custom_gtfs)
  genomes <- all_files[names(all_files) == "genome"]
  genomes <- c(genomes, custom_sequences)


  mix_species_name <- paste(species, collapse = "_x_")
  mix_species_folder <- file.path(ref_dir, mix_species_name)
  dir.create(mix_species_folder)

  message("Merging genomes (fasta files)")
  mixed_genome_path <- file.path(mix_species_folder, paste0(mix_species_name, ".fasta"))
  genome_cat <- paste("cat", paste(genomes, collapse = " "), ">", mixed_genome_path)
  system(genome_cat)
  indexFa(mixed_genome_path)

  message("Merging annotation (gtf files)")
  mixed_gtf_path <- file.path(mix_species_folder, paste0(mix_species_name, ".gtf"))
  gtf_cat <- paste("cat", paste(gtfs[1], collapse = " "), ">", mixed_gtf_path)
  system(gtf_cat)
  for (gtf in gtfs[-1]) {
    # Remove metadata header for remaining
    line_index <- first_non_metadata_line_annotation(gtf)
    tail_call <- paste0("tail -n +", line_index)
    gtf_cat <- paste(tail_call, paste(gtf, collapse = " "), ">>", mixed_gtf_path)
    system(gtf_cat)
  }


  output <- c(gtf = mixed_gtf_path, genome = mixed_genome_path)
  saveRDS(output, file.path(mix_species_folder, "outputs.rds"))

  annotation <- getGenomeAndAnnotation("temp", mix_species_folder, GTF = output["gtf"], genome = output["genome"])
  annotation <- c(annotation) # Contaminants from main species

  makeTxdbFromGenome(output["gtf"], output["genome"], organism = primary_organism, optimize = TRUE)
  if (gene_symbols == "primary_only") {
    gene_symbols <- file.path(primary_dir, "gene_symbol_tx_table.fst")
    if (file.exists(gene_symbols)) {
      file.copy(gene_symbols, file.path(mix_species_folder, "gene_symbol_tx_table.fst"), overwrite = TRUE)
      symbols <- fst::read_fst(gene_symbols, as.data.table = TRUE)
      txdb <- paste0(output["gtf"], ".db")
      txdb <- loadTxdb(txdb)
      tx <- loadRegion(txdb, "tx")
      tx <- tx[!(names(tx) %in% symbols$ensembl_tx_name)]
      symbols_append <- data.table(ensembl_gene_id = txNamesToGeneNames(names(tx), txdb),
                                   external_gene_name = as.character(NA),
                                   ensembl_tx_name = names(tx),
                                   uniprot_id = as.character(NA))
      symbols <- rbindlist(list(symbols, symbols_append))
      fst::write_fst(symbols, file.path(mix_species_folder, "gene_symbol_tx_table.fst"))

    }

    canonical_isoforms <- file.path(primary_dir, "canonical_isoforms.txt")
    if (file.exists(canonical_isoforms)) {
      file.copy(canonical_isoforms, file.path(mix_species_folder, "canonical_isoforms.txt"), overwrite = TRUE)
      isoforms <- fread(canonical_isoforms, header = FALSE)[[1]]
      isoforms_to_append <- symbols[!(ensembl_gene_id %in% symbols[ensembl_tx_name %in% isoforms,]$ensembl_gene_id),][!duplicated(ensembl_gene_id)]$ensembl_tx_name
      isoforms <- c(isoforms, isoforms_to_append)
      fwrite(data.table(isoforms), file.path(mix_species_folder, "canonical_isoforms.txt"), col.names = FALSE)
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

first_non_metadata_line_annotation <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  i <- 1
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) stop("There is no data in this gff/gtf!")# end of file (only comments)
    if (!grepl("^#", line)) return(i)      # first non-header line
    i <- i + 1
  }
}

#' Make gtf from fasta sequence files
#'
#' For each fasta, make a gtf, by default type is "cds_only", i.e. whole sequence is
#' cds.
#' @param fasta_files character, paths to fasta files
#' @param type character, default "cds_only", how ranges are defined TODO: implement
#' IRangesList input for non 1 cds start.
#' @param strand character, default "+", ranges are + stranded, else user defined of
#' length 1 or length of fasta files.
#' @return paths to gtfs, or NULL if length of fasta_files is 0.
#' @export
make_gtf_from_sequences <- function(fasta_files, type = "cds_only", strand = "+") {
  if (length(fasta_files) == 0) return(NULL)
  stopifnot(all(tools::file_ext(fasta_files) == "fasta"))
  stopifnot(type %in% c("cds_only"))
  stopifnot(length(strand) %in% c(1, length(fasta_files)))
  gtfs <- c()
  for (fasta_file in fasta_files) {
    message("Making gtf of fasta file: \n", fasta_file)
    indexFa(fasta_file)
    ranges <- scanFaIndex(fasta_file)
    strand(ranges) <- "+"
    seqnames = as.character(seqnames(ranges))
    gene_index <- seq_along(seqnames)
    mcols(ranges) <- DataFrame(source = "mNGSp", type = "gene",
                               score = as.numeric(NA), phase = as.integer(NA),
                               gene_id = paste0(seqnames, "_Gene_", gene_index),
                               transcript_id = paste0(seqnames, "_tx_", gene_index))
    ranges <- rep(ranges, each = 4)
    ranges$type <- c("gene", "transcript", "exon", "cds")
    gtf <- sub("fasta$", "gtf", fasta_file)
    rtracklayer::export(ranges, gtf, format = "gtf")
    gtfs <- c(gtfs, gtf)
  }
  return(gtfs)
}

