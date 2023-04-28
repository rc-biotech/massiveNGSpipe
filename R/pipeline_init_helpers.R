path_config <- function(experiment, assembly_name, config) {
  # Create ORFik experiment config manually, without separate folders
  # for different library strategies.
  conf <- config.exper(experiment, assembly_name, "", config[["config"]])
  # Fix bad naming from config.exper
  conf <- gsub("//", "/", conf); conf <- gsub("_$", "", conf)
  names(conf) <- gsub(" $", "", names(conf))
  sapply(conf[1:3], fs::dir_create)
  return(conf)
}

#' Get all required annotation for an organism
#'
#' Get: Fasta genome, transcriptome annotation (gff/gtf),
#' contaminants, ORFik speedup objects, gene symbols,
#' uniprot ids (to be added).
#' For local annotation, set genome and GTF arguments to full
#' path of respective objects.
#' @inheritParams ORFik::getGenomeAndAnnotation
#' @param output.dir path, file.path(ORFik::config()["ref"], organism)
#' @return the annotation paths as vector
#' @export
get_annotation <- function(organism,
      output.dir = file.path(ORFik::config()["ref"], organism),
      genome = TRUE, GTF = TRUE, db = "ensembl",
      gene_symbols = TRUE) {
  # Download if not existing, else just load
  ensembl_db <- try(
    annotation <- getGenomeAndAnnotation(
      organism = organism,
      genome = genome, GTF = GTF,
      phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
      output.dir = output.dir,
      assembly_type = "primary_assembly", optimize = TRUE,
      pseudo_5UTRS_if_needed = 100, notify_load_existing = FALSE,
      gene_symbols = TRUE
    )
  )
  if (is(ensembl_db, "try-error")) {
    message(ensembl_db)
    stop("Failed download of organism: ", organism,
         "if organism is not in ensembl, try db = 'refseq'")
  }
  return(annotation)
}
