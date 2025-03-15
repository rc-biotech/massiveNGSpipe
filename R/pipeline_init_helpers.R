#' Init all studies to pipeline objects
#' @inheritParams curate_metadata
#' @param simple_progress_report logical, default TRUE. Display current progress.
#' @param only_complete_genomes logical, default FALSE. If TRUE, will only init the subset
#' with complete genome/annotation directories. Will fail if 0 species are done.
#' @param dbs character, prioritized order of genome databases to search, default:
#' c("ensembl", "refseq"). Options are c("ensembl", "refseq", "genbank") in any wanted
#' priority order.
#' @return a list of pipelines
#' @export
pipeline_init_all <- function(config, complete_metadata = config$complete_metadata,
                              simple_progress_report = TRUE, gene_symbols = TRUE,
                              only_complete_genomes = FALSE, dbs = c("ensembl", "refseq")) {
  if (!file.exists(complete_metadata)) stop("You have not create a successful metadata table yet!")
  final_list <- fread(complete_metadata)[KEEP == TRUE,]

  if (only_complete_genomes) {
    message("- Subsetting to only complete genomes")
    complete_genomes <- list.genomes(reference.folder = config$config["ref"])$name
    total_genomes <- length(unique(final_list$ScientificName))
    total_samples <- nrow(final_list)
    final_list <- final_list[ScientificName %in% gsub("_", " ", (stringr::str_to_title(complete_genomes))), ]

    message("Complete genomes ratio:", round((length(complete_genomes) / total_genomes)*100, 2), "%")
    message("Used samples ratio:", round((nrow(final_list) / total_samples)*100, 2), "%")
  }
  if (nrow(final_list) == 0) stop("complete metadata table has 0 rows, ",
                                  "did you forget to set the 'KEEP' column to TRUE in",
                                  " 'config$temp_metadata'?")
  if (is.null(final_list$ScientificName))
    stop("Complete metadata must contain organism information in column 'ScientificName'")

  # Fetch all organism references ++
  reference_list <- get_all_annotation_and_index(final_list, config,
                                                 gene_symbols, dbs)

  # For each accession run init
  accessions <- unique(final_list$study_accession)
  pipelines <- lapply(accessions, function(accession)
    tryCatch(
      pipeline_init(final_list[final_list$study_accession == accession,],
                    accession,  config, reference_list),
      error=function(e) {warning(e); return(NULL)}))
  names(pipelines) <- accessions

  crashed_studies <- pipelines %in% list(NULL)
  if (any(crashed_studies)) {
    pipelines <- pipelines[!crashed_studies]
  }


  if (simple_progress_report) progress_report(pipelines, config, FALSE)
  return(pipelines)
}

get_all_annotation_and_index <- function(final_list, config,
                                         gene_symbols = TRUE,
                                         dbs = c("ensembl", "refseq")) {
  message("- Fetching all organism required")
  unique_organisms <- unique(final_list$ScientificName)
  reference_list <- lapply(unique_organisms,
                           function(organism) {
                             get_annotation_and_index(organism,
                                                      file.path(config$config["ref"], reference_folder_name(organism)),
                                                      gene_symbols, dbs)
                             })
  names(reference_list) <- unique_organisms
  message("-------- Done --------")
  return(reference_list)
}

get_annotation_and_index <- function(organism,
                                     dir = reference_folder_name(organism,
                                                                 TRUE),
                                     gene_symbols = TRUE, dbs = c("ensembl", "refseq")) {
  message("---- ", organism)
  # Ensembl
  for (db in dbs) {
    annotation <- get_annotation(organism, dir, gene_symbols = gene_symbols,
                                 db = db)
    if (!is.null(annotation)) break
  }

  index <- NULL
  if (!is.null(annotation)) {
    get_file_size_gb <- function(file_path) {
      file_size_bytes <- file.info(file_path)$size
      file_size_gb <- file_size_bytes / (1024^3)  # Convert bytes to GB
      return(as.numeric(file_size_gb))  # Ensure it's numeric
    }

    # Try to estimate safe ram size
    max_ram <- get_file_size_gb(annotation["genome"])*15
    index <- ORFik::STAR.index(annotation, notify_load_existing = FALSE,
                               max.ram = max_ram)
  }
  return(list(annotation = annotation, index = index))
}

#' Get all required annotation for an organism
#'
#' Get: Fasta genome, transcriptome annotation (gff/gtf),
#' contaminants, ORFik speedup objects, gene symbols,
#' uniprot ids (to be added).
#' For local annotation, set genome and GTF arguments to full
#' path of respective objects.
#' @inheritParams ORFik::getGenomeAndAnnotation
#' @param output.dir path, default: reference_folder_name(organism, TRUE)
#' @return the annotation paths as vector
#' @export
get_annotation <- function(organism,
      output.dir = reference_folder_name(organism, TRUE),
      genome = TRUE, GTF = TRUE, db = "ensembl",
      gene_symbols = TRUE) {
  # Download if not existing, else just load
  db_result <- try(
    annotation <- getGenomeAndAnnotation(
      organism = organism, db = db,
      genome = genome, GTF = GTF,
      phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
      output.dir = output.dir, optimize = TRUE,
      pseudo_5UTRS_if_needed = 100, notify_load_existing = FALSE
    )
  )

  if (is(db_result, "try-error") & db == "refseq") {
    message("--- GTF failed, trying gff3 format")
    db_result <- try(
      annotation <- getGenomeAndAnnotation(
        organism = organism, db = db,
        genome = genome, GTF = GTF,
        phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
        output.dir = output.dir, optimize = TRUE,
        pseudo_5UTRS_if_needed = 100, notify_load_existing = FALSE,
        refseq_genbank_format = "gff3"
      )
    )
  }
  if (is(db_result, "try-error")) {
    message(db_result)
    warning("Failed download of organism by db = '", db, "' of: ", organism)
    return(NULL)
  }

  if (gene_symbols) {
    txdb <- loadTxdb(paste0(annotation["gtf"], ".db"))
    get_symbols(txdb, file.path(dirname(annotation["gtf"]), "gene_symbol_tx_table.fst"))
  }
  return(annotation)
}

#' @import fst biomartr
get_symbols <- function(txdb, path = file.path(dirname(ORFik:::getGtfPathFromTxdb(txdb)), "gene_symbol_tx_table.fst"),
                        org = organism(txdb), force = FALSE, verbose = FALSE,
                        uniprot_id = FALSE) {
  if (file.exists(path) & !force) {
    if (verbose) message("Loading pre-existing symbols from file")
    return(invisible(NULL))
  }
  not_supported_by_biomart <- file.path(dirname(path), "no_symbols_found.rds")
  if (file.exists(not_supported_by_biomart))
    return(invisible(NULL))

  organism <- org
  name <- paste0(tolower(substr(organism,
                                1, 1)), gsub(".* ", replacement = "", organism),
                 "_gene_ensembl")

  vertebrate_list <- biomartr::getDatasets("ENSEMBL_MART_ENSEMBL")$dataset
  non_vertebrate <- !(name %in% vertebrate_list)

  if (non_vertebrate) {
    ensembl_genomes_marts <- setDT(biomartr::getMarts())[grep("_mart", mart)]
    name <- paste0(tolower(substr(organism,
                                  1, 1)), gsub(".* ", replacement = "", organism), "_eg_gene")
    valid_mart <- NULL
    for (mart in ensembl_genomes_marts$mart) {
      test <- try(biomartr::getDatasets(mart, TRUE), silent = TRUE)
      is_in_biomart <- !is(test, "try-error") & (name %in% test$dataset)
      if (is_in_biomart) {
        valid_mart <- try(biomaRt::useEnsemblGenomes(mart, name))
        if (is(dt, "try-error")) message("Error was catched and ignored, report on github")
        break
      }
    }
    failed_fetch <- TRUE
    if (is_in_biomart) {
      dt <- try(geneToSymbol(txdb, include_tx_ids = TRUE, ensembl = valid_mart,
                         force = force, uniprot_id = uniprot_id))
      if (is(dt, "try-error")) {
        message("Error was catched and ignored, most likely no gene names exist,
                if not the case report on massiveNGSpipe github issue page")
      } else failed_fetch <- FALSE
    }
    if (failed_fetch) {
      warning("Species not support symbols on ensembl mart, ignoring symbols")
      saveRDS(FALSE, not_supported_by_biomart)
    }
  } else dt <- geneToSymbol(txdb, include_tx_ids = TRUE, force = force,
                            uniprot_id = uniprot_id)
  if (is(dt, "data.table")) fst::write_fst(dt, path)
  return(invisible(NULL))
}

update_path_per_sample <- function(df, libtype_df, rel_dir = "pshifted_merged", ext = ".ofst", libtype_to_merge) {
  if (length(libtype_df) == 0) stop("No libs of libtype: '", libtype_to_merge, "' in experiments given")
  if (length(libtype_df) != 1) stop("Only single libtype experiments supported for merging")
  if (libtype_df == "") stop("Libtype of experiment must be defined!")
  df@listData$filepath <- file.path(dirname(filepath(df, "default")),
                                    rel_dir, paste0(libtype_df, ext))
  df@listData$rep <- seq(nrow(df))
  return(df)
}

#' Create a pipeline for the specified study.
#'
#' @param study a data.table of output from ORFik::download.SRA.metadata
#' @param study_accession any accession accepted by
#' \code{ORFik::download.SRA.metadata}. It is needed to know
#' if you used GEO, PRJ, SRP etc as the search query.
#' @param config Configured directories for pipeline as a list
#' @return a list of pipeline objects, one for each study,
#' subsetted by organism per study.
pipeline_init <- function(study, study_accession, config, reference_list) {
  # For each organism in the study, create an ORFik experiment path config
  organisms <- list()
  for (organism in unique(study$ScientificName)) {
    assembly_name <- reference_folder_name(organism)
    experiment <- paste(study_accession, assembly_name, sep = "-")
    conf <- path_config(experiment, assembly_name, config)
    org_ref <- reference_list[[organism]]
    organisms[[organism]] <- list(
      conf = conf, annotation = org_ref$annotation,
      index = org_ref$index, experiment = experiment
    )
  }
  return(list(
    accession = study_accession, study = study, organisms = organisms
  ))
}

get_kingdoms <- function(file) {
  stopifnot(file.exists(file))
  d <- fread(file, header = TRUE, select = c("organism_name", "group"))[!duplicated(organism_name),]
  d[, subgroup := gsub(".*_", "", group)]
  d[, group := gsub("( |_).*", "", group)]
  return(d)
}

get_kingdom_for_species <- function(species, db = "ensembl",
                                    file = file.path(readRDS("~/.biomartr_cache_dir.rds"), "AssemblyFilesAllKingdoms_refseq.txt"),
                                    regex = FALSE) {
  if (db == "ensembl") {
    kingdoms_temp <- sapply(species, function(x) {
      res <- biomartr:::ensembl_assembly_hits(x)
      if (is.logical(res)) return(as.character(NA))
      return(res$division[1])
    })

    kingdoms <- gsub("Ensembl", "", kingdoms_temp)
    kingdoms[species %in% c("Aedes aegypti", "Saccharomyces cerevisiae")] <- c("Vertebrates", "Fungi")
    kingdoms[kingdoms %in% "Vertebrates"] <- "v/Animalia"
    return(kingdoms)
  }

  if (regex) {
    return(get_kingdoms(file)[grep(paste(species, collapse = "|"), organism_name)])
  } else return(get_kingdoms(file)[organism_name %in% species][chmatch(species, organism_name)]
                [is.na(organism_name), organism_name := species[which(is.na(organism_name))]])
}

path_config <- function(experiment, assembly_name, config, type = ifelse(config$preset == "Ribo-seq", "", config$preset)) {
  # Create ORFik experiment config manually, without separate folders
  # for different library strategies.
  conf <- config.exper(experiment, assembly_name, type, config[["config"]])
  # Fix bad naming from config.exper
  conf <- gsub("//", "/", conf); conf <- gsub("_$", "", conf)
  names(conf) <- gsub(" .*", "", names(conf))
  names(conf) <- gsub(" $", "", names(conf))
  sapply(conf[1:3], fs::dir_create)
  return(conf)
}
