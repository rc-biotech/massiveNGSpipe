#' Init all studies to pipeline objects
#' @inheritParams curate_metadata
#' @param simple_progress_report logical, default TRUE. Display current progress.
#' @param only_complete_genomes logical, default FALSE. If TRUE, will only init the subset
#' with existing and complete genome/annotation directories. Will still run STAR index if not existing.
#' Mixed speces must have directory with all lowercase and "_" seperator. While in the sheet
#' each organism start with capital letter (e.g. Homo sapiens x Mus musculus), directory would be
#' homo_sapiens_x_mus_musculus. Will fail if 0 species are done.
#' @param only_complete_indices logical, default FALSE. If TRUE will only use complete index species,
#' ignored if only_complete_genomes is FALSE.
#' @param dbs character, prioritized order of genome databases to search, default:
#' c("ensembl", "refseq"). Options are c("ensembl", "refseq", "genbank") in any wanted
#' priority order.
#' @param verbose_load_annotation logical, default TRUE. Will list name of each organism loaded
#' @return a list of pipelines
#' @export
pipeline_init_all <- function(config, complete_metadata = config$complete_metadata,
                              progress_report = TRUE, gene_symbols = TRUE,
                              only_complete_genomes = FALSE, only_complete_indices = FALSE,
                              dbs = c("ensembl", "refseq"),
                              show_status_per_exp = TRUE, verbose_load_annotation = TRUE) {
  # Load all valid metadata
  final_list <- init_and_load_complete_metadata(config, complete_metadata,
                                                only_complete_genomes,
                                                only_complete_indices)
  # Fetch all organism references ++
  reference_list <- get_all_annotation_and_index(final_list, config,
                                                 gene_symbols, dbs,
                                                 verbose_load_annotation)
  # For each accession run pipeline_init
  pipelines <- finalize_pipeline_objects(final_list, config, reference_list)

  if (progress_report) progress_report(pipelines, config,
                                       show_status_per_exp = show_status_per_exp)
  return(pipelines)
}

init_and_load_complete_metadata <- function(config, complete_metadata = config$complete_metadata,
                                            only_complete_genomes = FALSE, only_complete_indices = FALSE) {
  if (!file.exists(complete_metadata)) stop("You have not create a successful metadata table yet!")
  final_list <- fread(complete_metadata)[KEEP == TRUE,]
  final_list <- blacklist_filter(final_list, blacklist_read(config, sync = FALSE))
  if (isTRUE(attr(final_list, "removed") > 0))
    message("- Blacklist: skipping ", attr(final_list, "removed"), " blacklisted samples")
  final_list <- platform_guard(final_list, config)

  if (only_complete_genomes) {
    message("- Subsetting to only complete genomes")
    genomes <- list.genomes(reference.folder = config$config["ref"])
    if (only_complete_indices) {
      message("-- Subsetting to only complete STAR indices")
      genomes <- genomes[STAR_index == TRUE,]
    }
    complete_genomes <- genomes$name
    total_genomes <- length(unique(final_list$ScientificName))
    total_samples <- nrow(final_list)
    complete_genomes_formated <- gsub("_", " ", (stringr::str_to_title(complete_genomes)))
    mixed_genomes <- grep("_x_", complete_genomes)
    complete_genomes_formated[mixed_genomes] <- complete_genomes[mixed_genomes]
    alternative_mixed_format <- unlist(lapply(strsplit(complete_genomes[mixed_genomes], "_x_"), function(x) {
      x <- stringr::str_to_title(x)
      gsub("_", " ", paste(x, collapse = " x "))
    }))
    complete_genomes_formated <- c(complete_genomes_formated, alternative_mixed_format)
    # match through the folder name (same mapping as reference_folder_name()), so strain keys such as
    # "Zymomonas mobilis ZM4" -> zymomonas_mobilis_zm4 are found; title-casing the folder name could not round-trip
    final_list <- final_list[ScientificName %in% complete_genomes_formated |
                             reference_folder_name(ScientificName) %in% complete_genomes,]

    message("Complete genomes ratio:", round((length(complete_genomes) / total_genomes)*100, 2), "%")
    message("Used samples ratio:", round((nrow(final_list) / total_samples)*100, 2), "%")
  }
  if (nrow(final_list) == 0) stop("complete metadata table has 0 rows, ",
                                  "did you forget to set the 'KEEP' column to TRUE in",
                                  " 'config$temp_metadata'?")
  if (is.null(final_list$ScientificName))
    stop("Complete metadata must contain organism information in column 'ScientificName'")
  return(final_list)
}

get_all_annotation_and_index <- function(final_list, config, gene_symbols = TRUE,
                                         dbs = c("ensembl", "refseq"),
                                         verbose = TRUE) {
  message("- Fetching all organism required / Loading valid existing")
  unique_organisms <- sort(unique(final_list$ScientificName))
  reference_list <- lapply(unique_organisms,
                           function(organism) {
                             get_annotation_and_index(organism,
                                                      file.path(config$config["ref"], reference_folder_name(organism)),
                                                      gene_symbols, dbs, verbose = verbose)
                             })
  names(reference_list) <- unique_organisms
  message("-------- Done --------")
  return(reference_list)
}

get_annotation_and_index <- function(organism,
                                     dir = reference_folder_name(organism,
                                                                 TRUE),
                                     gene_symbols = TRUE, dbs = c("ensembl", "refseq"),
                                     verbose = TRUE) {
  if (verbose) message("---- ", organism)
  # A curated assembly with a known database (reference_db) is downloaded from that database only
  pinned <- reference_db(organism)
  if (!is.na(pinned)) dbs <- pinned
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
    max_ram <- get_file_size_gb(annotation["genome"])*20
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
#' @param assembly assembly to download, default \code{reference_assembly(organism)}: the curated reference strain for
#'  organisms where Ensembl returns several strains (bacteria), else the organism name (first assembly, ORFik default)
#' @export
get_annotation <- function(organism,
      output.dir = reference_folder_name(organism, TRUE),
      genome = TRUE, GTF = TRUE, db = "ensembl",
      gene_symbols = TRUE, assembly = reference_assembly(organism)) {
  # Download if not existing, else just load
  db_result <- try(
    annotation <- getGenomeAndAnnotation(
      organism = organism, db = db,
      genome = genome, GTF = GTF, assembly = assembly,
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
        genome = genome, GTF = GTF, assembly = assembly,
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
  dt <- NULL
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
      dt <- try(symbols_from_gtf(ORFik:::getGtfPathFromTxdb(txdb)), silent = TRUE)
      if (is(dt, "try-error")) {
        warning("Species had no supported symbols on ensembl mart or in the annotation file,
                ignoring symbols!")
        saveRDS(FALSE, not_supported_by_biomart)
      }
    }
  } else dt <- geneToSymbol(txdb, include_tx_ids = TRUE, force = force,
                            uniprot_id = uniprot_id)
  if (is(dt, "data.table")) fst::write_fst(dt, path)
  return(invisible(NULL))
}

symbols_from_gtf <- function(gtf_path, gtf = import(gtf_path)) {
  stopifnot(file.exists(gtf_path))
  stopifnot(is(gtf, "GRanges"))
  dt <- as.data.table(gtf)
  stopifnot(c("gene_id", "transcript_id", "gene") %in% colnames(dt))

  dt <- dt[type != "gene" & !duplicated(transcript_id)][, .(gene_id, transcript_id, gene, protein_id)]
  colnames(dt) <- c("ensembl_gene_id", "external_gene_name", "ensembl_tx_name", "uniprot_id")
}

update_path_per_sample <- function(df, libtype_df, rel_dir = "pshifted_merged",
                                   ext = ".ofst", libtype_to_merge, check_exist = TRUE) {
  if (length(libtype_df) == 0) stop("No libs of libtype: '", libtype_to_merge, "' in experiments given")
  if (length(libtype_df) != 1) stop("Only single libtype experiments supported for merging")
  if (libtype_df == "") stop("Libtype of experiment must be defined!")


  df@listData$filepath <- file.path(libFolder(df, mode = "all"),
                                    rel_dir, paste0(libtype_df, ext))
  df@listData$rep <- seq(nrow(df))
  if (check_exist) {
    file_paths <- df$filepath
    if (!all(file.exists(file_paths))) {
      stop("Missing ", rel_dir, " tracks for exps:",
           paste(file_paths[!file.exists(file_paths)], collapse = ", "))
    }
  }
  return(df)
}

#' For each accession run pipeline_init
#' @noRd
finalize_pipeline_objects <- function(final_list, config, reference_list) {
  accessions <- unique(final_list$study_accession)
  pipelines <- lapply(accessions, function(accession)
    tryCatch(
      pipeline_init(final_list[final_list$study_accession == accession,],
                    accession,  config, reference_list),
      error=function(e) {warning(e); return(NULL)}))
  names(pipelines) <- accessions

  crashed_studies <- lengths(pipelines) == 0
  if (any(crashed_studies)) {
    warning("Some studies aborted on init, will subset to working studies")
    pipelines <- pipelines[!crashed_studies]
  }
  return(pipelines)
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
  conf_not_exist <- conf[1:3]
  conf_not_exist <- conf_not_exist[!file.exists(conf_not_exist)]
  sapply(conf_not_exist, fs::dir_create)
  return(conf)
}

#' Build a BPPARAM for one pipeline step from config
#'
#' @param config the mNGSp config object
#' @param step character, must be a name in \code{config$threads}
#' @param workers numeric, default \code{config$threads[[step]]}. Pass
#' an explicit override to use a dynamically-computed worker count
#' instead of the static config value for this one call (e.g.
#' \code{\link{memory_safe_worker_count}()} for a memory-heavy step) --
#' every other behavior (thread_type, logging) still comes from
#' \code{config} as usual.
#' @return a BPPARAM object
#' @noRd
bpparam_from_config <- function(config, step, workers = config$threads[[step]]) {
  stopifnot(is(config$threads, "list"))
  stopifnot(step %in% names(config$threads))
  stopifnot(length(workers) == 1 && is.numeric(workers))
  stopifnot(is(config$thread_type, "function"))
  parallel_conf <- config$parallel_conf
  stopifnot(is.list(parallel_conf) && !is.null(names(parallel_conf)))
  stopifnot(all(names(parallel_conf) %in% c("log", "logdir", "jobname", "stop.on.error")))

  # Warn (never block -- this can be a deliberate choice) when this
  # step's own worker count, multiplied by threads$main's own
  # stage-group worker count, is disproportionate to detected cores.
  # Generalizes the lesson behind threads_blas_cap (pshifted/
  # valid_pshift/pcounts specifically, see its own roxygen,
  # R/pipeline_config.R) to ANY current or future step: each forked
  # worker can independently spin up its own full-width BLAS/OpenMP
  # thread pool regardless of how many sibling workers exist, and the
  # aggregate across threads$main x this step's own workers can
  # exhaust a container's cgroup pids.max outright -- confirmed live,
  # 2026-10-08 (pthread_create() -> EAGAIN "Resource temporarily
  # unavailable", reproduced in a synthetic repro at matching scale).
  if (step != "main" && !is.null(config$threads$main) &&
      workers > 1 && config$threads$main > 1) {
    total_estimate <- config$threads$main * workers
    cores <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    if (!is.na(cores) && total_estimate > cores * 4) {
      warning("config$threads$main (", config$threads$main, ") x threads$", step, " (", workers,
              ") = ", total_estimate, " potential concurrent workers, well above detected cores (",
              cores, "). Each one may ALSO spin up its own BLAS/OpenMP thread pool independently, ",
              "which can exhaust a container's cgroup pids.max outright (confirmed live, 2026-10-08) ",
              "-- consider lowering one of these, or capping via a threads_blas_cap-style bound.",
              call. = FALSE)
    }
  }

  if (step != "main") {
    if (!is.null(parallel_conf$logdir) & !is.na(parallel_conf$logdir)) {
      parallel_conf$logdir <- file.path(parallel_conf$logdir, step)
      parallel_conf$jobname <- step
    }
  }

  if (!is.null(parallel_conf$logdir) & !is.na(parallel_conf$logdir) &
      !dir.exists(parallel_conf$logdir)) {
    dir.create(parallel_conf$logdir, FALSE, TRUE)
  }

  return(bpparam_with_options(config, workers, func = config$thread_type,
                              parallel_conf = parallel_conf))
}

bpparam_with_options <- function(config, workers = config$threads$default,
                                 func = config$thread_type,
                                 parallel_conf = config$parallel_conf) {
  is_serial <- workers == 1 | identical(func, BiocParallel::SerialParam)
  if (is_serial) return(SerialParam(log = parallel_conf$log,
                                    logdir = parallel_conf$logdir,
                                    jobname = parallel_conf$jobname,
                                    stop.on.error = parallel_conf$stop.on.error,
                                    progressbar = FALSE))

  func(workers, log = parallel_conf$log,
       logdir = parallel_conf$logdir,
       jobname = parallel_conf$jobname,
       stop.on.error = parallel_conf$stop.on.error)
}

#' Curated reference assembly per organism
#'
#' Ensembl Bacteria returns many strains per species and ORFik's default (\code{assembly = organism}) takes the first,
#' which gave wrong-strain references (B. subtilis subsp. spizizenii TU-B-10 for strain 168 data, 2026-10-14). Organisms
#' listed here get the standard laboratory reference strain; others keep the default. Strain-specific studies (e.g.
#' Zymomonas mobilis ZM4) get their own organism key and reference folder instead.
#' @param organism scientific name(s)
#' @return assembly accession(s), or the organism name when not curated
#' @export
reference_assembly <- function(organism) {
  curated <- c("Escherichia coli" = "GCA_000005845",        # K-12 MG1655
               "Bacillus subtilis" = "GCA_000009045",       # subsp. subtilis str. 168
               # fungal batch 2026-10-16 (see reference_db() for the database of each)
               "Saccharomyces mikatae" = "GCF_947241705.1",          # IFO 1815
               "Saccharomyces kudriavzevii" = "GCF_947243775.1",     # newer than Ensembl IFO1802 v1.0
               "Saccharomyces arboricola" = "GCF_000292725.1",       # H-6 (RefSeq, not flagged reference)
               "Saccharomyces jurei" = "GCA_900290405.1",            # GenBank only
               "Saccharomyces uvarum" = "GCA_947243795.1",           # GenBank only
               "Cryptococcus neoformans var. grubii" = "GCF_000149245.1",  # H99 (CNA3); not JEC21
               "Pyricularia oryzae" = "GCF_000002495.2",             # 70-15 (MG8), classic reference
               "Fusarium graminearum" = "GCF_000240135.3",           # PH-1
               "Penicillium oxalicum" = "GCF_001723175.1",
               "Naumovozyma castellii" = "GCF_000237345.1")       # CBS 4309
  out <- unname(curated[organism]); out[is.na(out)] <- organism[is.na(out)]
  out
}

#' Database to download the curated assembly from (organisms without an entry use the normal db order)
#' @param organism scientific name(s)
#' @return character, "refseq" / "genbank" or NA
#' @noRd
reference_db <- function(organism) {
  acc <- reference_assembly(organism)
  ifelse(grepl("^GCF_", acc), "refseq", ifelse(grepl("^GCA_", acc) & !organism %in% c("Escherichia coli", "Bacillus subtilis"),
                                              "genbank", NA_character_))
}

#' Exclude runs whose sequencing platform the short-read pipeline cannot process
#'
#' Long-read runs (Oxford Nanopore direct RNA, PacBio) and colour-space runs (ABI SOLiD: FASTQ in 0/1/2/3 colour calls,
#' not bases) are valid metadata (KEEP stays TRUE in the sheets) but must not enter the short-read base-space
#' trim/collapse/STAR presets. They are skipped at init with reason \code{requires_longread_processing} or
#' \code{requires_colorspace_decoding} and listed in \code{<project>/excluded_runs_platform.csv} for a future preset
#' (codex-rna / user 2026-10-14: EndoC-betaH1 PRJNA703954 direct RNA; SOLiD PRJNA603782 colour-space FASTQ confirmed).
#' @param final_list metadata rows (needs Platform; returned unchanged without it)
#' @param config pipeline config
#' @param longread,colorspace platforms excluded
#' @return final_list without the excluded runs
#' @noRd
platform_guard <- function(final_list, config, longread = c("OXFORD_NANOPORE", "PACBIO_SMRT"), colorspace = "ABI_SOLID") {
  if (!"Platform" %in% names(final_list)) return(final_list)
  pl <- toupper(final_list$Platform)
  ex <- pl %in% c(longread, colorspace)
  if (any(ex)) {
    out <- final_list[ex, intersect(c("Run", "BioProject", "ScientificName", "LIBRARYTYPE", "Platform", "Model"), names(final_list)), with = FALSE]
    out[, reason := ifelse(pl[ex] %in% colorspace, "requires_colorspace_decoding", "requires_longread_processing")]
    if (!is.null(config$project)) fwrite(out, file.path(config$project, "excluded_runs_platform.csv"))
    message("- Platform guard: skipping ", sum(ex), " long-read/colour-space runs (", paste(unique(out$Platform), collapse = ", "),
            "); see excluded_runs_platform.csv")
  }
  final_list[!ex]
}
