#' Curate metadata for massive_NGS_pipe
#'
#' Runs an interactive while loop until your metadata is curated and have this property:
#' Per accession, your study must have a unique set of rows per organism.
#' Such that ORFik experiment can understand what is what:
#'  (what is Wild type vs experimental), which are replicates etc.
#'  The uniqueness is defined from column LIBRARYTYPE to column INHIBITOR
#' @param accessions character vector, all candidate accession numbers, allowed types:\cr
#' - Biorpoject ID (Only numbers)\cr
#' - Bioproject accession (PRJ)\cr
#' - SRA study (SRA)\cr
#' - ENA study (ERA)\cr
#' - GEO study (GSE)\cr
#' - NULL (not "NULL", but NULL), don't add new accession, only look if you
#' current google sheet or local csv is ready.
#' @inheritParams pipeline_config
#' @inheritParams run_pipeline
#' @param organisms character vector, default "all" (Use all organisms found).
#' Else binomial latin name with capital letter for genus: "Homo sapiens" etc.
#' @param complete_metadata path, default file.path(project_dir, "RFP_FINAL_LIST.csv").
#' The list of final candidates that are checked and have unique rows per bioproject
#' @param LibraryStrategy character vector, default:
#' c("RNA-Seq", "miRNA-Seq", "OTHER"). Which library Strategies to keep.
#' Default is the total set of normal RNA sequencing. (RNA-seq, miRNA-seq,
#' Ribo-seq, CHIP-seq, CAGE etc.) Change if you only want a subset, or some
#' rare type.
#' @param libtypes character vector, default "RFP". Which NGS library types to
#' keep from all LIBRARYLAYOUT values. These names are made by mNGSp
#' automatically. If you made "disome" as a LIBRARYLAYOUT category, then
#' c("RFP", "disome") will keep both, but exclude any other.
#' @param LibraryLayouts character vector, default c("SINGLE", "PAIRED"),
#' either or both of: c("SINGLE", "PAIRED")
#' @param Platforms character vector, default: "ILLUMINA". The sequencer technologies allowed.
#' @param step_mode logical, default FALSE. If TRUE, actives browser() and lets you go through
#' each step in debug mode for full control.
#' @param fix_loop logical, default TRUE. If TRUE, will create a loop in
#' the end that runs until your metadata is valid. If you
#' only want to updata csv / google docs with new data and not try to complete
#' now, set this to FALSE.
#' @param open_editor logical default TRUE. If google URL is defined, open the
#' sheet. If not, open the DataEditR bord.
#' @param only_curated logical FALSE, if TRUE. Only validate (or fail) based on
#' inserted accessions (not the full list). Ignored if accessions is NULL.
#' @param update_google_sheet logical TRUE, Ignored if google_url is NULL.
#' If FALSE, will not update sheet after check. This
#' @return logical, TRUE if you had unique rows per accession per organism
#' @import googlesheets4 data.table ORFik BiocParallel fs stringr DataEditR
#' @export
curate_metadata <- function(accessions, config, organisms = "all",
                            google_url = config$google_url,
                            complete_metadata = config$complete_metadata,
                            LibraryLayouts = c("SINGLE", "PAIRED"),
                            LibraryStrategy = c("Ribo-seq","RNA-Seq", "miRNA-Seq", "OTHER"),
                            libtypes = "RFP",
                            Platforms = "ILLUMINA",
                            step_mode = FALSE, open_editor = interactive(),
                            fix_loop = TRUE, only_curated = FALSE,
                            update_google_sheet = TRUE) {
  if (!interactive() & step_mode)
    stop("In non interactive mode you can not run step_mode = TRUE!")
  stopifnot(length(fix_loop) == 1); stopifnot(is(fix_loop, "logical"))
  add_new_accessions <- !is.null(accessions)
  if (!update_google_sheet & add_new_accessions)
    message("update_google_sheet is FALSE and add_new_accessions is TRUE, is this intentional?")
  if (step_mode) browser()

  if (add_new_accessions) {
    add_new_data(accessions, config, organisms,
                 google_url, complete_metadata,
                 LibraryStrategy,
                 LibraryLayouts, Platforms,
                 open_editor)
  }

  # Step6: Now, Check if it is valid (if not repeat step with new csv)
  while(fix_loop) {
    readline(prompt = "You think metadata is ready?\n Press enter when ready: ")
    if (!is.null(google_url)) {
      message("- Reading google sheet")
      sheet <- read_sheet_safe(google_url)
    } else sheet <- fread(config$temp_metadata)
    finished <- pipeline_validate_metadata(sheet, config, accessions = accessions,
                                           only_curated = only_curated,
                                           libtypes = libtypes)
    if (!is.null(google_url) & update_google_sheet) {
      message("Uploading updated version to google sheet:")
      write_sheet(read.csv(config$temp_metadata),
                  ss = google_url,
                  sheet = 1)
    } else message("Updated your local csv metadata file")
    if (finished) break
  }
  return(fix_loop)
}

#' Download all metadata for all accessions
#' @inheritParams curate_metadata
#' @param force logical, default FALSE. Force redownload of all, if something failed
#' @param max_attempts numeric, default 7. How many attempts to download a specific study beforing
#' failing? This call flood SRA with calls, so set it higher if you use a lot
#' of samples. Will use a random increasing wait time between study attempts.
#' @param max_attempts_loop numeric, default 7. How many attempts for full download loop beforing
#' failing? This call flood SRA with calls, so set it higher if you use a lot
#' of studies. Will use a random increasing wait time between study attempts. See log files for more error info.
#' @return a data.table of all metadata for all accessions
pipeline_metadata <- function(accessions, config, force = FALSE, max_attempts = 7, max_attempts_loop = 1,
                              BPPARAM = if (length(accessions) > 5) {bpparam()} else SerialParam()) {
  do <- ifelse(config$mode == "online", "downloading", "processing")
  message("- Metadata ", do, "...")
  accessions <- gsub(" $|^ ", "", accessions)
  loop_attempts <- 0
  while (loop_attempts < max_attempts_loop) {
    if (loop_attempts > 0) {
      message("Unstable connection / SRA server is haltig queries (loop attempt: ", loop_attempts, ")")
      Sys.sleep(loop_attempts * 1.5)
    }
    all_SRA_metadata <- bplapply(accessions, function(study_accession, config, force) {
      # Avoid requirest overflow and still keep it fast
      overflooding_error <- "Error in open.connection(x, \"rb\") : HTTP error 429.\n"
      attempts <- 0
      save_file <- file.path(config[["metadata"]], paste0("SraRunInfo_", study_accession, ".csv"))
      while(attempts < max_attempts) {
        if ((attempts > 0 | force) & file.exists(save_file)) file.remove(save_file)
        if ((force | !file.exists(save_file))) Sys.sleep((attempts * 0.5) + runif(1, min = 0.1, max = 1))
        res <- try(
          cbind(study_accession = study_accession,
                ORFik::download.SRA.metadata(study_accession, abstract = "no",
                                             outdir = config[["metadata"]], auto.detect = TRUE))
          , silent = TRUE)
        if (!is(res, "try-error")) {
          break
        } else attempts <- attempts + 1
      }
      print(attempts)
      if (is(res, "try-error")) {
        if (res == overflooding_error) {
          warnings("Network overflow with ", attempts, " for ", study_accession)
        } else {
          warnings("Failed to access study:", study_accession)
          message("Here is error")
          print(res)
        }
        return(data.table())
      }
      safety_check <- c("LIBRARYTYPE", "REPLICATE", "STAGE", "CONDITION",
                        "INHIBITOR")
      are_there <- safety_check %in% colnames(res)
      if (!all(are_there)) {
        check_add <- data.table(LIBRARYTYPE = "", REPLICATE = "", STAGE = "",
                                CONDITION = "", INHIBITOR = "")
        res <- cbind(res, check_add[, ..are_there])
      }
      return(res)
    }, config = config, force = force, BPPARAM = BPPARAM)
    loop_attempts <- loop_attempts + 1
  }


  if (any(lengths(all_SRA_metadata) == 0)) {
    message("Studies with download errors:")
    print(accessions[lengths(all_SRA_metadata) == 0])
    stop("Unstable connection / SRA server is haltig queries, try again later!")
  }
  rbind_fill <- FALSE
  if (any(length(table(lengths(all_SRA_metadata))) != 1)) {
    if (all(length(table(lengths(all_SRA_metadata))) == 2)) {
      rbind_fill <- TRUE
    } else
      stop("Malformed data columns, debug or use 'force' = TRUE to redownload all")
  }

  message("-- Metadata ", do, " done")
  return(rbindlist(all_SRA_metadata, fill = rbind_fill))
}

pipeline_metadata_filter <- function(all_SRA_metadata, organisms = "all",
                                     LibraryLayouts = "SINGLE",
                                     LibraryStrategy = c("RNA-Seq", "miRNA-Seq", "OTHER"),
                                     Platforms = "ILLUMINA",
                                     removeAllNACols = TRUE) {
  stopifnot(all(LibraryLayouts %in% c("SINGLE", "PAIRED")))
  all_metadata_RFP <- all_SRA_metadata
  unique_organisms <- unique(all_SRA_metadata$ScientificName)
  if (any(organisms == "all")) {
    organisms <- unique_organisms
  } else if (!all(organisms %in% unique_organisms)) stop("You subseted to an organism not in existing in any study!")
  # Some experiments have mixed organisms:
  cat("-- Total samples to begin with", "\n")
  print(nrow(all_metadata_RFP))
  cat("-- All organisms to begin with", "\n")
  print(table(all_metadata_RFP$ScientificName))
  # Row filters:
  cat("-- Filtering unwanted sample classifications", "\n")
  filtered_RFP <- all_metadata_RFP[ScientificName %in% organisms &
                                     LibraryLayout %in% LibraryLayouts &
                                     Platform %in% Platforms,]

  if (length(grep("is currently private", filtered_RFP$sample_title)) > 0) {
    message("Private (not open for public) samples detected, they were removed!")
    filtered_RFP <- filtered_RFP[-grep("is currently private", sample_title),]
  }

  invalid_library_strategy <- !(filtered_RFP$LibraryStrategy %in% LibraryStrategy)
  if (any(invalid_library_strategy))
    filtered_RFP <- filtered_RFP[!invalid_library_strategy,]
  cat("-- Numer of samples filtered out", "\n")
  print(nrow(all_metadata_RFP) - nrow(filtered_RFP))
  if (nrow(filtered_RFP) == 0) stop("all metadata got filtered out! Check filtering conditions and try again.")

  # Now remove columns not wanted
  if (removeAllNACols) {
    # Remove all NA columns:
    empty_cols <- unlist(lapply(filtered_RFP, function(x)!all(is.na(x))))
    force_keep_cols <- c("sample_source", "sample_title", "CONDITION")
    stopifnot(all(names(force_keep_cols) %in% colnames(filtered_RFP)))
    empty_cols[force_keep_cols] <- TRUE
    filtered_RFP <- filtered_RFP[,which(empty_cols),with=F]
    if (is.logical(filtered_RFP$CONDITION)) {
      filtered_RFP$CONDITION <- as.character(filtered_RFP$CONDITION)
      filtered_RFP$CONDITION <- ""
    }
  }

  return(filtered_RFP)
}

pipeline_metadata_annotate <- function(filtered_RFP) {
  # Use ORFik auto detection of library type:
  # First check from sample
  filtered_RFP$LIBRARYTYPE <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                   ORFik:::libNames(), "auto")
  # Then check from source
  filtered_RFP[LIBRARYTYPE == "",]$LIBRARYTYPE <- ORFik:::findFromPath(filtered_RFP[LIBRARYTYPE == "",]$sample_source,
                                                                       ORFik:::libNames(), "auto")
  # Then check from type
  if (!is.null(filtered_RFP$LibraryName)) {
    filtered_RFP[LIBRARYTYPE == "",]$LIBRARYTYPE <-
      ORFik:::findFromPath(filtered_RFP[LIBRARYTYPE == "",]$LibraryName,
                           ORFik:::libNames(), "auto")
  }

  # Then check strategy
  filtered_RFP[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LIBRARYTYPE <-
    ORFik:::findFromPath(filtered_RFP[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LibraryStrategy,
                         ORFik:::libNames(), "auto")

  filtered_RFP$REPLICATE <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                 ORFik:::repNames(), "auto")

  filtered_RFP$BATCH <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                             ORFik:::batchNames(), "auto")
  filtered_RFP$TIMEPOINT <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                 ORFik:::stageNames(), "auto")
  filtered_RFP$TISSUE <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                              ORFik:::tissueNames(), "auto")
  # browser()
  filtered_RFP[TISSUE == "",]$TISSUE <- ORFik:::findFromPath(filtered_RFP[TISSUE == "",]$sample_source,
                                                             ORFik:::tissueNames(), "auto")
  filtered_RFP$CELL_LINE <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                 ORFik:::cellLineNames(), "auto")
  filtered_RFP[CELL_LINE == "",]$CELL_LINE <- ORFik:::findFromPath(filtered_RFP[CELL_LINE == "",]$sample_source,
                                                                   ORFik:::cellLineNames(), "auto")
  filtered_RFP[(TISSUE == "") & CELL_LINE != "",]$TISSUE <- ORFik:::findFromPath(filtered_RFP[(TISSUE == "") & CELL_LINE != "",]$CELL_LINE,
                                                                                 ORFik:::cellLineNames(TRUE), "auto")
  # Correct replicate names
  duplicated_samples <- any(duplicated(filtered_RFP$sample_title) & !is.na(filtered_RFP$sample_title))
  if (duplicated_samples) {
    filtered_RFP[, has_dups := any(duplicated(sample_title) | all(REPLICATE == "")), by= .(study_accession)]
    filtered_RFP[has_dups == TRUE, REPLICATE := seq(.N), by= .(study_accession, sample_title)]
    filtered_RFP$has_dups <- NULL
  }
  # Knock outs
  filtered_RFP$GENE <- ""
  combs <- expand.grid(c("Δ", "KO"), c(" .*", "_.*", "$"))
  KO_combs <- paste(paste0(combs[,1], combs[,2]), collapse = "|")
  KO_hits <- grep(KO_combs, filtered_RFP$sample_title)
  filtered_RFP[KO_hits,]$CONDITION <- "KO"
  filtered_RFP[KO_hits,]$GENE <- gsub(KO_combs, "", filtered_RFP[KO_hits,]$sample_title)
  # Exclude KO hits column
  filtered_RFP$CONDITION[-KO_hits] <- ORFik:::findFromPath(filtered_RFP$sample_title[-KO_hits],
                                                 ORFik:::conditionNames(), "auto")
  filtered_RFP$FRACTION <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                ORFik:::fractionNames(), "auto")
  filtered_RFP$INHIBITOR <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                 ORFik:::inhibitorNames(), "auto")
  filtered_RFP[LIBRARYTYPE == "" & INHIBITOR != "", LIBRARYTYPE := "RFP"]
  # Add columns needed for checks
  filtered_RFP[, `:=`(KEEP = "", UNIQUE = "", CHECKED = "", name = "", not_unique = "")]
  cat("Studies kept:", "\n")
  cat(length(unique(filtered_RFP$Submission)), "\n")
  cat("Library types detected", "\n")
  print(table(filtered_RFP$LIBRARYTYPE)) # It can't find all (bad info)
  filtered_RFP[]
  return(filtered_RFP)
}
#' Validate that metadata is ready to run
#'
#' Only checks validity of columns that are set to KEEP = TRUE.
#' Validity is check by seeing that for each study, per organism,
#' That there is a unique set of sample columns that identify that sample
#' @param dt a data.table of all metadata for all studies
#' @inheritParams curate_metadata
#' @param output_file a path to store final Ribo csv, default config$complete_metadata
#' @param backup_file a path to store backup csv, default: config$backup_metadata,
#' store all runs, even non Ribo-seq, such that it can be used if wanted.
pipeline_validate_metadata <- function(dt, config,
                                       accessions = NULL,
                                       only_curated = FALSE,
                                       output_file = config$complete_metadata,
                                       backup_file = config$backup_metadata,
                                       next_round_file = config$temp_metadata,
                                       libtypes = "RFP") {
  message("-- Validating metadata")
  message("- Total input samples: ", nrow(dt))
  if (nrow(dt) == 0) {
    message("Table is empty returning")
    return(FALSE)
  }
  add_backup(dt, backup_file)
  # Step1 Upload your csv to project folder
  #stopifnot(all(files$CHECKED %in% c(TRUE, FALSE)));stopifnot(all(files$DISCARD %in% c(TRUE, FALSE)));stopifnot(all(files$ASSIGNED_TO %in% c("PREETI", "KORNEL", "TESHOME", "HÅKON")))
  #files <- files[DISCARD == FALSE & CHECKED == TRUE,];
  files <- dt
  files[LIBRARYTYPE == "RPF", LIBRARYTYPE := "RFP"] # Force RFP naming for safety!
  # files <- files[KEEP %in% c(TRUE, NA),]; nrow(files)
  # files <- files[LIBRARYTYPE %in% c(libtypes),]
  if (config$mode == "local") {
    duplicate_stats <- files[,sum(duplicated(Run)), by = .(study_accession, ScientificName)]
    duplicates <- sum(duplicate_stats$V1)
    if (duplicates > 0) {
      message("You have duplicated runs in you metadata!")
      message("For mode = local, Runs are presumed to be file paths and should
             therefor never be duplicated within the (study, organism) pair!")
      print(duplicate_stats)
      fwrite(files, next_round_file)
      return(FALSE)
    }
  }
  if (config$mode == "online" && length(unique(files$Run)) != length(files$Run)) {
    message("You have duplicated runs in you metadata!")
    message("For mode = online, Runs are presumed to be SRR numbers and should
             therefor never be duplicated!")
    print(files[duplicated(files$Run),]$Run)
    fwrite(files, next_round_file)
    return(FALSE)
  }
  # Check if metadata is now valid
  files <- metadata_is_valid(files)
  if (!is.null(accessions) & only_curated) {
    files_to_check <- files[KEEP == TRUE & study_accession %in% accessions,]
  } else files_to_check <- files[KEEP == TRUE,]
  any_not_unique <- sum(files_to_check$not_unique) != 0
  if (any_not_unique) {
    message("Sorry, still not unique, try again, may the force be with you!")
    message("Number of not unique: ", sum(files_to_check$not_unique))
    top30 <- ifelse(sum(files_to_check$not_unique) > 30, "(showing first 30)", "")
    message("Run ids of non unique: ", top30)
    print(head(files_to_check[not_unique == TRUE, ]$Run, 30))
  } else {
    export_sucessful_metadata(files, libtypes, output_file,
                              next_round_file)
    return(TRUE)
  }
  # If failed, check why, and save in bottom and do a new round of manual
  files <- metadata_columns_cleanup(files)
  fwrite(files, next_round_file)
  return(FALSE)
}

metadata_is_valid <- function(files) {
  if ("GENE" %in% colnames(files)) {
    files[, name := paste(LIBRARYTYPE, BATCH, REPLICATE, TIMEPOINT, TISSUE, CELL_LINE,
                          CONDITION, GENE, FRACTION, INHIBITOR, sep = "_")]
  } else {
    files[, name := paste(LIBRARYTYPE, BATCH, REPLICATE, TIMEPOINT, TISSUE, CELL_LINE,
                          CONDITION, FRACTION, INHIBITOR, sep = "_")]
  }

  files[, not_unique := duplicated(name), by = .(BioProject, ScientificName)]
  return(files)
}

#' Add new data to the temp_metadata
add_new_data <- function(accessions, config, organisms = "all",
                         google_url = config$google_url,
                         complete_metadata = config$complete_metadata,
                         LibraryStrategy = c("RNA-Seq", "miRNA-Seq", "OTHER"),
                         LibraryLayouts = c("SINGLE", "PAIRED"), Platforms = "ILLUMINA",
                         open_editor = interactive(),
                         open_google_sheet = interactive(),
                         organism_name_cleanup = TRUE) {
  if (file.exists(config$blacklist)) {
    blacklist <- fread(config$blacklist)
    if (any(blacklist$id %in% accessions)) {
      print(blacklist$id[blacklist$id %in% accessions])
      stop("You tried to add an id in blacklist, see above for which!")
    }
  }
  # Step1: Get all metadata
  all_SRA_metadata <- pipeline_metadata(accessions, config)
  # Step2: Filter out what you do not want
  filtered_SRA_metadata <- pipeline_metadata_filter(all_SRA_metadata, organisms,
                                                    LibraryLayouts, LibraryStrategy,
                                                    Platforms)
  # Step3: Try to auto annotate
  all_SRA_metadata_RFP <- pipeline_metadata_annotate(filtered_SRA_metadata)
  # Step4: Fix organism names
  if (organism_name_cleanup) {
    all_SRA_metadata_RFP[, ScientificName := organism_name_cleanup(ScientificName)]
  }
  # Step5: Merge with existing finished metadata
  complete_metadata_dt <- data.table()
  new_studies_count <- nrow(all_SRA_metadata_RFP)
  if (!is.null(google_url)) {
    complete_metadata_dt <- read_sheet_safe(google_url)
  } else if (file.exists(complete_metadata)) complete_metadata_dt <- fread(config$temp_metadata)
  if (nrow(complete_metadata_dt) > 0) { # IF you have started before
    new_studies <- !(paste0(all_SRA_metadata_RFP$study_accession, "___",
                          all_SRA_metadata_RFP$ScientificName) %in%
                   paste0(complete_metadata_dt$study_accession, "___",
                          complete_metadata_dt$ScientificName))
    new_studies_count <- sum(new_studies)
    suppressWarnings(all_SRA_metadata_RFP[, STAGE := NULL])
    suppressWarnings(all_SRA_metadata_RFP[, Sex := NULL])
    all_SRA_metadata_RFP <- suppressWarnings(rbindlist(list(complete_metadata_dt,
                              all_SRA_metadata_RFP[new_studies,]), fill = TRUE))
  }
  message("Total samples before saving: ", nrow(all_SRA_metadata_RFP))
  message("New samples to annotate: ", new_studies_count)
  fwrite(all_SRA_metadata_RFP, config$temp_metadata)
  # Step6: Store as csv and open in google sheet and fix
  # Create a new sheet or use existing one
  #google_url <- gs4_create(name = "RFP_next_round_manual2.csv")
  google_url_defined <- !is.null(google_url)
  if (google_url_defined) {
    if (sum(new_studies) > 0) {
      local_google_copy <- file.path(config$project, "temp_google_local.csv")
      fwrite(all_SRA_metadata_RFP, local_google_copy)
      write_sheet(read.csv(local_google_copy),
                  ss = google_url,
                  sheet = 1)
    }
  }

  if(open_editor & interactive()) {
    if (google_url_defined & open_google_sheet & config$mode == "online") {
      browseURL(google_url)
    } else {
      message("Update data for unique rows and press syncronize,",
              " then press 'done'!")
      DataEditR::data_edit(config$temp_metadata, read_fun = "fread") %>%
        fwrite(config$temp_metadata)
    }
  }
  return(invisible(NULL))
}

match_bam_to_metadata <- function(bam_dir, study, paired_end) {
  bam_files <- ORFik:::findLibrariesInFolder(dir = bam_dir,
                                             types = c("bam", "bed", "wig", "ofst"),
                                             pairedEndBam = paired_end)

  bam_files_base <- ORFik:::remove.file_ext(bam_files, basename = T)
  matches <- match(study$Run, bam_files_base)
  names_have_info_extensions <- anyNA(matches)
  if (names_have_info_extensions) {
    # TODO: Make this all more clear and failproof
    bam_files_base <- gsub("_Aligned.*", "", bam_files_base)
    bam_files_base <- sub(".*trimmed_", "", bam_files_base)
    bam_files_base <- sub(".*trimmed2_", "", bam_files_base)
    bam_files_base <- sub(".*collapsed_", "", bam_files_base)
    stopifnot(all(bam_files_base != ""))
    matches <- match(study$Run, bam_files_base)
    if (anyNA(matches)) {
      bam_files_base <- gsub("_[0-9]$", "", bam_files_base)
      matches <- match(study$Run, bam_files_base)
      if (anyNA(matches)) {
        bam_files_base <- gsub(".*_", "", bam_files_base)
        matches <- match(study$Run, bam_files_base)
      }
    }

  }
  bam_files <- bam_files[matches]

  if (nrow(study) != length(bam_files)  | anyNA(bam_files)) {
    if ((nrow(study) > length(bam_files)) | anyNA(bam_files)) {
      message("Error: you are missing some bam files compared to study metadata")
      print(bam_files)
      stop("Missing bam files!")
    }  else stop("Not correct number of bam files in folder compared to study metadata")
  } else if (length(bam_files) != length(bam_files_base)) {
    warning("you have more bam files in folder compared to study metadata")
    print(bam_files)
  }
  if (length(bam_files) == 0) stop("Could not find SRR runs in aligned folder for: ")
  return(bam_files)
}

cleanup_metadata_for_exp <- function(study) {
  # Do some small correction to info and merge
  remove <- "^_|_$|^NA_|_NA$|^NA$|^_$|^__$|^___$"
  # Stage
  stage <- paste0(study$CELL_LINE, "_",study$TISSUE)
  stage <- gsub(paste0(remove, "|NONE_|_NONE"), "", stage)
  stage <- gsub(paste0(remove, "|NONE_|_NONE"), "", stage) # Twice
  # Condition
  condition <- paste0(study$CONDITION)
  if (!is.null(study$GENE)) condition <- paste(condition, study$GENE, sep = "_")
  condition <- gsub(remove, "", condition)
  condition <- gsub(remove, "", condition)
  condition <- gsub(remove, "", condition)
  condition <- gsub(remove, "", condition)
  # Fraction
  fraction <- paste(study$FRACTION,study$TIMEPOINT, study$BATCH, sep = "_")
  fraction <- gsub(remove, "", fraction)
  study$INHIBITOR[is.na(study$INHIBITOR)] <- ""
  add_inhibitor_to_fraction <-
    !all(study$INHIBITOR %in% c("chx", "CHX"))
  if (add_inhibitor_to_fraction) {
    fraction <- paste0(fraction, "_",study$INHIBITOR)
  }
  fraction <- gsub(remove, "", fraction); fraction <- gsub(remove, "", fraction)
  fraction <- gsub(remove, "", fraction); fraction <- gsub(remove, "", fraction)
  # PAIRED END
  paired_end <- study$LibraryLayout == "PAIRED"
  if (any(paired_end)) {
    message("Only running single end for now, make fix for this to work normally")
    paired_end <- FALSE
  }
  return(list(condition = condition, stage = stage, fraction = fraction,
              paired_end = paired_end))
}

add_backup <- function(dt, backup_file) {
  # Create a backup file for complete metadata
  # Only append new if already exists
  if (file.exists(backup_file)) {
    message("Appending new studies from current input to the backup file")
    backup <- fread(backup_file)
    new_backup <- rbindlist(list(backup,
                                 dt[!(study_accession %in% backup$study_accession),]), fill = TRUE)
    fwrite(new_backup, backup_file)
  } else fwrite(dt, backup_file)
  return(invisible(NULL))
}

metadata_columns_cleanup <- function(files) {
  # TODO: Split into cleanup and info
  # Check authors
  a <- unlist(files[, table(AUTHOR)])
  stopifnot(length(a[names(a) == ""]) == 0)
  # Check INHIBITOR
  files[INHIBITOR == "CHX", INHIBITOR := "chx"]
  files[INHIBITOR == "LTM", INHIBITOR := "ltm"]
  files[INHIBITOR == "Harringtonine", INHIBITOR := "harr"]
  a <- unlist(files[, table(INHIBITOR)])
  stopifnot(length(a[names(a) == ""]) == 0)
  files[INHIBITOR == "", INHIBITOR := "chx"]
  cat("# of Samples by Authors:\n")
  sort(table(files$AUTHOR), decreasing = TRUE) #unique(files[AUTHOR == "",]$BioProject)
  cat("# of Samples by Inhibitor:\n")
  table(files$INHIBITOR)
  cat("# of Samples by Condition:\n")
  table(files$CONDITION)
  table(files$FRACTION)
  table(files$TISSUE)#unique(files[TISSUE == "",]$BioProject)
  table(files$CELL_LINE)#unique(files[CELL_LINE == "Huh",]$BioProject)
  table(files$BATCH)
  return(files)
}

export_sucessful_metadata <- function(files, libtypes, output_file,
                                      next_round_file) {
  file_to_keep <- files[KEEP == TRUE & (LIBRARYTYPE %in% c(libtypes)),]
  message("Congratulations!")
  message("Data is unique!")
  message("Number of samples marked as KEEP: (",
          nrow(file_to_keep), "/", nrow(files), ")")
  message("Saving to: ", output_file)
  message("- Done")
  invalid_libtypes <- files[KEEP == TRUE & !(LIBRARYTYPE %in% c(libtypes)),]
  if (nrow(invalid_libtypes) > 0) {
    warning("Some samples were marked as keep with empty LIBRARYTYPE!",
            " Number of samples ignored this way: ", nrow(invalid_libtypes))

  }
  fwrite(files, next_round_file)
  fwrite(file_to_keep, output_file)
  return(invisible(NULL))
}

fill_with_random <- function(table_in, config, libtypes_keep = "all", checked_by = "auto",
                             start_auto_at_frac = max(suppressWarnings(max(as.integer(gsub("^auto_", "", grep("^auto_", table$FRACTION, value = TRUE)))) + 1), 1),
                             start_auto_at_time = max(suppressWarnings(max(as.integer(gsub("^auto_", "", grep("^auto_", table$TIMEPOINT, value = TRUE)))) + 1), 1),
                             upload_to_google = !is.null(config$google_url)) {
  stopifnot(("KEEP" %in% colnames(table_in)) && is(table_in$KEEP, "logical"))
  if (libtypes_keep == "all") libtypes_keep <- unique(table_in$LIBRARYTYPE)
  table <- copy(table_in)
  message("Samples that will be auto-named: ", nrow(table[LIBRARYTYPE %in% libtypes_keep & is.na(KEEP),]))
  table[LIBRARYTYPE %in% libtypes_keep & is.na(KEEP), CHECKED := checked_by]
  table[LIBRARYTYPE %in% libtypes_keep & is.na(KEEP), KEEP := TRUE]


  table <- metadata_is_valid(table)
  non_unique <- nrow(table[KEEP == TRUE & not_unique,]); non_unique

  if (non_unique > 0)
    suppressWarnings(table[LIBRARYTYPE %in% libtypes_keep & CHECKED == "auto" & not_unique & (is.na(FRACTION) || FRACTION == ""), FRACTION := paste0("auto_", seq(.N) + start_auto_at_frac)])
  table <- metadata_is_valid(table)
  non_unique <- nrow(table[KEEP == TRUE & not_unique,]); non_unique
  if (non_unique > 0)
    suppressWarnings(table[LIBRARYTYPE %in% libtypes_keep & CHECKED == "auto" & not_unique & (is.na(TIMEPOINT) || TIMEPOINT == ""), TIMEPOINT := paste0("auto_", seq(.N) + start_auto_at_time)])

  table <- metadata_is_valid(table)
  non_unique <- nrow(table[KEEP == TRUE & not_unique,]); non_unique
  if (non_unique > 0) warning("Some samples could not be auto named, check output and fix manually!")

  if (!is.null(config)) {
    message("Saving results to disc")
    local_google_copy <- file.path(config$project, "temp_google_local.csv")
    fwrite(table, config$temp_metadata)
    fwrite(table, local_google_copy)
    if (upload_to_google && !is.null(config$google_url)) {
      message("Saving results to google drive sheet")
      googlesheets4::write_sheet(read.csv(local_google_copy),
                                 ss = config$google_url,
                                 sheet = 1)
    }
  }

  return(table)
}

organism_name_cleanup <- function(organisms) {
  organisms[grepl("coronavirus 2|sars cov 2", organisms, ignore.case = TRUE)] <- "Sars cov2"

  organisms <- sub(" substr\\..*", "", organisms)
  organisms <- sub(" (K.12|BY4741|H37Rv|PAO1|Go1)", "", organisms, ignore.case = TRUE)
  organisms <- sub(" str\\..*", "", organisms)
  organisms <- sub(" subsp\\..*", "", organisms)
  # unique(organisms[sapply(gregexpr("\\s", organisms), function(x) sum(x > 0)) >= 2])
  return(organisms)
}


#' Create template study csv info
#' @param dir a directory with fast files
#'  (fasta or fastq, uncompressed/compressed)
#' @param exp_name experiment name, a prefix for directory name.
#'  The relative path without the -organism part.
#' @param files character vector of full paths,
#'  default: \code{list.files(dir, pattern = "\\.fast.*", full.names = TRUE)}
#' @param organism Scientific latin name of organism, give single if all
#' are equal, or specify for each file if multiple. Example: "Homo sapiens"
#' @param sample_title Define a valid sample name per sample, like,
#'  'WT_ribo_cyclohexamide_rep1'
#' @param paired character, either "SINGLE" or "PAIRED"
#' @param Model Sequencing machine name, default: "Illumina Genome Analyzer",
#' give single if all are equal, or specify for each file if multiple.
#' @param Platform Sequencing machine platform, default: "ILLUMINA",
#' give single if all are equal, or specify for each file if multiple.
#' @return a data.table with required columns:\cr
#' study_accession: the basename of folder for experiment\cr
#' RUN: The sample name (For paired end, still 1 row, with basename of sample)\cr
#' ...
#' @export
local_study_csv <- function(dir, exp_name, files = list.files(dir, pattern = "\\.fast.*", full.names = TRUE),
                            organism, sample_title, paired = FALSE,
                            Model = "Illumina Genome Analyzer",
                            Platform = "ILLUMINA",
                            LibraryName = "",
                            LibraryStrategy = "RNA-Seq",
                            LibrarySource = "TRANSCRIPTOMIC",
                            avgLength = NA,
                            SampleName = c(""), sample_source = c(""),
                            MONTH = substr(Sys.time(), 6,7),
                            YEAR = substr(Sys.time(), 1,4),
                            AUTHOR = Sys.info()["user"]) {
  #if (any(paired %in% "PAIRED")) stop("Only single end supported for now!")
  message("Using files:")
  print(basename(files))
  message("total files: ", length(files))
  stopifnot(length(files) > 0)
  stopifnot(all(paired %in% c("SINGLE", "PAIRED")))
  if (is.character(avgLength) && length(avgLength) == 1 && avgLength == "auto") {
    avgLength <- sapply(files, function(file) round(mean(width(readDNAStringSet(file, format = c("fasta", "fastq")[1 + grepl("\\.fastq|\\.fq", file)], nrec = 100)))))
  }

  file_basenames <- gsub("\\.fast.*", "", basename(files))
  file_basenames_split <- file_basenames
  file_basenames_unique <- unique(gsub("_\\.[1-9]\\.fast", "", file_basenames))
  size_MB <- floor(file.size(files)/1024^2)
  dt_temp <- data.table(Run = file_basenames_unique,
                        spots = NA, bases = NA, avgLength = avgLength,
                        size_MB, LibraryName, LibraryStrategy,
                        LibrarySource, LibraryLayout = paired,
                        Submission = exp_name,
                        Platform, Model, ScientificName = organism,
                        SampleName, MONTH, YEAR, AUTHOR, sample_source,
                        sample_title)
  dt_temp[, BioProject := as.integer(as.factor(Submission))]

  auto_detect <- TRUE # Always TRUE for now
  if (auto_detect) {
    dt_temp <- ORFik:::metadata.autnaming(dt_temp)
  }

  return(dt_temp)
}

