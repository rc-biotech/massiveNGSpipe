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
#' @inheritParams run_pipeline
#' @param organisms character vector, default "all" (Use all organisms found).
#' Else binomial latin name with capital letter for genus: "Homo sapiens" etc.
#' @param google_url character, default config$google_url. If not NULL,
#' will use the google sheet to check for updated metadata.
#' @param complete_metadata path, default file.path(project_dir, "RFP_FINAL_LIST.csv").
#' The list of final candidates that are checked and have unique rows per bioproject
#' @param LibraryLayouts character vector, default c("SINGLE", "PAIRED"),
#' either or both of: c("SINGLE", "PAIRED")
#' @param Platforms character vector, default: "ILLUMINA". The sequencer technologies allowed.
#' @param step_mode logical, default FALSE. If TRUE, actives browser() and lets you go through
#' each step in debug mode for full control.
#' @param fix_loop logical, default TRUE. If TRUE, will create a loop in
#' the end that runs until your metadata is valid. If you
#' only want to updata csv / google docs with new data and not try to complete
#' now, set this to FALSE.
#' @return logical, TRUE if you had unique rows per accession per organism
#' @import googlesheets4 data.table ORFik BiocParallel fs stringr
#' @export
curate_metadata <- function(accessions, config, organisms = "all",
                            google_url = config$google_url,
                            complete_metadata = config$complete_metadata,
                            LibraryLayouts = c("SINGLE", "PAIRED"), Platforms = "ILLUMINA",
                            step_mode = FALSE, open_google_sheet = interactive(),
                            fix_loop = TRUE) {
  if (!interactive() & step_mode)
    stop("In non interactive mode you can not run step_mode = TRUE!")
  stopifnot(length(fix_loop) == 1); stopifnot(is(fix_loop, "logical"))
  if (step_mode) browser()
  add_new_accessions <- !is.null(accessions)
  if (add_new_accessions) {
    add_new_data(accessions, config, organisms,
                 google_url, complete_metadata,
                 LibraryLayouts, Platforms,
                 open_google_sheet)
  }

  # Step6: Now, Check if it is valid (if not repeat step with new csv)
  while(fix_loop) {
    readline(prompt = "You think metadata is ready?\n Press enter when ready: ")
    if (!is.null(google_url)) {
      message("- Reading google sheet")
      sheet <- read_sheet_safe(google_url)
    } else sheet <- fread(config$temp_metadata)
    finished <- pipeline_validate_metadata(sheet, config)
    if (!is.null(google_url)) {
      message("Uploading updated version to google sheet:")
      write_sheet(read.csv(config$temp_metadata),
                  ss = google_url,
                  sheet = 1)
    } else message("Updated your local csv metadata file")
    if (finished) break
  }
  #}
  return(fix_loop)
}

#' Download all metadata for all accessions
#' @inheritParams curate_metadata
#' @param force logical, default FALSE. Force redownload of all, if something failed
#' @param max_attempts numeric, default 7. How many attempts to download beforing
#' failing? This call flood SRA with calls, so set it higher if you use a lot
#' of samples.
#' @return a data.table of all metadata for all accessions
pipeline_metadata <- function(accessions, config, force = FALSE, max_attempts = 7) {
  accessions <- gsub(" $|^ ", "", accessions)
  all_SRA_metadata <- bplapply(accessions, function(study_accession, config, force) {
    # Avoid requirest overflow and still keep it fast
    overflooding_error <- "Error in open.connection(x, \"rb\") : HTTP error 429.\n"
    attempts <- 0
    save_file <- file.path(config[["metadata"]], paste0("SraRunInfo_", study_accession, ".csv"))
    while(attempts < max_attempts) {
      if ((attempts > 0 | force) & file.exists(save_file)) file.remove(save_file)
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
    safety_check <- c("LIBRARYTYPE", "REPLICATE", "STAGE", "CONDITION", "INHIBITOR")
    are_there <- safety_check %in% colnames(res)
    if (!all(are_there)) {
      check_add <- data.table(LIBRARYTYPE = "", REPLICATE = "", STAGE = "",
                              CONDITION = "", INHIBITOR = "")
      res <- cbind(res, check_add[, ..are_there])
    }
    return(res)
  }
    , config = config, force = force)
  if (any(length(table(lengths(all_SRA_metadata))) != 1))
    stop("Malformed data, use 'force' = TRUE")
  return(rbindlist(all_SRA_metadata))
}

pipeline_metadata_filter <- function(all_SRA_metadata, organisms = "all",
                                     LibraryLayouts = "SINGLE",
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

  if (length(grep("is currently private", filtered_RFP$sample_title)) > 0)
    filtered_RFP <- filtered_RFP[-grep("is currently private", sample_title),]
  if (length(grep("miRNA-Seq", filtered_RFP$LibraryStrategy)) > 0)
    filtered_RFP <- filtered_RFP[-grep("miRNA-Seq", LibraryStrategy),]
  cat("-- Numer of samples filtered out", "\n")
  print(nrow(all_metadata_RFP) - nrow(filtered_RFP))
  # Now remove columns not wanted
  if (removeAllNACols) {
    # Remove all NA columns:
    filtered_RFP <- filtered_RFP[,which(unlist(lapply(filtered_RFP, function(x)!all(is.na(x))))),with=F]
    # Keep sample_title and sample_source
    if (is.null(filtered_RFP$sample_source)) filtered_RFP$sample_source <- ""
    if (is.null(filtered_RFP$sample_title)) filtered_RFP$sample_title <- ""
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
  filtered_RFP$BATCH <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                             ORFik:::batchNames(), "auto")
  filtered_RFP$REPLICATE <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                 ORFik:::repNames(), "auto")
  filtered_RFP$TIMEPOINT <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                 ORFik:::stageNames(), "auto")
  filtered_RFP$TISSUE <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                              ORFik:::tissueNames(), "auto")
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
  #filtered_RFP <- filtered_RFP[,-c(2:5, 10:15, 22:23, 26, 27, 28)]
  filtered_RFP[]
  return(filtered_RFP)
}
#' Validate that metadata is ready to run
#'
#' @param dt a data.table of all metadata for all studies
#' @inheritParams curate_metadata
#' @param output_file a path to store final Ribo csv, default config$complete_metadata
#' @param backup_file a path to store backup csv, default: config$backup_metadata,
#' store all runs, even non Ribo-seq, such that it can be used if wanted.
pipeline_validate_metadata <- function(dt, config,
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
  project_dir <- config$project
  if (file.exists(backup_file)) {
    message("Appending new studies from current input to the backup file")
    backup <- fread(backup_file)
    new_backup <- rbind(backup, dt[!(study_accession %in% backup$study_accession),])
    fwrite(new_backup, backup_file)
  } else fwrite(dt, backup_file)
  # Step1 Upload your csv to project folder
  #stopifnot(all(files$CHECKED %in% c(TRUE, FALSE)));stopifnot(all(files$DISCARD %in% c(TRUE, FALSE)));stopifnot(all(files$ASSIGNED_TO %in% c("PREETI", "KORNEL", "TESHOME", "HÅKON")))
  #files <- files[DISCARD == FALSE & CHECKED == TRUE,];
  files <- dt
  files[LIBRARYTYPE == "RPF", LIBRARYTYPE := "RFP"] # Force RFP naming for safety!
  # files <- files[KEEP %in% c(TRUE, NA),]; nrow(files)
  # files <- files[LIBRARYTYPE %in% c(libtypes),]
  if (length(unique(files$Run)) != length(files$Run)) {
    message("You have duplicated runs in you metadata!")
    print(files[duplicated(files$Run),]$Run)
    fwrite(files, next_round_file)
    return(FALSE)
  }
  files <- metadata_is_valid(files)
  #
  any_not_unique <- sum(files[KEEP == TRUE,]$not_unique) != 0
  if (any_not_unique) {
    message("Sorry, still not unique, try again, may the force be with you!")
  } else {
    message("Congratulations!")
    message("Data is unique!")
    message("Saving to: ", output_file)
    message("- Done")
    fwrite(files, next_round_file)
    fwrite(files[KEEP == TRUE & (LIBRARYTYPE %in% c(libtypes)),], output_file)
    return(TRUE)
  }

  # If failed, check why, and save in bottom and do a new round of manual
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
  # Check batch
  #files[BATCH == "1", BATCH := "b1"]; files[BATCH == "2", BATCH := "b2"]; files[BATCH == "3", BATCH := "b3"];files[BATCH == "4", BATCH := "b4"]
  #stopifnot(!any(files$batch %in% c(as.character(seq(5)))))
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

  files[, not_unique := duplicated(name), by = BioProject]
  return(files)
}

add_new_data <- function(accessions, config, organisms = "all",
                         google_url = config$google_url,
                         complete_metadata = config$complete_metadata,
                         LibraryLayouts = c("SINGLE", "PAIRED"), Platforms = "ILLUMINA",
                         open_google_sheet = interactive()) {
  # Step1: Get all metadata
  all_SRA_metadata <- pipeline_metadata(accessions, config)
  # Step2: Filter out what you do not want
  filetered_SRA_metadata <- pipeline_metadata_filter(all_SRA_metadata, organisms,
                                                     LibraryLayouts, Platforms)
  # Step3: Try to auto annotate
  all_SRA_metadata_RFP <- pipeline_metadata_annotate(filetered_SRA_metadata)
  # Step4: Merge with existing finished metadata
  complete_metadata_dt <- data.table()
  new_studies_count <- nrow(all_SRA_metadata_RFP)
  if (!is.null(google_url)) {
    complete_metadata_dt <- read_sheet_safe(google_url)
  } else if (file.exists(complete_metadata)) complete_metadata_dt <- fread(complete_metadata)
  if (nrow(complete_metadata_dt) > 0) { # IF you have started before
    new_studies <- !(all_SRA_metadata_RFP$study_accession %in%
                       complete_metadata_dt$study_accession)
    new_studies_count <- sum(new_studies)
    suppressWarnings(all_SRA_metadata_RFP[, STAGE := NULL])
    suppressWarnings(all_SRA_metadata_RFP[, Sex := NULL])
    all_SRA_metadata_RFP <- suppressWarnings(rbindlist(list(complete_metadata_dt,
                              all_SRA_metadata_RFP[new_studies,]), fill = TRUE))
  }
  message("New samples to annotate: ", new_studies_count)
  fwrite(all_SRA_metadata_RFP, file.path(config$project, "RFP_pre_manual_annotation.csv"))
  # Step5: Store as csv and open in google sheet and fix
  # Create a new sheet or use existing one
  #google_url <- gs4_create(name = "RFP_next_round_manual2.csv")
  if (!is.null(google_url)) {
    if (sum(new_studies) > 0) {
      write_sheet(read.csv((file.path(config$project, "RFP_pre_manual_annotation.csv"))),
                  ss = google_url,
                  sheet = 1)
    }

    if (interactive() & open_google_sheet) browseURL(google_url)
  }
  return(invisible(NULL))
}

match_bam_to_metadata <- function(bam_dir, study, paired_end) {
  bam_files <- ORFik:::findLibrariesInFolder(dir = bam_dir,
                                             types = c("bam", "bed", "wig", "ofst"),
                                             pairedEndBam = paired_end)

  bam_files_base <- ORFik:::remove.file_ext(bam_files, basename = T)
  bam_files_base <- gsub("_Aligned.*", "", bam_files_base)
  bam_files_base <- sub(".*trimmed_", "", bam_files_base)
  bam_files_base <- gsub(".*_", "", bam_files_base)
  stopifnot(all(bam_files_base != ""))
  bam_files <- bam_files[match(study$Run, bam_files_base)]

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
