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
#' @inheritParams run_pipeline
#' @param organisms character vector, default "all" (Use all organisms found).
#' Else binomial latin name with capital letter for genus: "Homo sapiens" etc.
#' @param google_url character, default our google doc sheet. If not NULL,
#' will use the google sheet to check for updated metadata.
#' @param complete_metadata path, default file.path(project_dir, "RFP_FINAL_LIST.csv").
#' The list of final candidates that are checked and have unique rows per bioproject
#' @param LibraryLayouts character vector, default c("SINGLE", "PAIRED"),
#' either or both of: c("SINGLE", "PAIRED")
#' @param Platforms character vector, default: "ILLUMINA". The sequencer technologies allowed.
#' @param step_mode logical, default FALSE. If TRUE, actives browser() and lets you go through
#' each step in debug mode for full control.
#' @return logical, TRUE if you had unique rows per accession per organism
#' @import googlesheets4 data.table ORFik BiocParallel fs stringr
#' @export
curate_metadata <- function(accessions, config, organisms = "all",
                            google_url = "https://docs.google.com/spreadsheets/d/18Y-WDvV_w0kTT3Xap4M5GZWpg39ZK-gUzzV7fuKbMvo/edit#gid=769582544",
                            complete_metadata = file.path(project_dir, "RFP_FINAL_LIST.csv"),
                            LibraryLayouts = c("SINGLE", "PAIRED"), Platforms = "ILLUMINA",
                            step_mode = FALSE) {
  #metadata_done <- file.exists(complete_metadata)
  if (step_mode) browser()
  if (!interactive()) stop("In non interactive mode you can not continue without RFP_FINAL_LIST.csv!")
  # Step1: Get all metadata
  all_SRA_metadata <- pipeline_metadata(accessions, config)
  # Step2: Try to auto annotate
  all_SRA_metadata_RFP <- pipeline_metadata_annotate(all_SRA_metadata, organisms,
                                                     LibraryLayouts, Platforms)
  complete_metadata_dt <- data.table()
  if (!is.null(google_url)) {
    complete_metadata_dt <- read_sheet_safe(google_url)
  } else if (file.exists(complete_metadata)) complete_metadata_dt <- fread(complete_metadata)
  if (nrow(complete_metadata_dt) > 0) { # IF you have started before
    new_studies <- !(all_SRA_metadata_RFP$study_accession %in%
                       complete_metadata_dt$study_accession)
    message("New samples to annotate: ", sum(new_studies))
    suppressWarnings(all_SRA_metadata_RFP[, STAGE := NULL])
    suppressWarnings(all_SRA_metadata_RFP[, Sex := NULL])
    all_SRA_metadata_RFP <- suppressWarnings(rbind(complete_metadata_dt,
                                                   all_SRA_metadata_RFP[new_studies,]))
  }
  fwrite(all_SRA_metadata_RFP, file.path(config$project, "RFP_pre_manual_annotation.csv"))
  # Step3: Store as csv and open in google sheet and fix
  # Create a new sheet or use existing one
  #google_url <- gs4_create(name = "RFP_next_round_manual2.csv")
  if (!is.null(google_url)) {
    write_sheet(read.csv((file.path(config$project, "RFP_pre_manual_annotation.csv"))),
                ss = google_url,
                sheet = 1); browseURL(google_url)
  }
  # Step4: Now, Check if it is valid (if not repeat step with new csv)
  while(TRUE) {
    readline(prompt = "You think metadata is ready?\n Press enter when ready: ")
    if (!is.null(google_url)) {
      message("- Reading google sheet")
      sheet <- read_sheet_safe(google_url)
    } else sheet <- fread(file.path(config$project, "RFP_next_round_manual.csv"))
    finished <- pipeline_validate_metadata(sheet, config)
    # TODO: Add upload of finished also to google sheet
    if (!finished) {
      if (!is.null(google_url)) {
        message("Uploading updated version to google sheet:")
        write_sheet(read.csv(file.path(config$project, "RFP_next_round_manual.csv")),
                    ss = google_url,
                    sheet = 1)
      } else message("Fix your local csv metadata file")
    }
  }
  #}
  return(TRUE)
}

#' Download all metadata for all accessions
#' @inheritParams curate_metadata
#' @param force logical, default FALSE. Force redownload of all, if something failed
#' @return a data.table of all metadata for all accessions
pipeline_metadata <- function(accessions, config, force = FALSE) {
  accessions <- gsub(" $|^ ", "", accessions)
  all_SRA_metadata <- bplapply(accessions, function(study_accession, config, force) {
    # Avoid requirest overflow and still keep it fast
    overflooding_error <- "Error in open.connection(x, \"rb\") : HTTP error 429.\n"
    attempts <- 0
    save_file <- file.path(config[["metadata"]], paste0("SraRunInfo_", study_accession, ".csv"))
    while(attempts < 5) {
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



pipeline_metadata_annotate <- function(all_SRA_metadata, organisms = "all",
                                       LibraryLayouts = "SINGLE",
                                       Platforms = "ILLUMINA") {
  stopifnot(all(LibraryLayouts %in% c("SINGLE", "PAIRED")))
  all_metadata_RFP <- all_SRA_metadata
  unique_organisms <- unique(all_SRA_metadata$ScientificName)
  if (organisms == "all") {
    organisms <- unique_organisms
  } else if (!all(organisms %in% unique_organisms)) stop("You subseted to an organism not in existing in any study!")
  # Some experiments have mixed organisms:
  cat("-- Total samples to begin with", "\n")
  print(nrow(all_metadata_RFP))
  cat("-- All organisms to begin with", "\n")
  print(table(all_metadata_RFP$ScientificName))
  # Row filters:
  filtered_RFP <- all_metadata_RFP[ScientificName %in% organisms &
                                   LibraryLayout %in% LibraryLayouts &
                                   Platform %in% Platforms,]
  nrow(filtered_RFP)
  if (length(grep("is currently private", filtered_RFP$sample_title)) > 0)
    filtered_RFP <- filtered_RFP[-grep("is currently private", sample_title),]
  nrow(filtered_RFP)
  if (length(grep("miRNA-Seq", filtered_RFP$LibraryStrategy)) > 0)
    filtered_RFP <- filtered_RFP[-grep("miRNA-Seq", LibraryStrategy),]
  nrow(filtered_RFP)
  # Remove all NA columns:
  filtered_RFP <- filtered_RFP[,which(unlist(lapply(filtered_RFP, function(x)!all(is.na(x))))),with=F]
  # Remove unwanted columns
  filtered_RFP[, Tumor := NULL]; filtered_RFP[, spots_with_mates := NULL]
  suppressWarnings(filtered_RFP[, Subject_ID := NULL])
  # Use ORFik auto detection of library type:
  # First check from sample
  filtered_RFP$LIBRARYTYPE <- ORFik:::findFromPath(filtered_RFP$sample_title,
                                                   ORFik:::libNames(), "auto")
  # Then check from source
  filtered_RFP[LIBRARYTYPE == "",]$LIBRARYTYPE <- ORFik:::findFromPath(filtered_RFP[LIBRARYTYPE == "",]$sample_source,
                                                                       ORFik:::libNames(), "auto")
  # Then check from type
  filtered_RFP[LIBRARYTYPE == "",]$LIBRARYTYPE <- ORFik:::findFromPath(filtered_RFP[LIBRARYTYPE == "",]$LibraryName,
                                                                       ORFik:::libNames(), "auto")
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
#' @param output_file a path to store final Ribo csv
#' @param backup_file a path to store backup csv, store all runs, even non Ribo-seq, such that it can be used if wanted.
pipeline_validate_metadata <- function(dt, config,
                                       output_file = file.path(config$project, "RFP_FINAL_LIST.csv"),
                                       backup_file = file.path(config$project, "ALL_BACKUP_LIST.csv")) {
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
  files <- files[KEEP %in% c(TRUE, NA),]; nrow(files); table(files$LIBRARYTYPE)
  files[LIBRARYTYPE == "RPF", LIBRARYTYPE := "RFP"] # Force RFP naming for safety!
  files <- files[LIBRARYTYPE %in% c("RFP"),]; nrow(files) # unique(files[LIBRARYTYPE == "",]$BioProject)
  stopifnot(length(unique(files$Run)) == length(files$Run)) #files[duplicated(files$Run),]$Run

  if ("GENE" %in% colnames(files)) {
    files[, name := paste(LIBRARYTYPE, BATCH, REPLICATE, TIMEPOINT, TISSUE, CELL_LINE,
                          CONDITION, GENE, FRACTION, INHIBITOR, sep = "_")]
  } else {
    files[, name := paste(LIBRARYTYPE, BATCH, REPLICATE, TIMEPOINT, TISSUE, CELL_LINE,
                          CONDITION, FRACTION, INHIBITOR, sep = "_")]
  }

  files[, not_unique := duplicated(name), by = BioProject]
  #
  if (sum(files$not_unique) != 0) {
    warning("Sorry, still not unique, try again, may the force be with you!")
  } else {
    message("Congratulations!")
    message("Data is unique!")
    message("Saving to: ", output_file)
    message("- Done")
    fwrite(files, output_file)
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
  fwrite(files, file.path(project_dir, "RFP_next_round_manual.csv"))
  return(FALSE)
}
