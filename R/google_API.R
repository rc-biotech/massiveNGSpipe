read_sheet_safe <- function(google_url) {
  sheet <- as.data.table(googlesheets4::read_sheet(google_url))
  col_classes <- sapply(sheet, class)
  if (any(col_classes == "list")) warning("You have columns that were returned as list,
                                          check your google doc!")
  sheet$CONDITION <- as.character(sheet$CONDITION); sheet[CONDITION == "NULL",]$CONDITION <- ""
  sheet$TIMEPOINT <- as.character(sheet$TIMEPOINT); sheet[TIMEPOINT == "NULL",]$TIMEPOINT <- ""
  sheet$CELL_LINE <- as.character(sheet$CELL_LINE); sheet[CELL_LINE == "NULL",]$CELL_LINE <- ""
  return(sheet)
}

#' Get pre-existing or create a new google sheet for metadata
#' @inheritParams pipeline_config
#' @param id numeric, the row index of id. Default is 1.
#' @param sheet_name character, default basename(project_dir)
#' @return a string of full url.
default_sheets <- function(project_dir, id = 1,
                           sheet_name = basename(project_dir)) {
  if (is.null(id)) return(NULL)
  url_table <- file.path(project_dir, "google_api_url_table.csv")
  if (!file.exists(url_table)) {
    message("Creating new Google sheet:")
    new_sheet <- googlesheets4::gs4_create(sheet_name)
    new_sheet <- paste0("https://docs.google.com/spreadsheets/d/",
                        as.character(new_sheet))
  } else fread(url_table)[id,]$id
}

#' Get all fastq files from google drive
#'
#' @param drive_url character, url to drive folder (not individual files), example:
#' "https://drive.google.com/drive/folders/1N_Tyig3-IxMIe8Im24wBSNiS5vSohi91"
#' @param pattern character pattern of files, default "fastq\\.gz$"
#' @param output_dir character, output directory, will be create if not existing
#' @param email character, email of user, default: \email{hakontjeldnes@@gmail.com}
#' @param overwrite logical, default FALSE. If FALSE, delete existing files,
#' if TRUE dont download existing files.
#' @return invisible(NULL)
#' @export
download_fastq_google_drive <- function(drive_url = "https://drive.google.com/drive/folders/1N_Tyig3-IxMIe8Im24wBSNiS5vSohi91",
                                        format_pattern = "fastq\\.gz$",
                                        output_dir = "~/livemount/Bio_data/raw_data/autophagy_timecourse2/",
                                        email = "hakontjeldnes@gmail.com",
                                        overwrite = FALSE) {
  # Script
  drive_auth(email = email)
  drive <- as_id(drive_url)
  files <- try(drive_ls(path = drive, pattern = "fastq\\.gz$", recursive = FALSE))

  if (is(files, "try-error")) stop("Folder does not exists on given drive!")
  if (is.null(files) && nrow(files) == 0) stop("No files found in drive")

  message("Found valid files, tota of ", nrow(files), " files:")
  print(files)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("- Started downloading files")
  start <- Sys.time()
  message(start)
  for (file in seq(nrow(files))) {
    message("(", file, "/", nrow(files), ") ", files[file,]$name)
    out_full_path <- file.path(output_dir, files[file,]$name)
    if (file.exists(out_full_path) & !overwrite) {
      message("-- File already exists!")
    } else {
      status <-
        try(drive_download(files[file,]$id,  path = out_full_path,
                           overwrite = overwrite), silent = TRUE)
      if (is(status, "try-error")) {
        message("Download failed, error:")
        print(status)
        stop("")
      }
    }
  }
  message("Done")
  print(Sys.time() - start)
  return(invisible(NULL))
}
