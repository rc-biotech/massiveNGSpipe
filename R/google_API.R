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


#' Export ggplots to multiple platforms
#'
#' Supports: Disc, google drive and discord
#' @param plot the ggplot
#' @param file_prefix the name of plot without the file extension
#' @param formats c("jpg", "svg")
#' @param send_to_google_drive = FALSE,
#' @param send_to_discord = FALSE
#' @param width = 8
#' @param height = 6
#' @param dpi = 600
#' @param discord_connection = discord_connection_default_cached()
#' @param google_drive_dir = google_drive_dir_links(1)
#' @param formats_discord = formats[!(formats %in% "svg")]
#' @param formats_google = formats[!(formats %in% c("jpg", "png"))]
#' @param preview_image = FALSE, if TRUE, will open browser and display image, then wait for you to press enter,
#' then you can cancel if ratios are wrong etc.
#' @param discord_message = NULL, if character, sends the message to discord after images are sent.
#' @return invisible(NULL)
#' @export
plot_all_versions <- function(plot, file_prefix, formats = c("jpg", "svg"), send_to_google_drive = FALSE,
                              send_to_discord = FALSE, width = 8, height = 6, dpi = 600,
                              discord_connection = discord_connection_default_cached(),
                              google_drive_dir = google_drive_dir_links(1),
                              formats_discord = formats[!(formats %in% "svg")],
                              formats_google = formats[!(formats %in% c("jpg", "png"))],
                              preview_image = FALSE, discord_message = NULL) {
  if (send_to_discord){
    if (is.null(discord_connection[[1]])) stop("Must have discord connection set!")
  }
  if (send_to_google_drive && is.null(google_drive_dir)) stop("Must have google drive folder connection set!")
  stopifnot(length(formats) > 0)

  for (f in formats) {
    message(f)
    file <- paste0(file_prefix, ".", f)
    ggsave(file, plot = plot,
           dpi = dpi, height = height, width = width)
    if (preview_image && f == formats[1]) {
      browseURL(file)
      readline(prompt="Press [enter] to continue if you are happy with image ratios")
    }

    if (send_to_google_drive & (f %in% formats_google)) {
      googledrive::drive_upload(file, google_drive_dir, name = basename(file))
    }
    if (send_to_discord & (f %in% formats_discord)) {
      discordr::send_webhook_file(file, conn = discord_connection)
    }
  }

  if (!is.null(discord_message) & is(discord_connection, "character")) {
    discordr::send_webhook_message(discord_message)
  }
  return(invisible(NULL))
}

#' Get discord default webhook
#' @param channel_id numeric, id of channel to use.
#' @param connection_file character, file path to cached webhooks
#' @param verbose logical, FALSE. If TRUE, will display discordr info.
#' @return list of size 4, a webhook connection object
#' @import discordr
discord_connection_default_cached <- function(channel_id = 1,
                                              connection_file = "~/livemount/.cache/discordr/r_to_discord_config",
                                              verbose = FALSE) {
  if (is.null(getOption("default_discordr_connection"))) {
    conn <- discordr::import_discord_connections(connection_file)[[channel_id]]
    discordr::set_default_discord_connection(conn)
    return(conn)
  } else {getOption("default_discordr_connection")}
}

google_drive_dir_links <- function(id = 1, id_file = "~/livemount/.cache/gargle/drive_folder_links") {
  stopifnot(length(id) > 0)

  dt <- fread(id_file)
  ids <- id
  if (is.numeric(id)) {
    if (!(id %in% seq_along(dt$id))) stop("numeric id > max ids")
    res <- dt[ids,]$link
    names(res) <- dt[ids,]$id
    return(res)
  }
  if (is.character(id)) {
    if (!(id %in% dt$id)) {
      print(dt$id)
      stop("ID given is not in list, valid ids given above")
      return(dt[id == ids,])
    }
  } else stop("id must be numeric (index) or character (id name)!")
}
