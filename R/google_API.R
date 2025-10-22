read_sheet_safe <- function(google_url, gargle_email = gargle_mail_safe()) {
  googledrive::drive_auth(email = gargle_email)

  file <- tempfile(fileext = ".csv")
  drive_download(google_url, path = file, type = "csv", overwrite = TRUE)
  sheet <- fread(file)
  col_classes <- sapply(sheet, class)
  if (any(col_classes == "list")) warning("You have columns that were returned as list,
                                          check your google doc!")
  sheet$CONDITION <- as.character(sheet$CONDITION); sheet[CONDITION == "NULL",]$CONDITION <- ""
  sheet$TIMEPOINT <- as.character(sheet$TIMEPOINT); sheet[TIMEPOINT == "NULL",]$TIMEPOINT <- ""
  sheet$CELL_LINE <- as.character(sheet$CELL_LINE); sheet[CELL_LINE == "NULL",]$CELL_LINE <- ""
  return(sheet)
}

#' Sync your local temp to google sheet
#' @param config a mNGSp config with defined google_url and existing
#'  temp_metadata file
#' @return invisible(NULL)
#' @export
sync_sheet_safe <- function(config, validate_to_complete_local = FALSE,
                            gargle_email = gargle_mail_safe()) {
  stopifnot(!is.null(config$google_url))
  stopifnot(!is.null(config$temp_metadata) & file.exists(config$temp_metadata))
  sheet <- read_sheet_safe(config$google_url, gargle_email)
  fwrite(sheet, config$temp_metadata)
  message("Google sheet synced to file: \n", config$temp_metadata)
  if (validate_to_complete_local) curate_metadata(NULL, config, google_url = NULL)
  return(invisible(NULL))
}

write_sheet_safe <- function(file, google_url, sheet = 1) {
  message("Uploading updated version to google sheet:")
  write_sheet_safe_new(file, google_url, sheet = sheet - 1)
}

gargle_mail_safe <- function(email = NULL) {
  if (!is.null(email)) return(email)
  if (!is.null(gargle::gargle_oauth_email())) {
    return(gargle::gargle_oauth_email())
  }
  emails <- unique(sub(".*_", "", list.files("~/.cache/gargle/")))
  if (length(emails) == 0) return(NULL)
  message("Auto selecting email for google api: ", emails[1])
  if (length(emails) > 1) {
    message("Alternatives:", paste(emails[-1], collapse = ", "))
    message("Change default by doing: options(gargle_oauth_email = 'your_email')")
  }
  options(gargle_oauth_email = emails[1])
  return(emails[1])
}

write_sheet_safe_new <- function(file, google_url, sheet = 0,
                                 tmp_name = sprintf("tmp_paste_%s", format(Sys.time(), "%Y%m%d-%H%M%S")),
                                 gargle_email = gargle_mail_safe()) {
  message("Uploading updated version to google sheet:")
  stopifnot(requireNamespace("googlesheets4"),
            requireNamespace("googledrive"),
            requireNamespace("data.table"),
            requireNamespace("httr"),
            requireNamespace("jsonlite"),
            requireNamespace("gargle"))

  if (!is(file, "data.frame")) file <- fread(file)
  df <- file
  message("- Estimated upload time: ", round((nrow(file) / 15e3) * 30, 0), " seconds..")

  target_ss <- google_url
  target_sheet <- sheet

  # Token + email setup
  googlesheets4::gs4_auth(email = gargle_email)
  googledrive::drive_auth(email = gargle_email)
  # One token with BOTH Drive + Sheets scopes
  tok <- gargle::token_fetch(scopes = c(
    "https://www.googleapis.com/auth/drive",
    "https://www.googleapis.com/auth/spreadsheets"
  ))
  suppressMessages(googledrive::drive_auth(token = tok, cache = TRUE))
  suppressMessages(googlesheets4::gs4_auth(token = tok))

  # 1) Create a temporary spreadsheet from CSV (INCLUDES HEADERS)
  tmp_csv <- tempfile(fileext = ".csv")
  data.table::fwrite(df, tmp_csv, sep = ",", na = "", quote = TRUE,
                     bom = FALSE, col.names = TRUE)
  tmp_file <- googledrive::drive_upload(
    media  = tmp_csv, name = tmp_name, type = "spreadsheet", verbose = FALSE
  )
  unlink(tmp_csv)

  target_ss_id <- as.character(googlesheets4::as_sheets_id(target_ss))
  tmp_ss_id    <- googledrive::as_id(tmp_file)

  # Resolve target tab (name | sheetId | 0-based index)
  props <- googlesheets4::sheet_properties(target_ss_id)
  resolve_target <- function(sheet, props) {
    if (is.character(sheet)) {
      stopifnot(sheet %in% props$name)
      list(id = props$id[match(sheet, props$name)], name = sheet)
    } else if (is.numeric(sheet) && length(sheet) == 1) {
      if (sheet %in% props$id) {
        nm <- props$name[match(sheet, props$id)]
        list(id = as.integer(sheet), name = nm)
      } else {
        idx_order <- order(props$index) # 0 = leftmost
        pos <- as.integer(sheet) + 1L
        stopifnot(pos >= 1L, pos <= length(idx_order))
        ix <- idx_order[pos]
        list(id = props$id[ix], name = props$name[ix])
      }
    } else stop("target_sheet must be name, sheetId, or 0-based index.")
  }
  tgt <- resolve_target(target_sheet, props)
  target_sheet_id   <- as.integer(tgt$id)
  target_sheet_name <- tgt$name
  current_rows <- as.integer(props$grid_rows[match(target_sheet_id, props$id)])
  current_cols <- as.integer(props$grid_columns[match(target_sheet_id, props$id)])

  # Temp sheetId (CSV->Sheet has one tab)
  tmp_props    <- googlesheets4::sheet_properties(tmp_ss_id)
  tmp_sheet_id <- as.integer(tmp_props$id[1])

  # Auth header
  access_token <- tok$credentials$access_token
  auth_header <- httr::add_headers(
    Authorization = sprintf("Bearer %s", access_token),
    `Content-Type` = "application/json"
  )
  gPOST <- function(url, body) {
    resp <- httr::RETRY("POST", url, auth_header,
                        body = jsonlite::toJSON(body, auto_unbox = TRUE),
                        times = 3, pause_base = 1, terminate_on = c(400,401,403,404))
    httr::stop_for_status(resp)
    ct <- httr::http_type(resp)
    if (identical(ct, "application/json")) {
      jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE)
    } else TRUE
  }

  # 2) Copy temp sheet into target spreadsheet (server-side)
  copy_to_url  <- sprintf("https://sheets.googleapis.com/v4/spreadsheets/%s/sheets/%s:copyTo",
                          as.character(tmp_ss_id), tmp_sheet_id)
  copy_res     <- gPOST(copy_to_url, list(destinationSpreadsheetId = target_ss_id))
  new_temp_sheet_id   <- as.integer(copy_res$sheetId)
  new_temp_sheet_name <- copy_res$title

  # 3) Resize grid if needed (only ENLARGE; never shrink to preserve existing formatting)
  need_rows <- nrow(df) + 1L   # +1 for header row
  need_cols <- ncol(df)
  enlarge_requests <- list()
  if (!is.na(current_rows) && need_rows > current_rows) {
    enlarge_requests <- c(enlarge_requests, list(list(updateSheetProperties = list(
      properties = list(sheetId = target_sheet_id,
                        gridProperties = list(rowCount = need_rows)),
      fields = "gridProperties.rowCount"
    ))))
  }
  if (!is.na(current_cols) && need_cols > current_cols) {
    enlarge_requests <- c(enlarge_requests, list(list(updateSheetProperties = list(
      properties = list(sheetId = target_sheet_id,
                        gridProperties = list(columnCount = need_cols)),
      fields = "gridProperties.columnCount"
    ))))
  }
  if (length(enlarge_requests)) {
    gPOST(sprintf("https://sheets.googleapis.com/v4/spreadsheets/%s:batchUpdate", target_ss_id),
          list(requests = enlarge_requests))
  }

  # 4) Clear values on target (keep formatting), then copyPaste HEADERS+DATA
  batch_url  <- sprintf("https://sheets.googleapis.com/v4/spreadsheets/%s:batchUpdate", target_ss_id)
  n_rows <- nrow(df) + 1L  # include header
  n_cols <- ncol(df)
  batch_body <- list(
    requests = list(
      # Clear values (big enough rectangle), preserves formats
      list(updateCells = list(
        range  = list(sheetId = target_sheet_id,
                      startRowIndex = 0, startColumnIndex = 0,
                      endRowIndex = max(current_rows, n_rows),
                      endColumnIndex = max(current_cols, n_cols)),
        fields = "userEnteredValue"
      )),
      # Copy headers+data from the imported temp tab (rows 0..n_rows)
      list(copyPaste = list(
        source = list(sheetId = new_temp_sheet_id,
                      startRowIndex = 0, startColumnIndex = 0,
                      endRowIndex = n_rows, endColumnIndex = n_cols),
        destination = list(sheetId = target_sheet_id,
                           startRowIndex = 0, startColumnIndex = 0),
        pasteType = "PASTE_VALUES",
        pasteOrientation = "NORMAL"
      )),
      # Delete the imported temp tab in the target
      list(deleteSheet = list(sheetId = new_temp_sheet_id))
    )
  )
  gPOST(batch_url, batch_body)

  # 5) Trash the temporary spreadsheet file
  suppressWarnings(googledrive::drive_trash(tmp_file))

  invisible(TRUE)
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
#' @param output_dir character, output directory, will be create if not existing
#' @param pattern character pattern of files, default "fastq\\.gz$"
#' @param email character, email of user, default: \email{hakontjeldnes@@gmail.com}
#' @param overwrite logical, default FALSE. If FALSE, delete existing files,
#' if TRUE dont download existing files.
#' @param recursive logical, default TRUE. Search sub folders for format files.
#' @return invisible(NULL)
#' @importFrom googledrive drive_auth as_id drive_ls drive_download
#' @export
download_fastq_google_drive <- function(drive_url,
                                        output_dir,
                                        files = google_drive_list_files(drive_url, email, format_pattern, recursive),
                                        format_pattern = "fastq\\.gz$",
                                        email = gargle_mail_safe(),
                                        overwrite = FALSE,
                                        recursive = TRUE) {
  # Script
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

#' List google files on drive wrapper
#' @inheritParams download_fastq_google_drive
#' @return a dribble of googl drive information
#' @export
google_drive_list_files <- function(drive_url, email = gargle_mail_safe(),
                                    format_pattern = "fastq\\.gz$", recursive = TRUE) {
  drive_auth(email = email)
  drive <- as_id(drive_url)
  files <- try(drive_ls(path = drive, pattern = format_pattern, recursive = recursive))

  if (is(files, "try-error")) stop("Folder does not exists on given google drive!")
  if (is.null(files) && nrow(files) == 0) stop("No files found in drive")
  return(files)
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
#' @param google_overwrite logical, default FALSE. If TRUE, if file exists on googledrive, it will replace.
#' Else it makes a copy.
#' @return invisible(NULL)
#' @export
plot_all_versions <- function(plot, file_prefix, formats = c("jpg", "svg"), send_to_google_drive = FALSE,
                              send_to_discord = FALSE, width = 8, height = 6, dpi = 600,
                              discord_connection = discord_connection_default_cached(),
                              google_drive_dir = google_drive_dir_links(1),
                              formats_discord = formats[!(formats %in% "svg")],
                              formats_google = formats[!(formats %in% c("jpg", "png"))],
                              preview_image = FALSE, discord_message = basename(file_prefix),
                              google_overwrite = FALSE) {
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
      googledrive::drive_upload(file, google_drive_dir, name = basename(file), overwrite = google_overwrite)
    }
    if (send_to_discord & (f %in% formats_discord)) {
      discordr::send_webhook_file(file, conn = discord_connection)
    }
  }

  if (send_to_discord & !is.null(discord_message) & is(discord_message, "character")) {
    discordr::send_webhook_message(discord_message)
  }
  return(invisible(NULL))
}

#' Get discord default webhook
#'
#' Steps to make the file:\cr
#' 1. Go to wanted discord channel\cr
#' 2. Create webhook\cr
#' 3. Call:\cr
#' \code{}  \cr
#' @param channel_id numeric, id of channel to use.
#' @param connection_file character, default: discord_cache_file(),
#' file path to cached webhooks.
#' @param verbose logical, FALSE. If TRUE, will display discordr info.
#' @return list of size 4, a webhook connection object or NULL if file does not exist
#' @import discordr
#' @export
discord_connection_default_cached <- function(channel_id = 1,
                                              connection_file = discord_cache_file(),
                                              verbose = FALSE) {
  if (is.null(getOption("default_discordr_connection"))) {
    if (identical(file.exists(connection_file), TRUE)) {
      conn <- suppressMessages(discordr::import_discord_connections(connection_file)[[channel_id]])
      discordr::set_default_discord_connection(conn)
    } else {
      conn <- NULL
    }
    return(conn)
  } else {getOption("default_discordr_connection")}
}

discord_cache_file <- function() {
  candidates <- file.path(c("~", "~/livemount"), ".cache/discordr/r_to_discord_config")
  candidates <- candidates[file.exists(candidates)]
  return(head(candidates, 1))
}

#' Given a webhook, create cache file
#'
#' Will append if webhooks exist in file already.
#' @param filepath character path, default "~/.cache/discordr/r_to_discord_config"
#' @return Returns the path you specified as input if sucessful.
#' @export
export_discord_cache_file <- function(filepath = "~/.cache/discordr/r_to_discord_config") {
  dir.create(dirname(filepath), recursive = TRUE)
  discordr::export_discord_connection(discordr::create_discord_connection(),
  filepath = filepath)
  return(filepath)
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
