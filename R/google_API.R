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
  } else fread()[id,]$id
}
