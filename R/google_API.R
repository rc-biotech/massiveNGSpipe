read_sheet_safe <- function(google_url) {
  sheet <- as.data.table(googlesheets4::read_sheet(google_url))
  sheet$TIMEPOINT <- as.character(sheet$TIMEPOINT); sheet[TIMEPOINT == "NULL",]$TIMEPOINT <- ""
  sheet$CELL_LINE <- as.character(sheet$CELL_LINE); sheet[CELL_LINE == "NULL",]$CELL_LINE <- ""
  return(sheet)
}
