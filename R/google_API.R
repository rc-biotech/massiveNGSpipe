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

default_sheets <- function() {
  "https://docs.google.com/spreadsheets/d/18Y-WDvV_w0kTT3Xap4M5GZWpg39ZK-gUzzV7fuKbMvo/edit#gid=769582544"
}
