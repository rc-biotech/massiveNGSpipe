#' Apply reviewed cell-level metadata fixes to a project's sheet
#'
#' For curation batches: each fix names a Run, a column, the value it expects to find (\code{old}) and the new
#' value. The sheet is read fresh (the online sheet is the source of truth); a fix whose cell no longer holds
#' \code{old} (changed by someone else meanwhile) is skipped and reported, never overwritten. Only the targeted
#' cells may change (verified), the upload goes through \code{write_sheet_safe()} and is read back, and every
#' applied cell is appended to \code{<project>/metadata_fix_log.csv} with the batch note.
#' @param fixes data.table with Run, column, old, new (and optionally reason)
#' @param project project folder name, e.g. "NGS_pipeline"
#' @param base_dir folder of the projects
#' @param note short batch description for the log (e.g. a decision file name)
#' @param use_google logical, default TRUE (else edit the local next_round_manual.csv)
#' @param dry_run logical, default TRUE
#' @return the fixes with a status column (applied / would_apply / changed_meanwhile / run_missing / ambiguous_run)
#' @export
sheet_apply_fixes <- function(fixes, project = "NGS_pipeline", base_dir = "~/livemount/Bio_data",
                              note = "", use_google = TRUE, dry_run = TRUE) {
  fixes <- copy(as.data.table(fixes))
  stopifnot(all(c("Run", "column", "old", "new") %in% names(fixes)), nrow(fixes) > 0,
            !anyDuplicated(fixes[, .(Run, column)]))
  norm <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }
  fixes[, `:=`(old = norm(old), new = norm(new))]
  bad <- fixes[grepl("^[0-9]*\\.[0-9]+$", new)]
  if (nrow(bad)) stop("Bare decimal values would be corrupted by the sheet's locale (e.g. ", bad$new[1],
                      "); store them with a unit, e.g. 0.5h, 0.3uM, Triton_0.1pct")
  sheet <- modality_read(project, base_dir, use_google)
  missing_cols <- setdiff(unique(fixes$column), names(sheet))
  if (length(missing_cols)) stop("Sheet has no column(s): ", paste(missing_cols, collapse = ", "))
  updated <- copy(sheet)
  runs <- as.character(sheet$Run)
  fixes[, status := NA_character_]
  for (k in seq_len(nrow(fixes))) {
    i <- which(runs == fixes$Run[k])
    if (!length(i)) { set(fixes, k, "status", "run_missing"); next }
    if (length(i) > 1) { set(fixes, k, "status", "ambiguous_run"); next }
    cur <- norm(updated[[fixes$column[k]]][i])
    if (cur != fixes$old[k]) { set(fixes, k, "status", "changed_meanwhile"); next }
    col <- fixes$column[k]
    value <- fixes$new[k]
    if (!is.character(updated[[col]])) set(updated, j = col, value = as.character(updated[[col]]))
    set(updated, i = i, j = col, value = value)
    set(fixes, k, "status", if (dry_run) "would_apply" else "applied")
  }
  print(fixes[, .N, by = status])
  if (dry_run || !any(fixes$status == "applied")) return(fixes[])
  # only the targeted cells may differ
  ok_cells <- fixes[status == "applied", paste(Run, column)]
  for (col in names(sheet)) {
    d <- which(norm(sheet[[col]]) != norm(updated[[col]]))
    if (!length(d)) next
    stray <- setdiff(paste(runs[d], col), ok_cells)
    if (length(stray)) stop("Unexpected changed cells: ", paste(head(stray, 5), collapse = "; "), "; nothing written")
  }
  modality_write(project, updated, base_dir, use_google)
  back <- modality_read(project, base_dir, use_google)
  a <- fixes[status == "applied"]
  got <- vapply(seq_len(nrow(a)), function(k) norm(back[[a$column[k]]][as.character(back$Run) == a$Run[k]])[1], character(1))
  if (nrow(back) != nrow(sheet) || any(got != a$new)) stop("Read-back of ", project, " does not show the applied fixes; check backups")
  log_file <- file.path(base_dir, project, "metadata_fix_log.csv")
  entry <- cbind(data.table(time = format(Sys.time(), tz = "UTC", usetz = TRUE), project = project, note = note),
                 a[, intersect(c("Run", "column", "old", "new", "reason"), names(a)), with = FALSE])
  fwrite(entry, log_file, append = file.exists(log_file))
  message(project, ": applied ", nrow(a), " fixes (log ", log_file, ")")
  fixes[]
}

#' Values a Google Sheet locale corrupts
#'
#' A sheet whose locale uses a decimal comma turns uploaded text like "0.5" into the number 0,5, and long decimals
#' like 0.3333333333 can lose their separator entirely ("3333333333"). This flags values that look already damaged
#' (\code{^[0-9]+,[0-9]+$}, 8+ digit numbers) and values at risk on upload (bare decimal numbers such as "0.5").
#' Store such values with a unit (0.5h, 0.3uM, Triton_0.1pct) so they stay text.
#' @param dt metadata table
#' @param columns annotation columns to check
#' @return data.table (Run, column, value, problem)
#' @export
metadata_locale_check <- function(dt, columns = c("CONDITION", "TIMEPOINT", "FRACTION", "BATCH", "REPLICATE", "GENE",
                                                 "TISSUE", "CELL_LINE", "INHIBITOR")) {
  dt <- as.data.table(dt)
  rbindlist(lapply(intersect(columns, names(dt)), function(cl) {
    v <- trimws(as.character(dt[[cl]])); v[is.na(v)] <- ""
    damaged <- grepl("^[0-9]+,[0-9]+$", v) | grepl("^[0-9]{8,}$", v)
    at_risk <- grepl("^[0-9]*\\.[0-9]+$", v)
    i <- which(damaged | at_risk)
    if (!length(i)) return(NULL)
    data.table(Run = as.character(dt$Run[i]), column = cl, value = v[i],
               problem = ifelse(damaged[i], "locale_damaged", "decimal_at_risk"))
  }))
}
