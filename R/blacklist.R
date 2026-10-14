#' Read the study / sample blacklist
#'
#' The blacklist (\code{config$blacklist}, default \code{<project>/BLACKLIST.csv}) lists accessions that are not what
#' their library strategy says, e.g. studies without Ribo-seq that match "ribosome profiling" style searches, or single
#' samples of a study that are not Ribo-seq. Columns: \code{id} (study accession: BioProject / SRA study / GEO series),
#' \code{run_id} (blank = the whole study is blacklisted; a Run accession = only that sample), \code{reason}.
#' Older files with only \code{id} are read as whole-study rows. The online sheet is the reference: with
#' \code{sync = TRUE} and a sheet URL (\code{config$blacklist_url}, else \code{<project>/blacklist_google_url.csv}),
#' the sheet is downloaded and the local CSV replaced first.
#' @param config pipeline config
#' @param sync logical, download the online sheet first (default TRUE when a URL is known)
#' @param url blacklist Google sheet URL or NULL
#' @return data.table(id, run_id, reason), run_id NA for whole-study rows
#' @export
blacklist_read <- function(config, sync = !is.null(url), url = blacklist_url(config)) {
  path <- config$blacklist
  if (isTRUE(sync) && !is.null(url)) {
    online <- tryCatch(read_sheet_safe(url), error = function(e) {
      warning("Blacklist sheet sync failed, using local copy: ", conditionMessage(e)); NULL })
    if (!is.null(online) && "id" %in% names(online)) write_csv_atomic(blacklist_normalize(online), path)
  }
  if (is.null(path) || !file.exists(path)) return(blacklist_normalize(data.table(id = character())))
  blacklist_normalize(fread(path, colClasses = "character"))
}

#' @noRd
blacklist_url <- function(config) {
  if (!is.null(config$blacklist_url)) return(config$blacklist_url)
  f <- if (!is.null(config$blacklist)) file.path(dirname(config$blacklist), "blacklist_google_url.csv") else ""
  if (nzchar(f) && file.exists(f)) fread(f)$id[1] else NULL
}

#' @noRd
blacklist_normalize <- function(x) {
  x <- as.data.table(x)
  for (cl in c("id", "run_id", "reason")) if (!cl %in% names(x)) x[, (cl) := NA_character_]
  x[, (names(x)) := lapply(.SD, function(v) { v <- trimws(as.character(v)); v[v %in% c("", "NA")] <- NA; v })]
  x <- x[!is.na(id) | !is.na(run_id)]
  x <- x[!duplicated(x[, .(id, run_id)])]
  setcolorder(x, c("id", "run_id", "reason"))
  x[]
}

#' Remove blacklisted studies and samples from a metadata table
#' @param dt metadata data.table with Run and study columns
#' @param blacklist from \code{blacklist_read()}
#' @param study_cols columns matched against whole-study ids
#' @return dt without blacklisted rows; attribute "removed" = number of rows removed
#' @export
blacklist_filter <- function(dt, blacklist, study_cols = c("BioProject", "study_accession", "SRAStudy")) {
  if (!nrow(blacklist) || !nrow(dt)) return(dt)
  studies <- blacklist[is.na(run_id) & !is.na(id)]$id
  runs <- blacklist[!is.na(run_id)]$run_id
  hit <- if ("Run" %in% names(dt)) dt$Run %in% runs else rep(FALSE, nrow(dt))
  for (cl in intersect(study_cols, names(dt))) hit <- hit | dt[[cl]] %in% studies
  out <- dt[!hit]
  setattr(out, "removed", sum(hit))
  out
}

#' Add studies or samples to the blacklist
#'
#' Syncs from the online sheet first (online is the reference), appends new rows (existing id/run_id pairs are
#' skipped), writes the local CSV and, with \code{upload = TRUE}, the online sheet (whole-tab rewrite, read back).
#' Removing the samples from the metadata sheet is a separate step.
#' @param config pipeline config
#' @param id study accession(s)
#' @param run_id NA (whole study) or Run accession(s), recycled with id
#' @param reason short text
#' @param upload logical, write the online sheet (requires permission to edit it)
#' @return the updated blacklist, invisibly
#' @export
blacklist_add <- function(config, id, run_id = NA_character_, reason = NA_character_, upload = FALSE) {
  url <- blacklist_url(config)
  cur <- blacklist_read(config, sync = !is.null(url), url = url)
  new <- blacklist_normalize(data.table(id = id, run_id = run_id, reason = reason))
  new <- new[!paste(id, run_id) %in% paste(cur$id, cur$run_id)]
  out <- rbind(cur, new, fill = TRUE)
  write_csv_atomic(out, config$blacklist)
  message("Blacklist: added ", nrow(new), " rows (", sum(is.na(new$run_id)), " studies, ",
          sum(!is.na(new$run_id)), " samples)")
  if (upload && nrow(new)) {
    if (is.null(url)) stop("No blacklist sheet URL known; set config$blacklist_url")
    local <- tempfile(fileext = ".csv"); fwrite(out, local)
    write_sheet_safe(local, url, sheet = 1)
    back <- blacklist_normalize(read_sheet_safe(url))
    if (nrow(back) != nrow(out)) stop("Blacklist read-back mismatch: ", nrow(back), " vs ", nrow(out), " rows")
  }
  invisible(out)
}
