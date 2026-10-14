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

#' Remove samples (rows) from a project's metadata sheet
#'
#' For blacklisted samples. The sheet is read fresh (online is the reference), the rows of \code{runs} are dropped,
#' all other rows must stay identical (verified), the upload goes through \code{write_sheet_safe()} (local backup in
#' \code{<project>/metadata_backups}) and is read back. Removed rows are saved to
#' \code{<project>/metadata_removed_rows.csv} and logged in \code{metadata_fix_log.csv}.
#' @param runs Run accessions to remove
#' @inheritParams sheet_apply_fixes
#' @return data.table(Run, status), status removed / would_remove / run_missing
#' @export
sheet_remove_runs <- function(runs, project = "NGS_pipeline", base_dir = "~/livemount/Bio_data",
                              note = "", use_google = TRUE, dry_run = TRUE) {
  runs <- unique(as.character(runs))
  sheet <- modality_read(project, base_dir, use_google)
  hit <- as.character(sheet$Run) %in% runs
  res <- data.table(Run = runs, status = ifelse(runs %in% sheet$Run, if (dry_run) "would_remove" else "removed", "run_missing"))
  print(res[, .N, by = status])
  if (dry_run || !any(hit)) return(res[])
  updated <- sheet[!hit]
  as_txt <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
  for (cl in names(sheet)) stopifnot(identical(as_txt(updated[[cl]]), as_txt(sheet[[cl]][!hit])))
  modality_write(project, updated, base_dir, use_google)
  back <- modality_read(project, base_dir, use_google)
  if (nrow(back) != nrow(updated) || any(as.character(back$Run) %in% runs))
    stop("Read-back of ", project, " still shows removed runs or a different row count; check backups")
  removed <- cbind(data.table(time = format(Sys.time(), tz = "UTC", usetz = TRUE), note = note), sheet[hit])
  f <- file.path(base_dir, project, "metadata_removed_rows.csv")
  fwrite(removed, f, append = file.exists(f))
  log_file <- file.path(base_dir, project, "metadata_fix_log.csv")
  fwrite(data.table(time = removed$time, project = project, note = note, Run = as.character(sheet$Run[hit]),
                    column = "ROW", old = "present", new = "removed", reason = "blacklisted sample"),
         log_file, append = file.exists(log_file))
  message(project, ": removed ", sum(hit), " rows (copies in ", f, ")")
  res[]
}

#' Drop samples from a processed ORFik experiment
#'
#' Rewrites the experiment CSV without \code{runs} (backup kept next to it) and removes the study's flags from the
#' first step that combines samples (\code{merged_lib}) onwards, so the next pipeline run rebuilds the study pool,
#' coverage formats and count tables from the remaining samples. Per-sample files of the dropped runs stay on disk.
#' @param exp experiment name
#' @param runs Run accessions to drop
#' @param config pipeline config
#' @param exp_dir directory of experiment CSVs
#' @param dry_run logical, default TRUE
#' @return number of samples dropped, invisibly
#' @export
experiment_remove_runs <- function(exp, runs, config, exp_dir = ORFik::config()["exp"], dry_run = TRUE) {
  exp_file <- file.path(exp_dir, paste0(exp, ".csv"))
  df <- read.experiment(exp_file, validate = FALSE, output.env = new.env())
  drop <- ORFik::runIDs(df) %in% runs
  message(exp, ": ", sum(drop), " of ", nrow(df), " samples to drop")
  if (dry_run || !any(drop)) return(invisible(sum(drop)))
  if (all(drop)) stop(exp, ": all samples would be dropped; remove the experiment instead")
  file.copy(exp_file, paste0(exp_file, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  lines <- readLines(exp_file)
  hdr <- grep('^"?libtype"?,', lines)[1]
  stopifnot(!is.na(hdr), length(lines) - hdr == nrow(df))
  keep_rows <- hdr + which(!drop)
  writeLines(lines[c(seq_len(hdr), keep_rows)], exp_file)
  check <- read.experiment(exp_file, validate = FALSE, output.env = new.env())
  stopifnot(nrow(check) == sum(!drop), !any(ORFik::runIDs(check) %in% runs))
  steps <- names(config$flag)
  steps <- steps[which(steps == "merged_lib"):length(steps)]
  remove_flag_all_exp(config, steps = steps, exps = exp)
  message(exp, ": dropped ", sum(drop), " samples; flags reset from merged_lib (", paste(steps, collapse = ", "), ")")
  invisible(sum(drop))
}
