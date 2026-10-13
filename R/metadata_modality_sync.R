#' Move metadata rows to the sheet of their library modality
#'
#' Each pipeline project (one per library modality) curates its samples in its own Google sheet
#' (and local \code{next_round_manual.csv}). Rows of other modalities that ended up in a project's
#' sheet (e.g. RNA-seq or disome rows in the Ribo-seq sheet) are moved to the sheet of the project
#' that processes that modality. Modalities without a pipeline are parked in a local csv
#' (\code{<other_project>/next_round_manual.csv}, same columns) until a pipeline exists.\cr
#' The online (Google) sheet is always the source of truth: with \code{use_google = TRUE} every read
#' is a fresh download of the sheet, and the local \code{next_round_manual.csv} is only overwritten
#' as a copy after a successful write, never used as input.\cr
#' Safety rules:
#' \itemize{
#'  \item online and local projects are separate families; rows never cross between them
#'  \item a row is only copied if its Run is not already in the destination
#'  \item if the destination already has the Run: the source row is removed only when its
#'   annotation (\code{compare_cols}) equals the destination row, otherwise it is a conflict and
#'   both stay
#'  \item a Run that appears with different LIBRARYTYPE in different sheets is a conflict, untouched
#'  \item rows are removed from a source only after a read-back confirmed them in the destination
#'  \item every sheet is backed up first; each write is a fresh read + append/remove + verified upload
#'  \item destination FINAL_LIST.csv files are not touched (the project owner validates new rows)
#' }
#' @param family "online" (NGS_pipeline*) or "local" (local_pipeline_*)
#' @param base_dir directory holding the project folders
#' @param homes named character, LIBRARYTYPE -> project folder that processes it
#' @param other_project project folder for modalities without a pipeline (local csv only)
#' @param compare_cols annotation columns that must agree for a duplicate to be dropped
#' @param on_conflict "keep" (default: conflicts stay where they are) or "archive": the home sheet's
#'  row wins; the stray copy is appended to \code{<other_project>/conflict_archive.csv} (with the reason)
#'  and removed from its sheet, but only if the Run remains in a sheet that is the home of its label there
#' @param dry_run logical, default TRUE: only return the plan
#' @param backup_dir directory for sheet backups and the plan
#' @param use_google logical, default TRUE: read/write the projects' Google sheets (else the local csv)
#' @return the plan (data.table), invisibly when applied
#' @export
sync_modalities <- function(family = c("online", "local"), base_dir = "~/livemount/Bio_data",
                            homes = modality_homes(family), other_project = modality_other_project(family),
                            compare_cols = c("LIBRARYTYPE", "BATCH", "REPLICATE", "TIMEPOINT", "TISSUE", "CELL_LINE",
                                             "CONDITION", "GENE", "FRACTION", "INHIBITOR"),
                            on_conflict = c("keep", "archive"),
                            dry_run = TRUE, backup_dir = tempfile("modality_sync_"), use_google = TRUE) {
  family <- match.arg(family); on_conflict <- match.arg(on_conflict)
  projects <- unique(unname(homes))
  stopifnot(length(projects) > 0, !other_project %in% projects)
  dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)
  sheets <- lapply(setNames(projects, projects), function(p) modality_read(p, base_dir, use_google))
  other_path <- file.path(base_dir, other_project, "next_round_manual.csv")
  other <- if (file.exists(other_path)) fread(other_path) else data.table()
  for (p in projects) fwrite(sheets[[p]], file.path(backup_dir, paste0(p, "_before.csv")))
  if (nrow(other)) fwrite(other, file.path(backup_dir, paste0(other_project, "_before.csv")))

  plan <- modality_plan(sheets, homes, other, other_project, compare_cols)
  if (on_conflict == "archive") plan <- modality_archive_eligible(plan, sheets, homes)
  fwrite(plan, file.path(backup_dir, "plan.csv"))
  modality_plan_summary(plan)
  if (dry_run) {
    message("Dry run: nothing written. Plan and backups: ", backup_dir)
    return(plan)
  }

  # 1) destinations: fresh read, append, verify, write
  moves <- plan[action %in% c("move", "park")]
  for (dest in setdiff(unique(moves$to), other_project)) {
    rows <- modality_rows(sheets, moves[to == dest])
    modality_append(dest, rows, base_dir, use_google)
  }
  if (nrow(moves[to == other_project])) {
    rows <- modality_rows(sheets, moves[to == other_project])
    dir.create(dirname(other_path), showWarnings = FALSE, recursive = TRUE)
    updated <- rbindlist(list(other, rows[!Run %in% other$Run]), fill = TRUE)
    write_csv_atomic(updated, other_path)
  }
  archive_path <- file.path(base_dir, other_project, "conflict_archive.csv")
  arch <- plan[action == "archive"]
  if (nrow(arch)) {
    rows <- modality_rows(sheets, arch)
    rows[, `:=`(archived_from = arch$from, archive_reason = arch$reason, archived_at = format(Sys.time(), tz = "UTC", usetz = TRUE))]
    old <- if (file.exists(archive_path)) fread(archive_path) else data.table()
    dir.create(dirname(archive_path), showWarnings = FALSE, recursive = TRUE)
    write_csv_atomic(rbindlist(list(old, rows), fill = TRUE), archive_path)
  }
  # 2) sources: remove only rows confirmed in their destination (or in the archive)
  remove <- plan[action %in% c("move", "park", "drop_duplicate")]
  confirmed <- modality_confirm(remove, base_dir, use_google, other_path, other_project)
  if (nrow(arch)) {
    archived <- fread(archive_path)
    confirmed <- c(confirmed, arch[paste(Run, from) %in% paste(archived$Run, archived$archived_from)]$Run)
    remove <- rbind(remove, arch)
  }
  for (src in unique(remove$from)) {
    runs <- remove[from == src & Run %in% confirmed]$Run
    if (length(runs)) modality_remove(src, runs, base_dir, use_google)
  }
  plan[, applied := action %in% c("move", "park", "drop_duplicate", "archive") & Run %in% confirmed]
  fwrite(plan, file.path(backup_dir, "plan_applied.csv"))
  message("Applied: ", sum(plan$applied), " rows; plan/backups in ", backup_dir)
  invisible(plan)
}

#' @rdname sync_modalities
#' @export
modality_homes <- function(family = c("online", "local")) {
  family <- match.arg(family)
  if (family == "online") {
    c(RFP = "NGS_pipeline", RNA = "NGS_pipeline_RNA", disome = "NGS_pipeline_disome", SSU = "NGS_pipeline_SSU")
  } else {
    c(RFP = "local_pipeline_ribo", RNA = "local_pipeline_RNA", disome = "local_pipeline_disome", SSU = "local_pipeline_SSU")
  }
}

#' @rdname sync_modalities
#' @export
modality_other_project <- function(family = c("online", "local")) {
  if (match.arg(family) == "online") "NGS_pipeline_other" else "local_pipeline_other"
}

#' The move plan: one row per foreign-modality sample (plus cross-sheet conflicts)
#' @noRd
modality_plan <- function(sheets, homes, other, other_project, compare_cols) {
  norm <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; trimws(x) }
  all_rows <- rbindlist(lapply(names(sheets), function(p) {
    s <- sheets[[p]]
    if (!nrow(s)) return(NULL)
    data.table(Run = as.character(s$Run), project = p, LIBRARYTYPE = norm(s$LIBRARYTYPE), row = seq_len(nrow(s)))
  }))
  # runs labelled with different modalities in different sheets: never touched
  labels <- all_rows[, .(n_types = uniqueN(LIBRARYTYPE)), by = Run]
  multi_type <- labels[n_types > 1]$Run
  plan <- rbindlist(lapply(names(sheets), function(p) {
    s <- sheets[[p]]
    if (!nrow(s)) return(NULL)
    lt <- norm(s$LIBRARYTYPE)
    home <- unname(homes[lt])
    foreign <- which(is.na(home) | home != p)
    if (!length(foreign)) return(NULL)
    rbindlist(lapply(foreign, function(i) {
      run <- as.character(s$Run[i]); to <- if (is.na(home[i])) other_project else home[i]
      if (run %in% multi_type) {
        return(data.table(Run = run, from = p, to = to, LIBRARYTYPE = lt[i], action = "conflict",
                          reason = "Run has different LIBRARYTYPE in different sheets"))
      }
      dest <- if (to == other_project) other else sheets[[to]]
      hit <- if (nrow(dest)) which(as.character(dest$Run) == run) else integer()
      if (!length(hit)) {
        return(data.table(Run = run, from = p, to = to, LIBRARYTYPE = lt[i],
                          action = if (to == other_project) "park" else "move", reason = ""))
      }
      cols <- intersect(compare_cols, intersect(names(s), names(dest)))
      differs <- cols[vapply(cols, function(cl) norm(s[[cl]][i]) != norm(dest[[cl]][hit[1]]), logical(1))]
      if (length(differs)) {
        data.table(Run = run, from = p, to = to, LIBRARYTYPE = lt[i], action = "conflict",
                   reason = paste("already in destination with different", paste(differs, collapse = ",")))
      } else {
        data.table(Run = run, from = p, to = to, LIBRARYTYPE = lt[i], action = "drop_duplicate",
                   reason = "identical annotation already in destination")
      }
    }))
  }))
  if (is.null(plan) || !nrow(plan))
    plan <- data.table(Run = character(), from = character(), to = character(), LIBRARYTYPE = character(),
                       action = character(), reason = character())
  plan[]
}

#' Conflicts become "archive" only if the Run stays in a sheet that is the home of its label there
#' @noRd
modality_archive_eligible <- function(plan, sheets, homes) {
  norm <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; trimws(x) }
  at_home <- unique(unlist(lapply(names(sheets), function(p) {
    s <- sheets[[p]]
    if (!nrow(s)) return(NULL)
    h <- unname(homes[norm(s$LIBRARYTYPE)])
    as.character(s$Run)[!is.na(h) & h == p]
  })))
  plan[action == "conflict" & Run %in% at_home, `:=`(action = "archive", reason = paste("home sheet wins;", reason))]
  plan[action == "conflict", reason := paste(reason, "(kept: Run has no home-sheet copy)")]
  plan[]
}

#' @noRd
modality_plan_summary <- function(plan) {
  if (!nrow(plan)) { message("All sheets already hold only their own modality."); return(invisible()) }
  print(plan[, .N, by = .(from, to, action)][order(from, to, action)])
  invisible()
}

#' @noRd
modality_url <- function(project, base_dir) {
  f <- file.path(base_dir, project, "google_api_url_table.csv")
  if (!file.exists(f)) stop("No google_api_url_table.csv for ", project)
  fread(f)$id[1]
}

#' @noRd
modality_read <- function(project, base_dir, use_google) {
  if (use_google) return(read_sheet_safe(modality_url(project, base_dir)))
  f <- file.path(base_dir, project, "next_round_manual.csv")
  if (file.exists(f)) fread(f) else data.table()
}

#' @noRd
modality_write <- function(project, dt, base_dir, use_google) {
  local <- file.path(base_dir, project, "next_round_manual.csv")
  write_csv_atomic(dt, local, backup_dir = file.path(base_dir, project, "metadata_backups"))
  if (use_google) write_sheet_safe(local, modality_url(project, base_dir), sheet = 1)
  invisible(dt)
}

#' Source rows for planned runs, in plan order
#' @noRd
modality_rows <- function(sheets, plan_part) {
  rbindlist(lapply(seq_len(nrow(plan_part)), function(i) {
    s <- sheets[[plan_part$from[i]]]
    s[as.character(Run) == plan_part$Run[i]][1]
  }), fill = TRUE)
}

#' Same column types as the destination, then append; existing rows must stay unchanged
#' @noRd
modality_append <- function(project, rows, base_dir, use_google) {
  current <- modality_read(project, base_dir, use_google)            # fresh read right before writing
  rows <- rows[!as.character(Run) %in% as.character(current$Run)]
  if (!nrow(rows)) return(invisible(current))
  extra <- setdiff(names(rows), names(current))
  if (nrow(current) && length(extra)) rows[, (extra) := NULL]
  for (cl in intersect(names(rows), names(current))) {
    if (!identical(class(rows[[cl]]), class(current[[cl]])))
      set(rows, j = cl, value = methods::as(rows[[cl]], class(current[[cl]])[1]))
  }
  updated <- rbindlist(list(current, rows), fill = TRUE)
  # existing rows must keep their values (compared as text: a type can widen for the new rows)
  as_txt <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
  changed <- names(current)[vapply(names(current), function(cl)
    !identical(as_txt(updated[[cl]][seq_len(nrow(current))]), as_txt(current[[cl]])), logical(1))]
  if (nrow(updated) != nrow(current) + nrow(rows) || length(changed))
    stop("Appending to ", project, " would change existing rows (", paste(changed, collapse = ", "),
         "); nothing written")
  modality_write(project, updated, base_dir, use_google)
  back <- modality_read(project, base_dir, use_google)
  if (nrow(back) != nrow(updated) || !all(as.character(rows$Run) %in% as.character(back$Run)))
    stop("Read-back of ", project, " does not contain the appended rows; sources left unchanged")
  message("- ", project, ": appended ", nrow(rows), " rows (", nrow(current), " -> ", nrow(back), ")")
  invisible(back)
}

#' Runs confirmed present in their destination (re-read now)
#' @noRd
modality_confirm <- function(remove, base_dir, use_google, other_path, other_project) {
  dests <- unique(remove$to)
  present <- unlist(lapply(dests, function(d) {
    runs <- if (d == other_project) { if (file.exists(other_path)) fread(other_path)$Run else character() } else
      modality_read(d, base_dir, use_google)$Run
    remove[to == d & Run %in% as.character(runs)]$Run
  }))
  unique(present)
}

#' @noRd
modality_remove <- function(project, runs, base_dir, use_google) {
  current <- modality_read(project, base_dir, use_google)            # fresh read right before writing
  keep <- !as.character(current$Run) %in% runs
  updated <- current[keep]
  modality_write(project, updated, base_dir, use_google)
  back <- modality_read(project, base_dir, use_google)
  if (nrow(back) != nrow(updated) || any(as.character(back$Run) %in% runs))
    stop("Read-back of ", project, " after removal is not as expected; check backups")
  message("- ", project, ": removed ", sum(!keep), " moved rows (", nrow(current), " -> ", nrow(back), ")")
  invisible(back)
}
