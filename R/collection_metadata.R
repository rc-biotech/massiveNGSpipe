#' Update the done-samples metadata of a pipeline project
#'
#' Writes \code{<config$project>/metadata_done_samples.csv}: the curated metadata
#' (\code{config$complete_metadata}) of every sample in every done experiment of the project,
#' in the order of the experiments. It is always computed from ALL done experiments of the
#' project, never from the subset a caller is processing, so running a subset of organisms
#' cannot shrink it.\cr
#' Safety: the file is written atomically (temporary file + rename), the previous version is
#' kept in \code{<config$project>/metadata_backups/}, and an update that would remove more than
#' \code{max_removed} samples stops unless \code{force = TRUE} (removed runs are listed).
#' Afterwards every run of the project's \code{all_samples-*} collections must be present
#' (RiboCrypt's megabrowser needs metadata for every megafst column).
#' @inheritParams pipeline_merge_org
#' @param done_experiments character, all done experiments of the project, default:
#'  \code{step_is_done_pipelines(config, "pcounts", pipelines)}
#' @param max_removed integer, default 0. Number of previously listed samples allowed to
#'  disappear without \code{force}.
#' @param force logical, default FALSE. Write even if more than \code{max_removed} samples are removed.
#' @param check_collections logical, default TRUE. Stop if a collection run is missing.
#' @param BPPARAM parallel parameters for reading experiments
#' @return the written data.table, invisibly
#' @export
update_done_samples_metadata <- function(config, pipelines = pipeline_init_all(config, gene_symbols = FALSE,
                                                                                only_complete_genomes = TRUE,
                                                                                show_status_per_exp = FALSE,
                                                                                verbose_load_annotation = FALSE),
                                         done_experiments = step_is_done_pipelines(config, "pcounts", pipelines),
                                         max_removed = 0, force = FALSE, check_collections = TRUE,
                                         BPPARAM = SerialParam()) {
  stopifnot(length(done_experiments) > 0)
  path <- done_samples_metadata_path(config)
  message("- Done samples metadata (", length(done_experiments), " done experiments): ", path)
  df_list <- bplapply(unname(done_experiments), function(e) read.experiment(e, validate = FALSE),
                      BPPARAM = BPPARAM)
  runs <- unlist(lapply(df_list, ORFik::runIDs), use.names = FALSE)
  libtype <- unique(unlist(lapply(df_list, libraryTypes)))
  if (length(libtype) != 1) stop("Only single libtype projects supported, got: ", paste(libtype, collapse = ", "))

  dt <- fread(config$complete_metadata)
  dt_complete <- dt[chmatch(runs, Run)]
  dt_complete <- dt_complete[!is.na(Run)]
  missing_meta <- setdiff(runs, dt$Run)
  if (length(missing_meta))
    warning(length(missing_meta), " done samples have no curated metadata: ", paste(head(missing_meta, 20), collapse = ", "))
  if (libtype == "RFP") {
    files <- unlist(lapply(df_list, function(df) filepath(df, "pshifted", base_folders = libFolder(df, mode = "all"))),
                    use.names = FALSE)
    dt_complete[, is_pshifted := grepl("_pshifted\\.ofst$", files[chmatch(Run, runs)])]
    if (any(!dt_complete$is_pshifted)) {
      warning("Some samples are not pshifted even though they are flagged as 'complete', they are listed below:")
      print(dt_complete[is_pshifted == FALSE, .(Run, BioProject, ScientificName)])
    } else message("All samples had existing pshifted files")
  }

  if (file.exists(path)) {
    old <- fread(path, select = "Run")
    removed <- setdiff(old$Run, dt_complete$Run)
    if (length(removed)) {
      message(length(removed), " previously listed samples are no longer done: ",
              paste(head(removed, 50), collapse = ", "), if (length(removed) > 50) ", ...")
      if (length(removed) > max_removed && !force)
        stop("Refusing to remove ", length(removed), " samples from ", path,
             " (max_removed = ", max_removed, "). Check the flags, or use force = TRUE.")
    }
    message(length(setdiff(dt_complete$Run, old$Run)), " new samples")
  }
  write_csv_atomic(dt_complete, path, backup_dir = file.path(config$project, "metadata_backups"))
  message("Saved all complete studies metadata to: ", path)
  if (check_collections)
    check_collection_runs_in_metadata(config, dt_complete$Run, libtype, unique(names(done_experiments)))
  return(invisible(dt_complete))
}

#' Extend done-samples metadata with cell line annotation
#'
#' Writes \code{<config$project>/metadata_done_samples_extended.csv} (the metadata used by the
#' RiboCrypt app): exactly the rows of \code{metadata_done_samples.csv}, plus cell line / tissue
#' annotation (Cancer_type, Cell_model, Cell_type, Organ_system, Sex, Life_stage) from the cell
#' line sheet, matched on "CELL_LINE | TISSUE". Cell_model is set to "cancer" when the sample
#' title/source mention a tumor/cancer/carcinoma/neuroblastoma/leukemia.\cr
#' Also writes \code{complete_metadata_with_cancer_status.csv} (human curated metadata with CANCER).\cr
#' The cell line sheet is only read (downloaded). Cell line/tissue pairs missing from it are
#' reported and get NA annotation; they are added to the sheet only with
#' \code{write_cell_line_sheet = TRUE} (a write to Google Sheets).
#' @param config a pipeline config
#' @param cell_line_sheet URL of the cell line annotation sheet
#' @param write_cell_line_sheet logical, default FALSE. Append new cell line/tissue pairs to the sheet.
#' @param sync_main_sheet logical, default FALSE. First download the main metadata sheet into the
#'  local curated metadata (\code{sync_sheet_safe}, read only on Google's side).
#' @param cell_lines data.table of the cell line sheet, default: downloaded from \code{cell_line_sheet}
#' @return the extended data.table, invisibly
#' @export
extend_done_samples_metadata <- function(config,
                                         cell_line_sheet = "https://docs.google.com/spreadsheets/d/1Jd7lLz6irWj_ZGMfRmaXPsVjy5-_S40n75dOydIs4hU/edit?gid=0#gid=0",
                                         write_cell_line_sheet = FALSE, sync_main_sheet = FALSE,
                                         cell_lines = NULL) {
  if (sync_main_sheet) sync_sheet_safe(config, TRUE, fix_loop = FALSE)
  done_path <- done_samples_metadata_path(config)
  if (!file.exists(done_path)) stop("Run update_done_samples_metadata() first, missing: ", done_path)
  a_live <- fread(done_path)
  a_live[, CT := paste0(CELL_LINE, " | ", TISSUE)]

  if (is.null(cell_lines)) cell_lines <- read_sheet_safe(cell_line_sheet)
  dd <- copy(as.data.table(cell_lines))
  for (col in intersect(c("CONDITION", "TIMEPOINT"), names(dd))) dd[, (col) := NULL]
  if (anyDuplicated(dd$Tissue...Cell.line))
    stop("Duplicated 'Tissue | Cell line' rows in the cell line sheet: ",
         paste(unique(dd$Tissue...Cell.line[duplicated(dd$Tissue...Cell.line)]), collapse = ", "))

  dt <- fread(config$complete_metadata)
  dt[, CT := paste0(CELL_LINE, " | ", TISSUE)]
  human <- dt[ScientificName == "Homo sapiens"]
  if (any(human$CT == " | "))
    warning(sum(human$CT == " | "), " human samples have neither CELL_LINE nor TISSUE: ",
            paste(head(human[CT == " | "]$Run, 20), collapse = ", "))
  new_pairs <- setdiff(sort(unique(human[CT != " | "]$CT)), dd$Tissue...Cell.line)
  if (length(new_pairs)) {
    message(length(new_pairs), " cell line/tissue pairs are missing in the cell line sheet: ",
            paste(new_pairs, collapse = "; "))
    if (write_cell_line_sheet) {
      upload <- rbindlist(list(dd, data.table(Tissue...Cell.line = new_pairs)), fill = TRUE)
      tmp <- tempfile(fileext = ".csv")
      fwrite(upload, tmp)
      message("Adding them to the cell line sheet (write_cell_line_sheet = TRUE)")
      write_sheet_safe(tmp, cell_line_sheet)
    } else message("Not written to the sheet (write_cell_line_sheet = FALSE); their annotation is NA.")
  }

  # human curated metadata with cancer status
  dd[, CANCER := Status..y.n. %in% "y"]
  a_res <- merge.data.table(human, dd[, .(CT = Tissue...Cell.line, CANCER)], by = "CT", all.x = TRUE)
  a_res[is.na(CANCER), CANCER := FALSE]
  a_res[grep("tumor|cancer", sample_title, ignore.case = TRUE), CANCER := TRUE]
  a_res[grep("tumor|cancer", sample_source, ignore.case = TRUE), CANCER := TRUE]
  write_csv_atomic(a_res, file.path(config$project, "complete_metadata_with_cancer_status.csv"))

  # done samples (all organisms) with cell line annotation, row order kept
  ann <- dd[, .(CT = Tissue...Cell.line, Cancer_type = Cancer.type, Cell_model = Biological.class,
                Cell_type, Organ_system, Sex, Life_stage)]
  a_live_res <- merge.data.table(a_live, ann, by = "CT", all.x = TRUE, all.y = FALSE, sort = FALSE)
  a_live_res <- a_live_res[chmatch(a_live$Run, a_live_res$Run)]
  cancer_words <- "tumor|cancer|carcinoma|neuroblastoma|leukemia"
  a_live_res[grep(cancer_words, sample_title, ignore.case = TRUE), Cell_model := "cancer"]
  a_live_res[grep(cancer_words, sample_source, ignore.case = TRUE), Cell_model := "cancer"]
  stopifnot(identical(a_live_res$Run, a_live$Run))
  path <- file.path(config$project, "metadata_done_samples_extended.csv")
  write_csv_atomic(a_live_res, path, backup_dir = file.path(config$project, "metadata_backups"))
  message("Saved extended done samples metadata (", nrow(a_live_res), " samples) to: ", path)
  return(invisible(a_live_res))
}

#' @noRd
done_samples_metadata_path <- function(config) file.path(config$project, "metadata_done_samples.csv")

#' Write a csv atomically, optionally keeping a dated backup of the previous version
#' @noRd
write_csv_atomic <- function(x, path, backup_dir = NULL) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp)
  if (!is.null(backup_dir) && file.exists(path)) {
    dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(path, file.path(backup_dir, paste0(sub("\\.csv$", "", basename(path)), "_",
                                                 format(file.mtime(path), "%Y%m%d_%H%M%S"), ".csv")),
              overwrite = FALSE, copy.date = TRUE)
  }
  if (!file.rename(tmp, path)) stop("Could not replace ", path)
  invisible(path)
}

#' Every run of the project's collections must have done-samples metadata
#' @noRd
check_collection_runs_in_metadata <- function(config, metadata_runs, libtype, organisms) {
  # the current per-organism collections, named as in pipeline_collection_org (snapshots excluded)
  wanted <- organism_collection_exp_name(organisms)
  if (libtype != "RFP") wanted <- paste0(wanted, "_", libtype)
  if (config$mode == "local") wanted <- paste0("local_", wanted)
  collections <- list.experiments(validate = FALSE, libtypeExclusive = libtype)
  collections <- collections[name %in% wanted]
  missing <- rbindlist(lapply(collections$name, function(e) {
    runs <- tryCatch(ORFik::runIDs(read.experiment(e, validate = FALSE)), error = function(err) character())
    miss <- setdiff(runs, metadata_runs)
    if (length(miss)) data.table(collection = e, Run = miss)
  }))
  if (nrow(missing)) {
    print(missing[, .(missing_runs = .N, example = Run[1]), by = collection])
    stop(nrow(missing), " collection runs lack done-samples metadata (RiboCrypt megabrowser would fail). ",
         "Rebuild those collections or check their flags.")
  }
  message("All ", nrow(collections), " collections are covered by the done-samples metadata")
  invisible(TRUE)
}
