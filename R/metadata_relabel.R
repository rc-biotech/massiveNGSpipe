#' Refresh the sample labels of a processed study after metadata fixes
#'
#' Sample labels (ORFik \code{bamVarName()}, e.g. RFP_IFNa_f4h_r1) come from the metadata via
#' \code{cleanup_metadata_for_exp()} (stage = CELL_LINE_TISSUE, condition = CONDITION_GENE,
#' fraction = FRACTION_TIMEPOINT_BATCH (+ INHIBITOR), rep = REPLICATE). They are stored in the experiment CSV and
#' in label-keyed products: study count tables (\code{QC_STATS/countTable_<region>.qs}, column names and colData)
#' and the frame QC table (\code{QC_STATS/Ribo_frames_all.csv} / \code{Ribo_frames_badzero.csv}), for both mapper
#' branches. Library files, covRle, bigWig, megafst and RiboCrypt metadata are keyed by Run and unaffected.\cr
#' This recomputes the labels from the curated metadata exactly as \code{pipeline_create_experiment()} does and,
#' if they changed, rewrites the experiment CSV in place (same files and Runs, header kept) and renames (does not
#' recount) the labels in those products, mapping old to new label by Run. Every change is backed up and logged
#' to \code{<config$project>/relabel_log.csv}. Collections pick up the new labels on the next
#' \code{pipeline_collection_org()} run.
#' @param exp experiment name (e.g. "PRJNA1393751-homo_sapiens")
#' @param config a pipeline config (its \code{complete_metadata} is the label source)
#' @param metadata data.table of curated metadata, default \code{fread(config$complete_metadata)}
#' @param exp_dir directory of experiment CSVs, default \code{ORFik::config()["exp"]}
#' @param dry_run logical, default TRUE: only report the label changes
#' @return data.table of label changes (Run, old, new), invisibly when applied
#' @export
refresh_study_labels <- function(exp, config, metadata = fread(config$complete_metadata),
                                 exp_dir = ORFik::config()["exp"], dry_run = TRUE) {
  exp_file <- file.path(exp_dir, paste0(exp, ".csv"))
  if (!file.exists(exp_file)) stop("No experiment CSV: ", exp_file)
  df <- read.experiment(exp_file, validate = FALSE, output.env = new.env())
  runs <- ORFik::runIDs(df)
  study <- as.data.table(metadata)[Run %in% runs]
  study <- study[ScientificName == organism(df)]
  if (nrow(study) != length(runs) || anyDuplicated(study$Run))
    stop(exp, ": curated metadata has ", nrow(study), " rows for ", length(runs), " experiment runs")
  study <- study[match(runs, Run)]
  lab <- cleanup_metadata_for_exp(study)
  new <- data.table(stage = as.character(lab$stage), rep = as.character(study$REPLICATE),
                    condition = as.character(lab$condition), fraction = as.character(lab$fraction))
  norm <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
  old <- data.table(stage = norm(df$stage), rep = norm(df$rep), condition = norm(df$condition), fraction = norm(df$fraction))
  for (cl in names(new)) set(new, j = cl, value = norm(new[[cl]]))
  if (identical(old, new)) {
    message(exp, ": labels unchanged")
    return(invisible(data.table(Run = character(), old = character(), new = character())))
  }
  df_new <- df
  # experiment `$<-` fails (S4 replacement); `[, col] <-` works. rep stays integer like ORFik's reader makes it.
  for (cl in c("stage", "condition", "fraction")) df_new[, cl] <- new[[cl]]
  df_new[, "rep"] <- type.convert(new$rep, as.is = TRUE)
  old_lab <- bamVarName(df); new_lab <- bamVarName(df_new)
  if (anyDuplicated(new_lab)) stop(exp, ": new labels are not unique: ", paste(unique(new_lab[duplicated(new_lab)]), collapse = ", "))
  changes <- data.table(Run = runs, old = old_lab, new = new_lab,
                        old_stage = old$stage, new_stage = new$stage, old_condition = old$condition, new_condition = new$condition,
                        old_fraction = old$fraction, new_fraction = new$fraction, old_rep = old$rep, new_rep = new$rep)
  changes <- changes[old_stage != new_stage | old_condition != new_condition | old_fraction != new_fraction | old_rep != new_rep]
  message(exp, ": ", nrow(changes), " sample labels change")
  if (dry_run) return(changes)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_dir <- file.path(config$project, "relabel_backups", paste0(exp, "_", stamp))
  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
  # 1) experiment CSV: rewrite the four label columns of the table part, keep the header lines
  file.copy(exp_file, backup_dir)
  relabel_experiment_csv(exp_file, runs, new)
  # 2) label-keyed products of both mapper branches
  short_old <- bamVarName(df, skip.libtype = TRUE); short_new <- bamVarName(df_new, skip.libtype = TRUE)
  sample_new <- bamVarName(df_new, skip.replicate = TRUE)
  for (unique_mode in c(FALSE, TRUE)) {
    d <- df; uniqueMappers(d) <- unique_mode
    q <- QCfolder(d)
    if (!dir.exists(q)) next
    bdir <- file.path(backup_dir, if (unique_mode) "unique_mappers" else "all_mappers")
    dir.create(bdir, showWarnings = FALSE)
    for (f in list.files(q, "^countTable_.*\\.qs$", full.names = TRUE)) {
      file.copy(f, bdir)
      se <- ORFik:::read_RDSQS(f)
      i <- match(colnames(se), old_lab)
      if (anyNA(i)) { warning("Unexpected columns in ", f, "; skipped"); next }
      colnames(se) <- new_lab[i]
      cd <- SummarizedExperiment::colData(se)
      if ("SAMPLE" %in% names(cd)) cd$SAMPLE <- sample_new[i]
      for (cl in intersect(c("stage", "condition", "fraction"), names(cd))) cd[[cl]] <- new[[cl]][i]
      if ("replicate" %in% names(cd)) cd$replicate <- type.convert(new$rep[i], as.is = TRUE)
      rownames(cd) <- new_lab[i]
      SummarizedExperiment::colData(se) <- cd
      ORFik:::save_RDSQS(se, f)
    }
    for (f in file.path(q, c("Ribo_frames_all.csv", "Ribo_frames_badzero.csv"))) {
      if (!file.exists(f)) next
      file.copy(f, bdir)
      x <- fread(f)
      i <- match(x$fraction, short_old)
      x[!is.na(i), fraction := short_new[i[!is.na(i)]]]
      fwrite(x, f)
    }
  }
  # 3) log
  log_file <- file.path(config$project, "relabel_log.csv")
  entry <- cbind(data.table(time = format(Sys.time(), tz = "UTC", usetz = TRUE), experiment = exp), changes)
  fwrite(entry, log_file, append = file.exists(log_file))
  stopifnot(identical(bamVarName(read.experiment(exp_file, validate = FALSE, output.env = new.env())), new_lab))
  message(exp, ": relabelled ", nrow(changes), " samples (backup ", backup_dir, ")")
  invisible(changes)
}

#' Rewrite stage/rep/condition/fraction of an ORFik experiment CSV by Run, keeping the header lines
#' @noRd
relabel_experiment_csv <- function(exp_file, runs, new) {
  lines <- readLines(exp_file)
  hdr <- grep('^"?libtype"?,', lines)[1]
  if (is.na(hdr)) stop("Cannot find the table header in ", exp_file)
  tab <- fread(text = lines[hdr:length(lines)], colClasses = "character")
  i <- match(runs, tab$Run)
  if (anyNA(i)) stop("Experiment CSV lacks runs: ", paste(runs[is.na(i)], collapse = ", "))
  for (cl in c("stage", "rep", "condition", "fraction")) set(tab, i = i, j = cl, value = new[[cl]])
  tmp <- paste0(exp_file, ".tmp")
  writeLines(lines[seq_len(hdr - 1L)], tmp)
  fwrite(tab, tmp, append = TRUE, col.names = TRUE, quote = TRUE)
  if (!file.rename(tmp, exp_file)) stop("Could not replace ", exp_file)
  invisible(exp_file)
}
