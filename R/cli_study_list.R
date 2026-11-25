# Study-list runner behind `mngsp plan|run|status` (ported from the Valen runner run_study_list_v5.R, 2026-11-24):
# plan per experiment, --reset, manual trim tables, detector provenance, run_pipeline().

#' dbGaP controlled-access studies
#'
#' BioProjects whose runs need dbGaP authorisation (downloads fail with 403). Skipped by
#' \code{mngsp run} unless \code{--no-skip-dbgap}.
#' @return character vector of BioProject accessions
#' @export
dbgap_studies <- function() c("PRJNA1028780", "PRJNA1158872", "PRJNA401960", "PRJNA429519", "PRJNA615003")

#' Read a study list (one accession per line, '#' lines and blanks ignored)
#' @param file path
#' @return unique character vector
#' @noRd
read_study_list <- function(file) {
  x <- trimws(readLines(file, warn = FALSE))
  unique(x[nzchar(x) & !startsWith(x, "#")])
}

#' Plan a study list: one row per experiment (study x organism)
#'
#' Read-only. Studies are matched by BioProject or study_accession of the FINAL list.
#' @param study_list path of the study list, or a character vector of accessions
#' @param config the mNGSp config
#' @param keep character, studies to leave untouched
#' @param trim_table data.table manual trim table (status "hold" studies are left out) or NULL
#' @param organisms data.table with BioProject, ScientificName (only these experiments) or NULL
#' @param exclude_runs character, runs never processed
#' @param exclude_studies character, studies never processed
#' @param skip_dbgap logical, leave out \code{dbgap_studies()}
#' @return list(plan, meta_rows, accessions, todo, missing, counts)
#' @export
study_list_plan <- function(study_list, config, keep = character(), trim_table = NULL, organisms = NULL,
                            exclude_runs = character(), exclude_studies = character(), skip_dbgap = TRUE) {
  acc <- if (length(study_list) == 1 && file.exists(study_list)) read_study_list(study_list) else unique(study_list)
  dbgap <- if (skip_dbgap) dbgap_studies() else character()
  hold <- if (!is.null(trim_table) && nrow(trim_table)) trim_table[status == "hold"]$BioProject else character()
  todo <- setdiff(acc, c(keep, dbgap, hold, exclude_studies))
  if (!file.exists(config$complete_metadata))
    stop("No FINAL list at ", config$complete_metadata, ": curate the studies first (mngsp curate <study list>)", call. = FALSE)
  final <- fread(config$complete_metadata)
  meta_rows <- final[(BioProject %in% todo | study_accession %in% todo) & !(Run %in% exclude_runs)]
  if (!is.null(organisms))
    meta_rows <- meta_rows[paste(BioProject, ScientificName) %in% paste(organisms$BioProject, organisms$ScientificName)]
  missing <- setdiff(todo, c(final$BioProject, final$study_accession))
  ref <- config$config["ref"]; bam <- config$config["bam"]
  # experiments are named by study_accession (a GEO series for many studies), not by BioProject: name the plan the same way,
  # otherwise --reset clears a non-existing <BioProject>-<organism> and the real <GSE>-<organism> keeps its flags and folder
  plan <- meta_rows[, .(runs = .N), by = .(study = BioProject, accession = study_accession, organism = ScientificName)]
  steps <- names(config$flag)
  if (nrow(plan)) {
    plan[, experiment := paste(accession, reference_folder_name(organism), sep = "-")]
    rf <- file.path(ref, reference_folder_name(plan$organism))
    plan[, reference_complete := file.exists(file.path(rf, "outputs.rds")) &
                                 file.exists(file.path(rf, "STAR_index", "genomeDir", "SA"))]
    plan[, folder_exists := dir.exists(file.path(bam, experiment))]
    done <- vapply(plan$experiment, function(e) vapply(steps, function(s) step_is_done(config, s, e), logical(1)),
                   logical(length(steps)))
    done <- matrix(done, nrow = length(steps))
    plan[, flags_set := as.integer(colSums(done))]
    plan[, flags_total := length(steps)]
    plan[, next_step := apply(done, 2, function(d) if (all(d)) "done" else steps[which(!d)[1]])]
    plan[, manual_trim := if (!is.null(trim_table) && nrow(trim_table))
                            study %in% trim_table[status == "apply"]$BioProject else FALSE]
  } else {
    plan <- data.table(study = character(), accession = character(), organism = character(), runs = integer(), experiment = character(),
                       reference_complete = logical(), folder_exists = logical(), flags_set = integer(),
                       flags_total = integer(), next_step = character(), manual_trim = logical())
  }
  counts <- c(study_list = length(acc), keep = sum(acc %in% keep), dbgap = sum(acc %in% dbgap),
              hold = sum(acc %in% hold), excluded = sum(acc %in% exclude_studies), missing = length(missing),
              experiments = nrow(plan), reference_missing = sum(!plan$reference_complete))
  list(plan = plan[], meta_rows = meta_rows, accessions = acc, todo = todo, missing = missing, counts = counts)
}

#' Reset experiments for a full restart
#'
#' Clears all flags and per-sample markers, error records of the last session, the study folder (processed data)
#' and the experiment CSV. Raw FASTQ is downloaded again for online studies. A production action: only
#' \code{mngsp run --apply --reset} calls it.
#' @param config the mNGSp config
#' @param experiments character, experiment names (<study>-<organism folder>)
#' @return invisible experiments
#' @export
study_list_reset <- function(config, experiments) {
  if (!length(experiments)) return(invisible(experiments))
  remove_flag_all_exp(config, steps = "all", exps = experiments, samples = TRUE)
  for (e in experiments) {
    try(suppressMessages(clear_pipeline_error(config, e, note = "mngsp run --reset")), silent = TRUE)
    d <- file.path(config$config["bam"], e)
    if (dir.exists(d)) unlink(d, recursive = TRUE)
    unlink(file.path(config$config["exp"], paste0(e, ".csv")))
  }
  message("Reset ", length(experiments), " experiments")
  invisible(experiments)
}

cli_read_inputs <- function(v) {
  list(keep = if (!is.null(v$keep)) read_study_list(v$keep) else character(),
       trim = if (!is.null(v$trim)) {
         mt <- fread(v$trim, colClasses = c(adapter = "character"))
         miss <- setdiff(c("BioProject", "adapter", "barcode5p_size", "barcode3p_size", "status"), names(mt))
         if (length(miss)) stop("--trim table ", v$trim, " lacks columns: ", paste(miss, collapse = ", "), call. = FALSE)
         mt
       } else NULL,
       organisms = if (!is.null(v$organisms)) {
         o <- fread(v$organisms)
         if (!all(c("BioProject", "ScientificName") %in% names(o)))
           stop("--organisms ", v$organisms, " needs columns BioProject,ScientificName", call. = FALSE)
         o
       } else NULL,
       exclude_runs = if (!is.null(v$exclude_runs)) {
         x <- fread(v$exclude_runs)
         if (!("Run" %in% names(x))) stop("--exclude-runs ", v$exclude_runs, " needs a Run column", call. = FALSE)
         x$Run
       } else character())
}

cli_plan_from_values <- function(r, config) {
  v <- r$values; inp <- cli_read_inputs(v)
  p <- study_list_plan(r$study_list, config, keep = inp$keep, trim_table = inp$trim, organisms = inp$organisms,
                       exclude_runs = inp$exclude_runs, exclude_studies = v$exclude_studies %||% character(),
                       skip_dbgap = !isFALSE(v$skip_dbgap))
  p$inputs <- inp
  p
}

cli_open_file_limit <- function() {
  x <- tryCatch(system("ulimit -n", intern = TRUE), error = function(e) NA_character_)[1]
  if (identical(x, "unlimited")) Inf else suppressWarnings(as.numeric(x))
}

cli_is_footprint <- function(config) libtype_long_to_short(config$preset) %in% footprint_library_types()

cli_list_tag <- function(list_file) tools::file_path_sans_ext(basename(list_file))

#' mngsp plan / run
#' @noRd
cli_cmd_run <- function(r, apply) {
  v <- r$values
  config <- cli_build_config(v)
  preflight <- NULL
  if (isTRUE(v$preflight) && cli_is_footprint(config))
    preflight <- preflight_detector(checkout = v$dev %||% (if (nzchar(Sys.getenv("MNGSP_DEV"))) Sys.getenv("MNGSP_DEV") else NA_character_),
                                    min_orfik = v$min_orfik)
  nofile <- cli_open_file_limit()
  if (!is.na(nofile) && nofile < v$min_open_files) {
    msg <- paste0("Open-file limit is ", nofile, " (< ", v$min_open_files, "); STAR needs more: run 'ulimit -n 524288' ",
                  "in the shell (the mngsp wrapper raises it up to the hard limit; see 'ulimit -Hn')")
    if (apply) stop(msg, call. = FALSE) else message("Note: ", msg)
  }
  p <- cli_plan_from_values(r, config)
  out_dir <- v$out_dir %||% config$project
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tag <- cli_list_tag(r$study_list)
  plan_file <- file.path(out_dir, paste0(v$file_prefix, "plan_", tag, ".csv"))
  fwrite(p$plan, plan_file)
  print(p$plan[, .(experiment, runs, reference = reference_complete, folder = folder_exists,
                   flags = paste0(flags_set, "/", flags_total), next_step, manual_trim)], nrows = 1000)
  k <- p$counts
  message("study list: ", k["study_list"], " | keep: ", k["keep"], " | dbGaP: ", k["dbgap"], " | hold (trim table): ",
          k["hold"], " | excluded: ", k["excluded"], " | not in FINAL list: ", k["missing"], " | experiments: ",
          k["experiments"], " (reference missing: ", k["reference_missing"], ") | mode: ",
          if (isTRUE(v$reset)) "full reset" else "resume", " | project: ", config$project)
  if (length(p$missing)) message("Not in the FINAL list (", v$metadata_hint, "): ", paste(p$missing, collapse = ", "))
  if (any(!p$plan$reference_complete)) {
    folders <- paste(unique(reference_folder_name(p$plan[reference_complete == FALSE]$organism)), collapse = " ")
    message("Missing references: ", gsub("{folders}", folders, v$reference_hint, fixed = TRUE),
            if (isTRUE(v$download_references)) " (--download-references: pipeline_init_all() fetches them)")
  }
  if (!apply) {
    message("Dry run: nothing changed", if (isTRUE(v$reset)) " (--reset needs --apply)", ". Plan: ", plan_file,
            ". Add --apply to run.")
    return(invisible(0L))
  }
  if (isTRUE(v$reset)) {
    stopifnot(!any(p$plan$study %in% p$inputs$keep))
    study_list_reset(config, unique(p$plan$experiment))
  }
  study_list_execute(p, config, v, r$study_list, preflight, out_dir)
}

#' Run a planned study list (after an optional reset)
#' @noRd
study_list_execute <- function(p, config, v, list_file, preflight, out_dir) {
  tag <- cli_list_tag(list_file); mt <- p$inputs$trim; meta_rows <- p$meta_rows; todo <- p$todo
  if (!nrow(meta_rows)) stop("Nothing to run: no FINAL list rows for the study list", call. = FALSE)
  meta <- file.path(out_dir, paste0("subset_", tag, v$subset_suffix, ".csv"))
  fwrite(meta_rows, meta)
  only_complete <- !isTRUE(v$download_references)
  pipelines <- pipeline_init_all(config, complete_metadata = meta, only_complete_genomes = only_complete,
                                 only_complete_indices = only_complete, gene_symbols = FALSE,
                                 show_status_per_exp = FALSE)
  # pipelines are named by study_accession (often a GSE), the list holds BioProjects: map between them
  acc_map <- unique(meta_rows[, .(study_accession, BioProject)])
  run <- names(pipelines)[names(pipelines) %in% c(todo, acc_map[BioProject %in% todo]$study_accession)]
  bp_of <- function(s) { b <- acc_map[match(s, study_accession)]$BioProject; ifelse(is.na(b), s, b) }
  # 2026-11-25: studies (BioProjects of the list) vs experiments (one per study x organism) counted separately (was "73 of 72")
  message("Studies with a complete reference: ", length(unique(bp_of(run))), " of ", length(todo),
          " (", length(run), " experiments)")
  if (!length(run)) stop("No study has a complete reference (", v$reference_hint, ")", call. = FALSE)
  confs <- function(s) lapply(pipelines[[s]]$organisms, function(o) o$conf)
  trim_dir <- function(conf) { d <- file.path(conf["bam"], "trim"); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
  # manual trim tables (lost with a reset): copied before trimming, they take precedence over detection
  if (!is.null(v$manual_dir)) for (s in run) for (conf in confs(s)) {
    src <- file.path(v$manual_dir, conf["exp"])
    if (!dir.exists(src)) next
    file.copy(list.files(src, full.names = TRUE), trim_dir(conf), overwrite = TRUE)
    message("Manual tables ", conf["exp"], ": ", paste(list.files(src), collapse = ", "))
  }
  gaps <- unlist(lapply(run, function(s) lapply(names(pipelines[[s]]$organisms), function(o) {
    conf <- pipelines[[s]]$organisms[[o]]$conf
    m <- manual_table_gaps(file.path(conf["bam"], "trim"), meta_rows[ScientificName == o & (study_accession == s | BioProject == s)]$Run)
    if (length(m)) paste0(conf["exp"], ": ", length(m), " runs not in barcodes_manual.csv") })))
  if (length(gaps)) message("barcodes_manual.csv lists only some runs (the others are auto-detected):\n  ", paste(gaps, collapse = "\n  "))
  # provenance per experiment: which code and detector trimmed it
  prov <- cbind(data.table(mngsp_version = as.character(utils::packageVersion("massiveNGSpipe")),
                           detector_mode = barcode_detector_mode(), preset = config$preset),
                if (!is.null(preflight)) as.data.table(preflight)[, !c("mngsp_version", "detector_mode"), with = FALSE])
  for (s in run) for (conf in confs(s))
    fwrite(cbind(experiment = conf["exp"], host = Sys.info()[["nodename"]], prov), file.path(trim_dir(conf), "detector_provenance.csv"))
  if (!is.null(preflight)) fwrite(preflight, file.path(out_dir, paste0(v$file_prefix, "preflight_", tag, ".csv")))
  if (!is.null(mt)) for (s in run[bp_of(run) %in% mt[status == "apply"]$BioProject]) {
    m <- mt[BioProject == bp_of(s)][1]
    for (conf in confs(s)) trim_dir(conf)
    adapters_manual_assign(pipelines[[s]], m$adapter)
    barcodes_manual_assign(pipelines[[s]], m$barcode5p_size, m$barcode3p_size)
    if ("mate" %in% names(mt) && !is.na(m$mate) && m$mate %in% c(1, 2)) {
      runs <- meta_rows[BioProject == bp_of(s) & LibraryLayout == "PAIRED"]$Run
      for (conf in confs(s))
        fwrite(data.table(Run = runs, mate = m$mate, reverse = toupper(as.character(m$reverse)) == "TRUE", adapter = m$adapter,
                          trim_front = m$barcode5p_size, trim_tail = m$barcode3p_size, source = basename(v$trim)),
               file.path(conf["bam"], "trim", "mate_override_manual.csv"))
      message("Mate override ", s, ": mate ", m$mate, ", reverse ", m$reverse, " for ", length(runs), " paired runs")
    }
    message("Manual trim ", s, ": adapter ", m$adapter, ", UMI 5' ", m$barcode5p_size, " / 3' ", m$barcode3p_size)
  }
  ok <- run_pipeline(pipelines[run], config, wait = v$wait)
  print(last_session_errors(config))
  invisible(if (isTRUE(ok)) 0L else 2L)
}

#' barcodes_manual.csv gap: runs of an experiment missing from an existing table (they are auto-detected)
#' @param trim_dir the experiment's trim folder
#' @param runs character, runs of the experiment
#' @return character, runs not in the table (empty without a table)
#' @noRd
manual_table_gaps <- function(trim_dir, runs) {
  f <- file.path(trim_dir, "barcodes_manual.csv")
  if (!file.exists(f)) return(character())
  setdiff(runs, fread(f)$Run)
}

#' mngsp status
#' @noRd
cli_cmd_status <- function(r) {
  v <- r$values
  config <- cli_build_config(v)
  if (!is.null(r$study_list)) {
    p <- cli_plan_from_values(r, config)
    errs <- tryCatch(pipeline_errors(config, index = v$session), error = function(e) NULL)
    tab <- p$plan[, .(experiment, runs, reference = reference_complete, flags = paste0(flags_set, "/", flags_total), next_step)]
    tab[, error := if (!is.null(errs) && nrow(errs)) errs$step[match(experiment, errs$exp)] else NA_character_]
    print(tab, nrows = 1000)
    message(sum(tab$next_step == "done"), " of ", nrow(tab), " experiments done; ", sum(!is.na(tab$error)),
            " with an error record (mngsp errors)",
            if (length(p$missing)) paste0("; not in the FINAL list: ", paste(p$missing, collapse = ", ")))
    return(invisible(0L))
  }
  path <- tryCatch(get_checklist_path(config, v$session), error = function(e) NA_character_)
  if (is.na(path) || !file.exists(path)) {
    message("No checklist yet for session ", v$session, " in ", config$project,
            " (no run started; 'mngsp status <study list>' shows flags per experiment)")
    return(invisible(0L))
  }
  last_session_progress(config, v$session)
  message("Checklist: ", path, " (stale after a hard kill)")
  invisible(0L)
}

#' mngsp errors
#' @noRd
cli_cmd_errors <- function(r) {
  v <- r$values
  config <- cli_build_config(v)
  e <- pipeline_errors(config, index = v$session)
  if (!nrow(e)) message("No error records in session ", v$session, " of ", config$project) else {
    for (i in seq_len(nrow(e))) {
      cat(sprintf("%s  [%s]  %s\n  %s\n", e$exp[i], e$step[i], format(e$time[i], "%Y-%m-%d %H:%M"), e$message[i]))
      if (!is.na(e$hint[i]) && nzchar(e$hint[i])) cat("  hint: ", e$hint[i], "\n", sep = "")
    }
    message(nrow(e), " error record(s). Fix the cause, then: mngsp retry <experiment> (a running pipeline retries it)")
  }
  if (v$events > 0) print(pipeline_error_events(config, n = v$events))
  invisible(0L)
}

#' mngsp retry
#' @noRd
cli_cmd_retry <- function(r) {
  v <- r$values
  config <- cli_build_config(v)
  exps <- r$positional
  if (isTRUE(v$all)) {
    if (length(exps)) stop("give experiments or --all, not both", call. = FALSE)
    exps <- pipeline_errors(config, index = v$session)$exp
    if (!length(exps)) { message("No error records in session ", v$session); return(invisible(0L)) }
  }
  if (!length(exps)) stop("name the experiments to retry (see mngsp errors), or --all", call. = FALSE)
  moved <- clear_pipeline_error(config, exps, index = v$session, note = v$note)
  message("Cleared ", length(moved), " error record(s) (moved to resolved/): ", paste(exps, collapse = ", "),
          ". A running pipeline retries them next round; else start 'mngsp run <list> --apply'.")
  invisible(0L)
}

#' mngsp curate
#' @noRd
cli_cmd_curate <- function(r) {
  v <- r$values
  if (is.null(v$google_sheet)) v$google_sheet <- "none"
  config <- cli_build_config(v)
  acc <- if (isTRUE(v$validate_only)) NULL else read_study_list(r$study_list)
  ok <- curate_metadata(acc, config, organisms = v$organism %||% "all", open_editor = FALSE, fix_loop = FALSE,
                        update_google_sheet = !is.null(config$google_url))
  if (isTRUE(ok) && file.exists(config$complete_metadata)) {
    fl <- fread(config$complete_metadata)
    acc <- read_study_list(r$study_list)
    message("FINAL list ", config$complete_metadata, ": ", fl[BioProject %in% acc | study_accession %in% acc, .N],
            " runs of the study list. Next: mngsp plan ", r$study_list)
  }
  message("If validation failed: edit ", config$temp_metadata, " (KEEP, LIBRARYTYPE, CONDITION, REPLICATE, ...) and run ",
          "'mngsp curate ", r$study_list, " --validate-only'")
  invisible(0L)
}

#' mngsp init
#' @noRd
cli_cmd_init <- function(r, cwd) {
  v <- r$values
  if (length(r$positional) != 1) stop("usage: mngsp init <dir>", call. = FALSE)
  dir <- cli_abs_path(r$positional, cwd)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  conf <- file.path(dir, "mngsp.conf"); lst <- file.path(dir, "studies.txt")
  for (f in c(conf, lst)) if (file.exists(f) && !isTRUE(v$force)) stop(f, " exists (use --force to overwrite)", call. = FALSE)
  writeLines(cli_config_template(project_dir = "project", preset = v$preset %||% "Ribo-seq", mode = v$mode %||% "online"), conf)
  writeLines(c("# One study per line: BioProject (PRJNA...), GEO series or SRA study for mode online,",
               "# folder names under the fastq folder for mode local. '#' lines are comments;",
               "# '#: key: value' lines are settings for this list only (e.g. '#: threads: 8').",
               v$example), lst)
  dir.create(file.path(dir, "project"), showWarnings = FALSE)
  message("Created ", conf, " and ", lst, " (example study ", v$example, ")\nNext:\n",
          "  cd ", dir, "\n  mngsp preflight\n  mngsp curate studies.txt --organism 'Homo sapiens'\n",
          "  mngsp plan studies.txt\n  mngsp run studies.txt --apply --download-references")
  invisible(0L)
}

#' Commented settings file with every option and its default
#' @noRd
cli_config_template <- function(project_dir = "project", preset = "Ribo-seq", mode = "online") {
  specs <- Filter(function(s) isTRUE(s$file_ok) && s$name != "study_list" && s$name != "config" &&
                    any(c("run", "status") %in% s$commands), cli_options())
  out <- c("# massiveNGSpipe command line settings (mngsp). 'key: value' per line; relative paths are relative to this file.",
           "# Read automatically for study lists in this folder; command-line flags win. 'mngsp config <list>' shows the result.",
           "", paste0("project_dir: ", project_dir), paste0("preset: ", preset), paste0("mode: ", mode))
  for (g in unique(vapply(specs, `[[`, "", "group"))) {
    out <- c(out, "", paste0("# ---- ", g))
    for (s in Filter(function(s) s$group == g, specs)) {
      if (s$name %in% c("project_dir", "preset", "mode")) next
      d <- cli_default_label(s)
      key <- if (s$name == "threads_<step>") "threads_trim" else s$name
      out <- c(out, paste0("# ", s$help), paste0("# ", key, ": ", if (is.null(d)) "" else d))
    }
  }
  out
}

#' mngsp config
#' @noRd
cli_cmd_config <- function(r) {
  v <- r$values
  specs <- Filter(function(s) "config" %in% s$commands && s$name != "threads_<step>", cli_options())
  keys <- unique(c(names(specs), grep("^threads_", names(v), value = TRUE)))
  tab <- data.table(setting = keys,
                    value = vapply(keys, function(k) {
                      if (!is.null(v[[k]])) paste(format(v[[k]]), collapse = ",") else {
                        d <- cli_default_label(cli_spec_lookup(k)); if (is.null(d)) "" else d } }, ""),
                    source = vapply(keys, function(k) if (!is.null(r$source[k]) && !is.na(r$source[k])) r$source[[k]]
                                    else "pipeline_config() default", ""))
  trunc <- function(x, n) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
  cat(sprintf("%-42s %-44s %s", tab$setting, trunc(tab$value, 44), tab$source), sep = "\n")
  if (length(r$files)) message("Settings files: ", paste(r$files, collapse = ", "))
  config <- cli_build_config(v)
  message("Project: ", config$project, " | preset ", config$preset, " (", config$mode, ") | FINAL list: ",
          config$complete_metadata, "\nData folders (ORFik::config()): ",
          paste(names(config$config), config$config, sep = " = ", collapse = ", "),
          "\nWorkers: ", paste(names(config$threads), unlist(config$threads), sep = "=", collapse = " "),
          " | barcode detector: ", barcode_detector_mode(),
          " | Google sheet: ", if (is.null(config$google_url)) "none" else "yes")
  invisible(0L)
}
