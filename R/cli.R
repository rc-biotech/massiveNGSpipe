# Command line interface: `mngsp <command> [options]` (inst/bin/mngsp calls mngsp_cli()).
# Layout: option spec (cli_options) -> parser (cli_parse) -> settings resolution with sources (cli_resolve) ->
# pipeline config (cli_build_config) -> command handlers (cli_cmd_*, study list logic in R/cli_study_list.R).

# ---- commands and options ---------------------------------------------------------------------------------------

cli_list_commands <- c("run", "plan", "status", "curate", "config")

cli_commands <- function() list(
  init      = list(usage = "mngsp init <dir> [--example <accession>] [--preset <p>] [--mode <m>]",
                   summary = "create a project folder with mngsp.conf and an example study list",
                   details = c("Writes <dir>/mngsp.conf (all settings, commented, with their defaults) and <dir>/studies.txt.",
                               "The project (FINAL list, flags, logs) is <dir>/project; data folders come from ORFik::config().")),
  curate    = list(usage = "mngsp curate <study list> [--organism <name>] [--validate-only]",
                   summary = "download metadata of the listed studies and validate it into the FINAL list",
                   details = c("Online mode downloads SRA run tables; rows that pass validation (KEEP = TRUE, unique labels) go to",
                               "the FINAL list. Otherwise edit <project>/next_round_manual.csv and rerun with --validate-only.",
                               "No Google sheet is used unless --google-sheet <url> is given.")),
  plan      = list(usage = "mngsp plan <study list> [options]",
                   summary = "dry run: experiments, references and flags of a study list (writes only the plan CSV)",
                   details = "Same as 'mngsp run' without --apply. --reset is shown in the plan, but nothing is reset."),
  run       = list(usage = "mngsp run <study list> [--apply] [--reset] [options]",
                   summary = "run the pipeline on a study list (dry run unless --apply)",
                   details = c("Without --apply: pre-flight, plan, nothing changed (exit 0).",
                               "With --apply: optional --reset, reference/metadata subset, manual trim tables, provenance,",
                               "then run_pipeline(). Exit status 0 = all done, 2 = finished with error records, 1 = failure.")),
  status    = list(usage = "mngsp status [<study list>] [--session <n>]",
                   summary = "progress: last session checklist, or a stage table per experiment of a study list"),
  errors    = list(usage = "mngsp errors [--session <n>] [--events <n>]",
                   summary = "error records of the last (or n-th last) session, with hints"),
  retry     = list(usage = "mngsp retry <experiment>... | --all [--session <n>] [--note <text>]",
                   summary = "clear error records so a running pipeline retries those experiments"),
  preflight = list(usage = "mngsp preflight [--dev <checkout>]",
                   summary = "system and barcode detector checks (open files, versions, synthetic barcode call, disk, memory)"),
  config    = list(usage = "mngsp config [<study list>] [options]",
                   summary = "show every resolved setting and where it came from (default, file, flag)"),
  resources = list(usage = "mngsp resources show | set <normal|conserve|minimal> --reason <text> [--hours <h>] [--max-mem-gb <n>] [--max-threads <n>] | refresh [--hours <h>] | clear [--set-by <who>]",
                   summary = "show or change the server-wide resource mode (an expiring cap on threads and memory)",
                   details = c("The mode file (default ~/livemount/Bio_data/NGS_pipeline/RESOURCE_MODE, env MNGSP_RESOURCE_MODE_FILE)",
                               "caps every load-based thread count of all pipelines at their next stage start. It expires after",
                               "--hours (at most 6); refresh it to keep it. No file = configured maxima with load-based counts.")),
  help      = list(usage = "mngsp help [<command>]", summary = "show help (also: mngsp <command> --help)"),
  version   = list(usage = "mngsp --version", summary = "print the massiveNGSpipe version")
)

cli_opt <- function(name, type, group, commands, help, target = "cli", arg = name, default = NULL, choices = NULL,
                    file_ok = TRUE) {
  list(name = name, type = type, group = group, commands = commands, help = help, target = target, arg = arg,
       default = default, choices = choices, file_ok = file_ok)
}

#' Option table of the command line interface
#'
#' One row per option: name (flags use dashes, config files underscores), type, the commands it applies to, where it
#' goes (cli = runner logic, pc = argument of \code{pipeline_config()}, field = \code{config$<name>} after
#' \code{pipeline_config()}, threads = \code{config$threads}), and its default. Options with target pc and no default
#' use the \code{pipeline_config()} default.
#' @return named list of option specs
#' @noRd
cli_options <- function() {
  L <- cli_list_commands; RUN <- c("run", "plan", "config"); PIPE <- c(RUN, "status", "curate", "preflight")
  PROJ <- c(PIPE, "errors", "retry")
  o <- list(
    # study list and settings files
    cli_opt("study_list", "file", "Input", L, "study list (alternative to the positional argument): one accession per line, '#' comments", file_ok = FALSE),
    cli_opt("config", "file", "Input", PROJ, "settings file ('key: value' lines); read after the automatic ones", file_ok = FALSE),
    cli_opt("no_config", "flag", "Input", PROJ, "ignore mngsp.conf, <list>.conf and '#:' lines in the study list", default = FALSE, file_ok = FALSE),
    # run control
    cli_opt("apply", "flag", "Run", "run", "execute; without it this is a dry run that only writes the plan", default = FALSE, file_ok = FALSE),
    cli_opt("reset", "flag", "Run", RUN, "full restart of the listed experiments: clears flags, per-sample markers, errors, the study folder and its experiment CSV (raw FASTQ is downloaded again)", default = FALSE, file_ok = FALSE),
    cli_opt("keep", "file", "Run", c(RUN, "status"), "file of studies to leave untouched (never reset, never run)"),
    cli_opt("trim", "file", "Run", RUN, "manual trim table (BioProject, adapter, barcode5p_size, barcode3p_size, mate, reverse, status, evidence): status 'apply' = written to <experiment>/trim/ before trimming, 'hold' = study not run"),
    cli_opt("skip_failed", "chr", "Run", RUN, "skip studies that failed in the last n sessions, or since a date (2026-11-20); 0 = retry all", target = "pc", arg = "skip_failed_sessions"),
    cli_opt("organisms", "file", "Run", c(RUN, "status"), "CSV with BioProject,ScientificName: only these experiments of multi-organism studies"),
    cli_opt("manual_dir", "dir", "Run", RUN, "manual trim tables per experiment (<dir>/<experiment>/{barcodes,adapters,mate_override}_manual.csv), copied into <experiment>/trim/ after a reset"),
    cli_opt("exclude_runs", "file", "Run", c(RUN, "status"), "CSV with a Run column: runs never processed (KEEP = FALSE runs, restricted access)"),
    cli_opt("exclude_studies", "list", "Run", c(RUN, "status"), "comma-separated studies never processed"),
    cli_opt("skip_dbgap", "lgl", "Run", c(RUN, "status"), "skip known dbGaP controlled-access studies (dbgap_studies())", default = TRUE),
    cli_opt("download_references", "flag", "Run", RUN, "download genomes and build STAR indices that are missing (pipeline_init_all); default: skip studies without a complete reference", default = FALSE),
    cli_opt("preflight", "lgl", "Run", RUN, "barcode detector pre-flight before footprint presets (Ribo-seq, disome, SSU); --no-preflight skips it", default = TRUE),
    cli_opt("dev", "dir", "Run", c(RUN, "preflight"), "checkout whose git history the pre-flight checks (default: $MNGSP_DEV)"),
    cli_opt("min_orfik", "chr", "Run", c(RUN, "preflight"), "minimum ORFik version of the pre-flight", default = "1.29.18"),
    cli_opt("min_open_files", "int", "Run", c(RUN, "preflight"), "minimum 'ulimit -n' for --apply (STAR opens many files)", default = 65535L),
    cli_opt("wait", "num", "Run", RUN, "seconds between run_pipeline() polling rounds", default = 10),
    cli_opt("out_dir", "path", "Run", RUN, "folder of the plan, subset metadata and pre-flight CSVs (default: the project folder)"),
    cli_opt("file_prefix", "chr", "Run", RUN, "prefix of the plan and pre-flight CSV names", default = "mngsp_"),
    cli_opt("subset_suffix", "chr", "Run", RUN, "suffix of subset_<list><suffix>.csv (the Valen sync scripts read '_v4')", default = ""),
    cli_opt("reference_hint", "chr", "Run", RUN, "message for missing references; {folders} = their folder names", default = "add --download-references, or copy them: {folders}"),
    cli_opt("metadata_hint", "chr", "Run", RUN, "message for studies missing from the FINAL list", default = "curate them first: mngsp curate <study list>"),
    # project
    cli_opt("project_dir", "path", "Project", PROJ, "project folder (FINAL list, flags, logs, errors)", target = "pc"),
    cli_opt("preset", "choice", "Project", PROJ, "pipeline preset", target = "pc", choices = c("Ribo-seq", "RNA-seq", "disome", "SSU")),
    cli_opt("mode", "choice", "Project", PROJ, "online = public SRA/ENA studies, local = your own FASTQ folders", target = "pc", choices = c("online", "local")),
    cli_opt("complete_metadata", "path", "Project", PROJ, "FINAL list CSV the pipeline runs", target = "pc"),
    cli_opt("google_sheet", "chr", "Project", PROJ, "Google sheet URL of the curated metadata, or 'none'. Default: the project's google_api_url_table.csv if it exists, else none (the CLI never creates sheets)"),
    cli_opt("discord", "lgl", "Project", PROJ, "Discord progress messages (when a discordr connection is cached)", default = TRUE),
    cli_opt("organism", "list", "Project", "curate", "curate only these organisms (comma-separated scientific names); default all"),
    cli_opt("validate_only", "flag", "Project", "curate", "only validate <project>/next_round_manual.csv (after editing it)", default = FALSE),
    # pipeline settings
    cli_opt("delete_raw_files", "lgl", "Pipeline", PIPE, "delete raw FASTQ after trimming", target = "pc"),
    cli_opt("delete_trimmed_files", "lgl", "Pipeline", PIPE, "delete trimmed FASTQ", target = "pc"),
    cli_opt("delete_collapsed_files", "lgl", "Pipeline", PIPE, "delete collapsed FASTA", target = "pc"),
    cli_opt("split_unique_mappers", "lgl", "Pipeline", PIPE, "also produce unique-mapper outputs", target = "pc"),
    cli_opt("all_mappers", "lgl", "Pipeline", PIPE, "produce all-mapper outputs", target = "pc"),
    cli_opt("fast_trim", "lgl", "Pipeline", PIPE, "single-pass C trim + collapse for single-end reads", target = "pc"),
    cli_opt("trim_parallel", "choice", "Pipeline", PIPE, "trim stage: several libraries at once (staged = fastest) or one at a time", target = "pc", choices = c("staged", "naive", "serial")),
    cli_opt("trim_workers", "int", "Resources", PIPE, "libraries trimmed at once (default: from idle CPUs and free memory)", target = "pc"),
    cli_opt("demultiplex", "lgl", "Pipeline", PIPE, "demultiplex pooled lanes first (see ?demultiplex)", target = "pc"),
    cli_opt("contam", "lgl", "Pipeline", PIPE, "contaminant depletion step", target = "pc"),
    cli_opt("adapter_detection", "choice", "Pipeline", PIPE, "adapter detection method", target = "pc", choices = c("scan", "fastqc")),
    cli_opt("barcode_detector", "choice", "Pipeline", PIPE, "barcode/UMI detector mode", default = "guarded", choices = c("guarded", "legacy")),
    cli_opt("auto_footprint_mate", "lgl", "Pipeline", PIPE, "paired footprint libraries: process one mate as single-end", target = "pc"),
    cli_opt("pshift_per_library", "lgl", "Pipeline", PIPE, "per-library P-shift fallback when the study-wide pass fails", target = "pc"),
    cli_opt("reuse_shifts_if_existing", "lgl", "Pipeline", PIPE, "reuse a valid shifting_table.rds", target = "pc"),
    cli_opt("accepted_lengths_rpf", "ints", "Pipeline", PIPE, "read lengths to P-shift, e.g. 20,21,25:33", target = "pc"),
    cli_opt("max_unprocessed_downloads", "int", "Pipeline", PIPE, "pause downloads while this many studies wait for trimming", target = "pc"),
    cli_opt("stop_downloading_new_data_at_drive_usage", "num", "Pipeline", PIPE, "pause downloads above this drive usage (%)", target = "pc"),
    cli_opt("compress_raw_data", "lgl", "Pipeline", PIPE, "compress raw FASTQ", target = "pc"),
    cli_opt("keep_contaminants", "lgl", "Pipeline", PIPE, "keep contaminant-aligned reads", target = "pc"),
    cli_opt("count_cache", "lgl", "Pipeline", PIPE, "reuse validated per-sample coverage counts", target = "pc"),
    cli_opt("ssd_scratch_dir", "chr", "Pipeline", PIPE, "local SSD scratch for alignment, or 'none'", target = "pc"),
    cli_opt("ssd_min_free_gb", "num", "Pipeline", PIPE, "free-space margin on the SSD scratch", target = "pc"),
    # threads and memory
    cli_opt("threads", "int", "Resources", PIPE, "workers per stage (threads_default)", target = "pc", arg = "threads_default"),
    cli_opt("threads_<step>", "int", "Resources", PIPE, "workers of one step, e.g. --threads-trim 8, --threads-collapse 16, --threads-pshifted 32, --threads-pcounts 1, --threads-ofst 8", target = "threads"),
    cli_opt("blas_cap", "int", "Resources", PIPE, "BLAS/OpenMP threads per worker and cap of pshifted/valid_pshift/pcounts workers", target = "pc", arg = "threads_blas_cap"),
    cli_opt("dynamic_thread_count", "lgl", "Resources", PIPE, "thread counts load-based within the configured maxima and capped by the resource mode file; --no-dynamic-thread-count uses the configured values", target = "pc"),
    cli_opt("thread_type", "choice", "Resources", PIPE, "BiocParallel backend (default multicore)", choices = c("multicore", "serial")),
    cli_opt("max_ram_wait_minutes", "num", "Resources", PIPE, "how long an alignment may wait for free memory", target = "pc"),
    cli_opt("ram_reserve_gb", "num", "Resources", PIPE, "memory every alignment must leave available (default: shared budget reserve, else 0)", target = "field"),
    cli_opt("bam_sort_memory_gb", "num", "Resources", PIPE, "cap of the per-job STAR BAM-sort memory", target = "field"),
    cli_opt("external_bam_sort_gb", "num", "Resources", PIPE, "budget of the external BAM sorter (default 16)", target = "field"),
    cli_opt("ram_wait_gap_minutes", "num", "Resources", PIPE, "a sample not rechecked for this long starts a new memory wait (default 30)", target = "field"),
    cli_opt("ram_recheck_seconds", "num", "Resources", PIPE, "recheck interval of memory-deferred samples (default 300)", target = "field"),
    # resources
    cli_opt("hours", "num", "Resources mode", "resources", "hours until the mode expires (max 6; refresh to keep it)", default = 2, file_ok = FALSE),
    cli_opt("reason", "chr", "Resources mode", "resources", "why the mode is set", file_ok = FALSE),
    cli_opt("max_mem_gb", "num", "Resources mode", "resources", "cap of the memory one step may plan to use", file_ok = FALSE),
    cli_opt("max_threads", "int", "Resources mode", "resources", "cap of the threads / workers of any step", file_ok = FALSE),
    cli_opt("set_by", "chr", "Resources mode", "resources", "writer name (set); clear only a mode set by this writer (clear)", file_ok = FALSE),
    cli_opt("mode_file", "chr", "Resources mode", "resources", "resource mode file (default: MNGSP_RESOURCE_MODE_FILE or the shared path)", file_ok = FALSE),
    # status / errors / retry / init
    cli_opt("session", "int", "Status", c("status", "errors", "retry"), "session index, newest first", default = 1L, file_ok = FALSE),
    cli_opt("events", "int", "Status", "errors", "also show the last n events of error_logs/events.jsonl", default = 0L, file_ok = FALSE),
    cli_opt("all", "flag", "Status", "retry", "clear every error record of the session", default = FALSE, file_ok = FALSE),
    cli_opt("note", "chr", "Status", "retry", "reason stored with the cleared records", default = "", file_ok = FALSE),
    cli_opt("example", "chr", "Init", "init", "example study written to studies.txt", default = "PRJNA779948", file_ok = FALSE),
    cli_opt("force", "flag", "Init", "init", "overwrite an existing mngsp.conf / studies.txt", default = FALSE, file_ok = FALSE)
  )
  names(o) <- vapply(o, `[[`, "", "name")
  o
}

cli_flag_name <- function(name) paste0("--", gsub("_", "-", name))

cli_spec_lookup <- function(name, specs = cli_options()) {
  if (!is.null(specs[[name]]) && name != "threads_<step>") { s <- specs[[name]]; s$key <- name; return(s) }
  if (grepl("^threads_[a-z][a-z_]*$", name)) {
    s <- specs[["threads_<step>"]]; s$key <- name; s$step <- sub("^threads_", "", name); return(s)
  }
  NULL
}

cli_nearest <- function(x, choices, max = 2) {
  d <- utils::adist(x, choices)[1, ]
  choices[d <= max & d == min(d)]
}

cli_as_logical <- function(x, strict = FALSE) {
  x <- tolower(trimws(as.character(x)))
  ok <- if (strict) c("true", "false") else c("true", "false", "t", "f", "yes", "no", "1", "0", "on", "off")
  if (length(x) != 1 || !(x %in% ok)) return(NA)
  x %in% c("true", "t", "yes", "1", "on")
}

cli_convert <- function(spec, value, where) {
  value <- trimws(value)
  bad <- function(what) stop(sprintf("%s: %s expects %s, got '%s'", where, cli_flag_name(spec$key %||% spec$name),
                                     what, value), call. = FALSE)
  switch(spec$type,
    flag = , lgl = { v <- cli_as_logical(value); if (is.na(v)) bad("TRUE or FALSE"); v },
    int = { v <- suppressWarnings(as.numeric(value)); if (is.na(v) || v != round(v) || v < 0) bad("a whole number"); as.integer(v) },
    num = { v <- suppressWarnings(as.numeric(value)); if (is.na(v)) bad("a number"); v },
    ints = {
      parts <- trimws(strsplit(value, ",")[[1]])
      v <- unlist(lapply(parts, function(p) {
        r <- suppressWarnings(as.integer(strsplit(p, ":", fixed = TRUE)[[1]]))
        if (!length(r) || anyNA(r) || length(r) > 2) bad("integers like 20,21,25:33")
        if (length(r) == 2) seq(r[1], r[2]) else r
      }))
      v
    },
    list = { v <- trimws(strsplit(value, ",")[[1]]); v[nzchar(v)] },
    choice = { if (!(value %in% spec$choices)) bad(paste0("one of ", paste(spec$choices, collapse = ", "))); value },
    value)
}

# ---- parsing --------------------------------------------------------------------------------------------------------

#' Parse mngsp command line arguments
#'
#' \code{--name value}, \code{--name=value}, dashes or underscores (\code{--manual-dir} = \code{--manual_dir}),
#' \code{--flag}, \code{--flag TRUE}, \code{--no-flag}. Positional arguments are the study list (run, plan, status,
#' curate, config), the experiments (retry) or the folder (init).
#' @param args character vector, as from \code{commandArgs(TRUE)}
#' @return list(command, positional, opts, help)
#' @noRd
cli_parse <- function(args) {
  args <- as.character(args)
  cmds <- cli_commands()
  if (!length(args) || args[1] %in% c("-h", "--help")) return(list(command = "help", positional = args[-1][1][!is.na(args[-1][1])], opts = list(), help = TRUE))
  if (args[1] %in% c("-v", "-V", "--version")) return(list(command = "version", positional = character(), opts = list(), help = FALSE))
  cmd <- args[1]
  if (!(cmd %in% names(cmds))) {
    if (startsWith(cmd, "-")) stop("the command comes first: mngsp <command> [options] (see mngsp --help)", call. = FALSE)
    near <- cli_nearest(cmd, names(cmds))
    stop("unknown command '", cmd, "'", if (length(near)) paste0(" (did you mean: ", paste(near, collapse = ", "), "?)"),
         "; see mngsp --help", call. = FALSE)
  }
  specs <- cli_options()
  rest <- args[-1]; pos <- character(); opts <- list(); help <- FALSE; i <- 1L
  while (i <= length(rest)) {
    a <- rest[i]
    if (a %in% c("-h", "--help")) { help <- TRUE; i <- i + 1L; next }
    if (a == "--") { pos <- c(pos, rest[-seq_len(i)]); break }
    if (startsWith(a, "--")) {
      body <- substring(a, 3); val <- NULL
      if (grepl("=", body, fixed = TRUE)) { val <- sub("^[^=]*=", "", body); body <- sub("=.*$", "", body) }
      name <- gsub("-", "_", body); negate <- FALSE
      spec <- cli_spec_lookup(name, specs)
      if (is.null(spec) && startsWith(name, "no_")) {
        spec <- cli_spec_lookup(substring(name, 4), specs)
        if (!is.null(spec) && spec$type %in% c("flag", "lgl")) negate <- TRUE else spec <- NULL
      }
      if (is.null(spec)) {
        near <- cli_nearest(name, setdiff(names(specs), "threads_<step>"))
        stop("unknown option ", a, if (length(near)) paste0(" (did you mean ", paste(cli_flag_name(near), collapse = ", "), "?)"),
             "; see mngsp ", cmd, " --help", call. = FALSE)
      }
      if (!(cmd %in% spec$commands))
        stop(cli_flag_name(spec$key), " does not apply to 'mngsp ", cmd, "' (it is used by: ",
             paste(spec$commands, collapse = ", "), ")", call. = FALSE)
      if (spec$type %in% c("flag", "lgl")) {
        if (negate) {
          if (!is.null(val)) stop(a, " takes no value", call. = FALSE)
          v <- FALSE
        } else {
          if (is.null(val) && i < length(rest) && !is.na(cli_as_logical(rest[i + 1L], strict = TRUE))) {
            val <- rest[i + 1L]; i <- i + 1L
          }
          v <- if (is.null(val)) TRUE else cli_convert(spec, val, "command line")
        }
      } else {
        if (is.null(val)) {
          if (i == length(rest) || (startsWith(rest[i + 1L], "--") && nchar(rest[i + 1L]) > 2))
            stop(cli_flag_name(spec$key), " needs a value", call. = FALSE)
          val <- rest[i + 1L]; i <- i + 1L
        }
        v <- cli_convert(spec, val, "command line")
      }
      opts[[spec$key]] <- v
      i <- i + 1L; next
    }
    if (startsWith(a, "-") && nchar(a) > 1) stop("unknown option ", a, "; options are long (--name)", call. = FALSE)
    pos <- c(pos, a); i <- i + 1L
  }
  list(command = cmd, positional = pos, opts = opts, help = help)
}

# ---- settings files and resolution ----------------------------------------------------------------------------------

#' Read a mngsp settings file
#'
#' Format: one \code{key: value} (or \code{key = value}) per line, \code{#} comments, keys as the long options with
#' dashes or underscores. With \code{inline = TRUE}, only \code{#: key: value} lines are read (settings inside a
#' study list). Relative paths are relative to the file's folder.
#' @param path file
#' @param inline logical
#' @param command the command the settings are for (keys of other commands are ignored)
#' @return list(values, source)
#' @noRd
cli_read_config <- function(path, inline = FALSE, command = "run") {
  lines <- readLines(path, warn = FALSE)
  lnum <- seq_along(lines)
  if (inline) {
    keep <- grepl("^#:", lines); lines <- sub("^#:\\s*", "", lines[keep]); lnum <- lnum[keep]
  }
  specs <- cli_options(); values <- list(); where <- character()
  for (k in seq_along(lines)) {
    l <- sub("\\s+#.*$", "", trimws(lines[k]))
    if (!nzchar(l) || startsWith(l, "#")) next
    if (!grepl("^[A-Za-z][A-Za-z0-9_<>.-]*\\s*[:=]", l))
      stop(path, ":", lnum[k], ": expected 'key: value', got '", l, "'", call. = FALSE)
    key <- gsub("-", "_", trimws(sub("\\s*[:=].*$", "", l)))
    val <- trimws(sub("^[^:=]*[:=]\\s*", "", l))
    val <- sub("^([\"'])(.*)\\1$", "\\2", val)
    spec <- cli_spec_lookup(key, specs)
    loc <- paste0(path, ":", lnum[k])
    if (is.null(spec)) stop(loc, ": unknown setting '", key, "' (see mngsp help)", call. = FALSE)
    if (!isTRUE(spec$file_ok))
      stop(loc, ": '", key, "' can only be given on the command line (", cli_flag_name(key), ")", call. = FALSE)
    if (!(command %in% spec$commands)) next
    v <- cli_convert(spec, val, loc)
    if (spec$type %in% c("file", "dir", "path")) v <- cli_abs_path(v, dirname(path))
    if (spec$key == "study_list") next
    values[[spec$key]] <- v; where[[spec$key]] <- loc
  }
  list(values = values, source = where)
}

cli_abs_path <- function(p, base = getwd()) {
  p <- path.expand(p)
  if (!grepl("^/", p)) p <- file.path(base, p)
  normalizePath(p, mustWork = FALSE)
}

#' Resolve settings: defaults < automatic files < --config < flags
#' @return list(command, positional, values, source, study_list, files)
#' @noRd
cli_resolve <- function(parsed, cwd = getwd()) {
  cmd <- parsed$command; opts <- parsed$opts; specs <- cli_options()
  pos <- parsed$positional
  list_file <- NULL
  if (cmd %in% cli_list_commands) {
    if (!is.null(opts$study_list) && length(pos)) stop("give the study list once (positional or --study-list)", call. = FALSE)
    if (length(pos) > 1) stop("'mngsp ", cmd, "' takes one study list, got: ", paste(pos, collapse = " "), call. = FALSE)
    list_file <- opts$study_list %||% (if (length(pos)) pos[1] else NULL)
    if (!is.null(list_file)) {
      list_file <- cli_abs_path(list_file, cwd)
      if (!file.exists(list_file)) stop("study list not found: ", list_file, call. = FALSE)
    }
    if (cmd %in% c("run", "plan", "curate") && is.null(list_file))
      stop("'mngsp ", cmd, "' needs a study list; see mngsp ", cmd, " --help", call. = FALSE)
  }
  # defaults of the options of this command
  values <- list(); source <- character()
  for (s in specs) if (cmd %in% s$commands && !is.null(s$default)) { values[[s$name]] <- s$default; source[[s$name]] <- "default" }
  layers <- list(); files <- character()
  if (!isTRUE(opts$no_config)) {
    cand <- if (!is.null(list_file)) c(file.path(dirname(list_file), "mngsp.conf"),
                                      paste0(tools::file_path_sans_ext(list_file), ".conf")) else file.path(cwd, "mngsp.conf")
    for (f in unique(cand)) if (file.exists(f)) { layers <- c(layers, list(cli_read_config(f, command = cmd))); files <- c(files, f) }
    if (!is.null(list_file) && any(grepl("^#:", readLines(list_file, warn = FALSE)))) {
      layers <- c(layers, list(cli_read_config(list_file, inline = TRUE, command = cmd))); files <- c(files, list_file)
    }
  }
  if (!is.null(opts$config)) {
    f <- cli_abs_path(opts$config, cwd)
    if (!file.exists(f)) stop("settings file not found: ", f, call. = FALSE)
    layers <- c(layers, list(cli_read_config(f, command = cmd))); files <- c(files, f)
  }
  for (l in layers) for (k in names(l$values)) { values[[k]] <- l$values[[k]]; source[[k]] <- l$source[[k]] }
  for (k in names(opts)) {
    v <- opts[[k]]
    s <- cli_spec_lookup(k, specs)
    if (s$type %in% c("file", "dir", "path")) v <- cli_abs_path(v, cwd)
    values[[k]] <- v; source[[k]] <- "command line"
  }
  for (k in names(values)) {
    s <- cli_spec_lookup(k, specs)
    if (s$type == "file" && !file.exists(values[[k]])) stop(cli_flag_name(k), ": file not found: ", values[[k]], " (", source[[k]], ")", call. = FALSE)
    if (s$type == "dir" && !dir.exists(values[[k]])) stop(cli_flag_name(k), ": folder not found: ", values[[k]], " (", source[[k]], ")", call. = FALSE)
  }
  if (cmd == "plan" && isTRUE(values$apply)) stop("'mngsp plan' is always a dry run; use 'mngsp run --apply'", call. = FALSE)
  list(command = cmd, positional = pos, values = values, source = source, study_list = list_file, files = files)
}

# ---- pipeline config ------------------------------------------------------------------------------------------------

#' Build the pipeline config from resolved CLI settings
#'
#' Settings without a value keep the \code{pipeline_config()} default. Differences: no Google sheet unless the
#' project already has one (or \code{google_sheet} is set), so the CLI never creates sheets.
#' @param values named list, resolved settings
#' @param data_config named character, ORFik data folders (fastq, bam, ref, exp)
#' @return the config list of \code{pipeline_config()}
#' @noRd
cli_build_config <- function(values, data_config = ORFik::config(), verbose = FALSE) {
  specs <- cli_options(); args <- list(config = data_config, verbose = verbose)
  for (k in names(values)) {
    s <- cli_spec_lookup(k, specs)
    if (identical(s$target, "pc")) args[[s$arg]] <- values[[k]]
  }
  if (!is.null(args$skip_failed_sessions)) {
    x <- suppressWarnings(as.numeric(args$skip_failed_sessions))
    if (!is.na(x)) args$skip_failed_sessions <- x
    else if (is.na(suppressWarnings(as.POSIXct(args$skip_failed_sessions, optional = TRUE))))
      stop("--skip-failed expects a number of sessions or a date, got '", args$skip_failed_sessions, "'", call. = FALSE)
  }
  if (identical(args$ssd_scratch_dir, "none")) args["ssd_scratch_dir"] <- list(NULL)
  project <- args$project_dir %||% file.path(dirname(data_config)[1], "NGS_pipeline")
  gs <- values$google_sheet
  if (is.null(gs)) {
    if (!file.exists(file.path(project, "google_api_url_table.csv"))) args["google_url"] <- list(NULL)
  } else if (tolower(gs) %in% c("none", "false", "")) args["google_url"] <- list(NULL) else args$google_url <- gs
  if (isFALSE(values$discord)) args["discord_webhook"] <- list(NULL)
  if (!is.null(values$thread_type))
    args$thread_type <- if (values$thread_type == "serial") BiocParallel::SerialParam else BiocParallel::MulticoreParam
  if (!is.null(values$barcode_detector)) options(massiveNGSpipe.barcode_detector = values$barcode_detector)
  config <- do.call(pipeline_config, args)
  if (!is.null(values$barcode_detector)) config$barcode_detector <- values$barcode_detector
  for (k in names(values)) {
    s <- cli_spec_lookup(k, specs)
    if (identical(s$target, "field")) config[[k]] <- values[[k]]
    if (identical(s$target, "threads")) {
      if (!(s$step %in% c(names(config$threads), names(config$flag), "ofst")))
        warning("--threads-", s$step, ": no step '", s$step, "' in this preset (steps: ",
                paste(names(config$threads), collapse = ", "), ")", call. = FALSE)
      config$threads[[s$step]] <- values[[k]]
    }
  }
  if (!identical(basename(config$project), "NGS_pipeline") || config$preset == "Ribo-seq") return(config)
  warning("preset ", config$preset, " uses the Ribo-seq project folder ", config$project,
          "; set --project-dir (e.g. NGS_pipeline_RNA) to keep modalities apart", call. = FALSE)
  config
}

# ---- help -----------------------------------------------------------------------------------------------------------

cli_default_label <- function(s) {
  if (!is.null(s$default)) return(paste(format(s$default), collapse = ","))
  if (identical(s$target, "pc")) {
    f <- formals(pipeline_config)[[s$arg]]
    if (is.null(f)) return(NULL)
    if (is.call(f) && identical(f[[1]], as.name("[")))   # c("online", "local")[1] -> "online"
      return(format(tryCatch(eval(f, baseenv()), error = function(e) deparse(f))))
    if (!is.call(f)) return(paste(format(f), collapse = ","))
    return(paste(gsub("\\s+", " ", deparse(f)), collapse = " "))
  }
  NULL
}

cli_wrap <- function(text, indent, width = 100) {
  out <- strwrap(text, width = width - indent)
  paste0(c("", rep(strrep(" ", indent), max(0, length(out) - 1))), out, collapse = "\n")
}

cli_help <- function(command = NULL) {
  cmds <- cli_commands()
  ver <- tryCatch(as.character(utils::packageVersion("massiveNGSpipe")), error = function(e) "dev")
  if (is.null(command) || !(command %in% names(cmds))) {
    w <- max(nchar(names(cmds)))
    lines <- c(paste0("mngsp: command line interface of massiveNGSpipe ", ver),
               "", "Usage: mngsp <command> [options]   (mngsp <command> --help for its options)", "", "Commands:",
               sprintf("  %-*s  %s", w, names(cmds), vapply(cmds, `[[`, "", "summary")), "",
               "Settings, later wins: pipeline_config() defaults < mngsp.conf next to the study list < <list>.conf",
               "  < '#: key: value' lines in the study list < --config <file> < command-line flags ('mngsp config <list>' shows them).",
               "", "Typical use:",
               "  mngsp init ~/mngsp_demo && cd ~/mngsp_demo",
               "  mngsp curate studies.txt --organism 'Homo sapiens'",
               "  mngsp plan studies.txt",
               "  mngsp run studies.txt --apply --download-references",
               "  mngsp status studies.txt; mngsp errors; mngsp retry <experiment>",
               "  mngsp run studies.txt --apply --reset          # re-trim from scratch")
    return(paste(lines, collapse = "\n"))
  }
  cm <- cmds[[command]]
  specs <- Filter(function(s) command %in% s$commands, cli_options())
  lines <- c(paste0("Usage: ", cm$usage), "", cli_wrap(cm$summary, 0), if (!is.null(cm$details)) c("", cm$details))
  if (length(specs)) {
    groups <- unique(vapply(specs, `[[`, "", "group"))
    for (g in groups) {
      lines <- c(lines, "", paste0(g, ":"))
      for (s in Filter(function(s) s$group == g, specs)) {
        meta <- switch(s$type, flag = "", lgl = " [TRUE|FALSE]", choice = paste0(" <", paste(s$choices, collapse = "|"), ">"),
                       ints = " <n,n:m>", list = " <a,b>", int = " <n>", num = " <x>", file = " <file>", dir = " <dir>",
                       path = " <path>", " <text>")
        d <- cli_default_label(s)
        txt <- paste0(s$help, if (!is.null(d) && s$type != "flag") paste0(" [default: ", d, "]"))
        flag <- paste0("  ", cli_flag_name(s$name), meta)
        lines <- c(lines, if (nchar(flag) > 34) c(flag, paste0(strrep(" ", 36), cli_wrap(txt, 36)))
                          else paste0(sprintf("%-36s", flag), cli_wrap(txt, 36)))
      }
    }
  }
  paste(lines, collapse = "\n")
}

# ---- entry point ----------------------------------------------------------------------------------------------------

#' massiveNGSpipe command line interface
#'
#' Entry point of the \code{mngsp} shell command (\code{inst/bin/mngsp}, put on the PATH with
#' \code{\link{install_cli}}). Subcommands: init, curate, plan, run, status, errors, retry, preflight, config.
#' Run \code{mngsp --help} or \code{mngsp <command> --help} for the options. Settings come from
#' \code{pipeline_config()} defaults, a \code{mngsp.conf} next to the study list, \code{<list>.conf},
#' \code{#: key: value} lines in the study list, \code{--config <file>} and flags (later wins).
#' @param args character vector of arguments, default \code{commandArgs(TRUE)}
#' @param cwd folder that relative paths are relative to
#' @return integer exit status, invisibly: 0 ok, 1 error, 2 finished with pipeline error records
#' @export
#' @examples
#' mngsp_cli("--help")
#' mngsp_cli(c("run", "--help"))
#' \dontrun{
#' mngsp_cli(c("plan", "studies.txt"))
#' mngsp_cli(c("run", "studies.txt", "--apply", "--threads", "24"))
#' }
mngsp_cli <- function(args = commandArgs(TRUE), cwd = getwd()) {
  status <- tryCatch({
    p <- cli_parse(args)
    if (p$command == "version") { cat("massiveNGSpipe", as.character(utils::packageVersion("massiveNGSpipe")), "\n"); return(invisible(0L)) }
    if (p$command == "help" || isTRUE(p$help)) {
      topic <- if (p$command == "help") p$positional[1] else p$command
      cat(cli_help(if (length(topic) && !is.na(topic)) topic else NULL), "\n", sep = "")
      return(invisible(0L))
    }
    r <- cli_resolve(p, cwd)
    switch(r$command,
           init = cli_cmd_init(r, cwd), curate = cli_cmd_curate(r), plan = cli_cmd_run(r, apply = FALSE),
           run = cli_cmd_run(r, apply = isTRUE(r$values$apply)), status = cli_cmd_status(r),
           errors = cli_cmd_errors(r), retry = cli_cmd_retry(r), preflight = cli_cmd_preflight(r),
           config = cli_cmd_config(r), resources = cli_cmd_resources(r))
  }, error = function(e) {
    message("mngsp: error: ", conditionMessage(e))
    1L
  })
  invisible(as.integer(status))
}

#' mngsp resources show|set|refresh|clear
#' @noRd
cli_cmd_resources <- function(r) {
  v <- r$values; pos <- r$positional
  action <- if (length(pos)) pos[1] else "show"
  path <- if (!is.null(v$mode_file)) path.expand(v$mode_file) else resource_mode_file()
  switch(action,
         show = resource_mode_show(path),
         set = {
           if (length(pos) < 2) stop("mngsp resources set <normal|conserve|minimal> --reason <text>", call. = FALSE)
           args <- list(mode = pos[2], hours = v$hours %||% 2, reason = v$reason %||% "",
                        max_mem_gb = v$max_mem_gb %||% NA, max_threads = v$max_threads %||% NA, path = path)
           if (!is.null(v$set_by)) args$set_by <- v$set_by
           do.call(resource_mode_set, args)
           resource_mode_show(path)
         },
         refresh = { resource_mode_refresh(hours = v$hours %||% 2, path = path); resource_mode_show(path) },
         clear = resource_mode_clear(only_set_by = v$set_by, path = path),
         stop("mngsp resources: unknown action '", action, "' (show, set, refresh, clear)", call. = FALSE))
  invisible(0L)
}
