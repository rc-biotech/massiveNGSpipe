# Pre-flight checks of `mngsp preflight` and `mngsp run` (detector part ported from the Valen preflight_detector.R).

#' Barcode detector pre-flight
#'
#' Refuses to trim footprint libraries with a detector that would cut wrongly: the detector mode in effect must be
#' "guarded"; the kappa rule, chance guard and adapter check must exist; a synthetic library (29-nt templated inserts
#' with 4 random nt at the 5' end, nothing at the 3' end) must be called 4/0; ORFik must be at least
#' \code{min_orfik}; and a git checkout (when given) must contain \code{required_commits}.
#' @param checkout path of a git checkout loaded with devtools (its history is checked), or NA for the installed
#'  package (its RemoteSha is reported)
#' @param min_orfik character, minimum ORFik version
#' @param required_commits character, commits the checkout must contain
#' @param seed integer, seed of the synthetic library
#' @return one-row data.frame report; stops when a check fails
#' @export
preflight_detector <- function(checkout = NA_character_, min_orfik = "1.29.18",
                               required_commits = c("fbd56c4", "3121b3a"), seed = 7L) {
  fail <- character()
  has_kappa <- exists("end_templating_profile", mode = "function")
  has_chance <- "chance_max" %in% names(formals(end_barcode_call)) && exists("known_adapter_constructs", mode = "function")
  mode <- barcode_detector_mode()
  if (!has_kappa) fail <- c(fail, "end_templating_profile() missing: detector older than be59114 (kappa rule)")
  if (!has_chance) fail <- c(fail, "end_barcode_call() without chance_max / known_adapter_constructs(): detector older than 3121b3a")
  if (mode != "guarded") fail <- c(fail, paste0("barcode detector mode is '", mode, "', must be 'guarded'"))
  call5 <- call3 <- NA_integer_
  if (has_kappa) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()), add = TRUE)
    set.seed(seed)
    b <- c("A", "C", "G", "T")
    tx <- vapply(1:40, function(i) paste(sample(b, 1500, TRUE), collapse = ""), "")
    n <- 30000L; t_i <- sample(40, n, TRUE); pos <- sample(1:1400, n, TRUE); len <- sample(28:31, n, TRUE)
    ins <- substring(tx[t_i], pos, pos + len - 1L)
    flank <- vapply(seq_len(n), function(i) paste(sample(b, 4, TRUE), collapse = ""), "")
    s <- paste0(flank, ins)
    call5 <- as.integer(sum(suppressMessages(end_barcode_call(s, "5p", sample = "preflight"))))
    call3 <- as.integer(sum(suppressMessages(end_barcode_call(s, "3p", full_length = 50L, sample = "preflight"))))
    if (!identical(call5, 4L) || !identical(call3, 0L))
      fail <- c(fail, sprintf("synthetic 4-nt random 5' flank: detector called %s/%s, expected 4/0", call5, call3))
  }
  orfik <- as.character(utils::packageVersion("ORFik"))
  if (utils::compareVersion(orfik, min_orfik) < 0) fail <- c(fail, paste0("ORFik ", orfik, " < ", min_orfik))
  mver <- as.character(utils::packageVersion("massiveNGSpipe"))
  commit <- NA_character_; has_commits <- NA
  if (!is.na(checkout) && file.exists(file.path(checkout, ".git"))) {
    commit <- tryCatch(system2("git", c("-C", checkout, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
                       error = function(e) NA_character_)
    anc <- vapply(required_commits, function(cm) tryCatch(
      system2("git", c("-C", checkout, "merge-base", "--is-ancestor", cm, "HEAD"), stdout = FALSE, stderr = FALSE) == 0,
      error = function(e) NA), logical(1))
    has_commits <- all(anc)
    if (isFALSE(has_commits)) fail <- c(fail, paste0("checkout ", checkout, " (", commit, ") lacks ",
                                                     paste(required_commits[!anc %in% TRUE], collapse = ", ")))
  } else {
    d <- utils::packageDescription("massiveNGSpipe")
    commit <- if (!is.null(d$RemoteSha)) substr(d$RemoteSha, 1, 7) else NA_character_
  }
  rep <- data.frame(time = format(Sys.time(), tz = "UTC", usetz = TRUE), mngsp_commit = commit,
                    contains_required_commits = has_commits, kappa_rule = has_kappa, chance_adapter_guard = has_chance,
                    mngsp_version = mver, detector_mode = mode, synthetic_call = paste0(call5, "/", call3),
                    orfik = orfik, ok = !length(fail))
  if (length(fail)) stop("PRE-FLIGHT FAILED (do not trim footprint libraries):\n  - ", paste(fail, collapse = "\n  - "), call. = FALSE)
  message("pre-flight OK: massiveNGSpipe ", mver, " (", commit, "), guarded detector, synthetic 4/0, ORFik ", orfik)
  rep
}

#' System checks of mngsp preflight
#' @return data.table check, value, ok, note
#' @noRd
cli_system_checks <- function(config, min_open_files = 65535L) {
  nofile <- cli_open_file_limit()
  stack <- tryCatch(system("ulimit -s", intern = TRUE), error = function(e) NA_character_)[1]
  mem <- tryCatch({
    m <- readLines("/proc/meminfo"); as.numeric(gsub("\\D", "", grep("^MemAvailable", m, value = TRUE))) / 1024^2
  }, error = function(e) NA_real_)
  disk <- function(d) tryCatch({
    x <- system2("df", c("-P", "-k", shQuote(d)), stdout = TRUE, stderr = FALSE)
    f <- strsplit(trimws(x[length(x)]), "\\s+")[[1]]; c(free_gb = as.numeric(f[4]) / 1024^2, used = as.numeric(sub("%", "", f[5])))
  }, error = function(e) c(free_gb = NA, used = NA))
  bam <- disk(if (dir.exists(config$config["bam"])) config$config["bam"] else dirname(config$config["bam"]))
  reserve <- config$ram_reserve_gb %||% 0
  rows <- list(
    list("open files (ulimit -n)", format(nofile), !is.na(nofile) && nofile >= min_open_files,
         "STAR needs many open files: 'ulimit -n 524288' (the mngsp wrapper raises it up to 'ulimit -Hn')"),
    list("C stack (ulimit -s)", stack, identical(stack, "unlimited") || isTRUE(as.numeric(stack) >= 65536),
         "64 MB avoids 'C stack usage too close to the limit' in run_pipeline (the wrapper sets it)"),
    list("LC_ALL", if (nzchar(Sys.getenv("LC_ALL"))) Sys.getenv("LC_ALL") else "(unset)", !nzchar(Sys.getenv("LC_ALL")),
         "start R with 'env -u LC_ALL' (the mngsp wrapper unsets it)"),
    list("ORFik data folders", paste(names(config$config), collapse = ","), all(dir.exists(config$config[c("fastq", "bam", "ref")])),
         paste("ORFik::config():", paste(config$config, collapse = " "))),
    list("processed data drive used %", format(bam[["used"]]), isTRUE(bam[["used"]] < config$stop_downloading_new_data_at_drive_usage),
         paste0("downloads pause above ", config$stop_downloading_new_data_at_drive_usage, "%; free GB: ", round(bam[["free_gb"]]))),
    list("MemAvailable GB", format(round(mem)), is.na(mem) || mem > reserve,
         paste0("alignments keep ", reserve, " GB free (ram_reserve_gb)")),
    list("project folder", config$project, dir.exists(config$project), "created by pipeline_config()"),
    list("FINAL list", config$complete_metadata, file.exists(config$complete_metadata), "mngsp curate <study list>"))
  data.table::rbindlist(lapply(rows, function(x) data.table::data.table(check = x[[1]], value = x[[2]], ok = x[[3]], note = x[[4]])))
}

#' mngsp preflight
#' @noRd
cli_cmd_preflight <- function(r) {
  v <- r$values
  config <- cli_build_config(v)
  sys <- cli_system_checks(config, v$min_open_files)
  print(sys)
  det <- tryCatch(preflight_detector(checkout = v$dev %||% (if (nzchar(Sys.getenv("MNGSP_DEV"))) Sys.getenv("MNGSP_DEV") else NA_character_),
                                     min_orfik = v$min_orfik),
                  error = function(e) { message(conditionMessage(e)); NULL })
  if (!is.null(det)) print(det)
  failed <- sys$check[!sys$ok %in% TRUE]
  if (length(failed)) message("Not ready: ", paste(failed, collapse = "; "))
  invisible(if (is.null(det) || length(failed)) 1L else 0L)
}
