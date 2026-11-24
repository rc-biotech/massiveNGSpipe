# Thread, worker and memory decisions of the whole package, plus the optional server-wide resource mode file.
#
# One rule (Håkon, 2026-11-24):
#   * Every thread / worker count goes through resource_threads(step, max = <configured value>). The configured value
#     is the MAXIMUM. With dynamic thread counts on (pipeline_config(dynamic_thread_count = TRUE), the default) the
#     count is load-based: idle CPUs (cores - 1-min load average) and MemAvailable above the reserve, bounded to
#     min..max (memory may go below min, down to 1).
#   * The resource mode file (RESOURCE_MODE, see resource_mode_file()) is a master regulator on top: when present and
#     valid it caps the load-based result further (never raises it). Absent, expired or invalid = no extra cap.
#   * dynamic_thread_count = FALSE: the configured maximum, always (load and mode file ignored).
# Step names are registered in resource_steps(); unknown steps get the "default" row, so new code can just call
# resource_threads("my_step", max = n).

.resource_env <- new.env(parent = emptyenv())

# ---- step registry -------------------------------------------------------------------------------------------------

#' Registered thread / worker steps
#'
#' One row per step: how many cores one worker (or thread) uses (\code{cores_per_worker}, 0 = no CPU check),
#' its memory in GB (\code{per_worker_gb}, NA = no memory check), the busy-host floor \code{min} (the CPU check never
#' goes below it; a memory shortage can) and what happens when /proc is unreadable (\code{fallback}: "max" or "1").
#' Unknown steps use the "default" row. \code{min} is capped at the requested maximum.
#' @return data.frame
#' @export
resource_steps <- function() {
  s <- function(step, cores_per_worker, per_worker_gb, min, fallback, where, what)
    data.frame(step = step, cores_per_worker = cores_per_worker, per_worker_gb = per_worker_gb, min = min,
               fallback = fallback, where = where, what = what, stringsAsFactors = FALSE)
  rbind(
    s("default",        1, NA,  4, "max", "bpparam_from_config()", "BiocParallel workers of any config$threads step without its own row"),
    s("main",           0, NA,  1, "max", "run_pipeline(): config$threads$main", "stage groups run at once (mostly waiting; no CPU check)"),
    s("trim",           1, NA,  4, "max", "pipeline_trim(): config$threads$trim", "BiocParallel workers of the trim stage helpers"),
    s("collapse",       1, NA,  4, "max", "pipeline_collapse(): config$threads$collapse", "collapse workers (after memory_safe_worker_count())"),
    s("collapse_fast",  1, NA,  4, "max", "pipeline_collapse(): fast_collapse_fastq()", "threads of the trimcollapse binary when collapsing"),
    s("pshifted",       1, NA,  4, "max", "config$threads$pshifted", "P-shift workers"),
    s("valid_pshift",   1, NA,  4, "max", "config$threads$valid_pshift", "P-shift validation workers"),
    s("pcounts",        1, NA,  4, "max", "config$threads$pcounts", "count table workers"),
    s("ofst",           1, NA,  1, "max", "ofst_conversion_param()", "BAM -> ofst sample workers (after ofst_worker_count())"),
    s("blas",           1, NA,  4, "max", "pipeline_config(), bpparam_from_config()", "BLAS/OpenMP threads per process (threads_blas_cap)"),
    s("trim_workers",   4,  8,  2, "1",   "trim_workers_dynamic()", "libraries trimmed at once (parallel trim)"),
    s("ofst_derive",    1,  2,  4, "dt",  "derive_threads_dynamic()", "ORFik ofst_derive threads (ORFik.derive_threads)"),
    s("star",           1, NA,  8, "max", "pipeline_align_one_organism(): STAR --runThreadN", "STAR mapping threads (ORFik max.cpus)"),
    s("star_sort",      1, NA,  6, "max", "star_sort_threads(): --outBAMsortingThreadN", "STAR BAM sort/compress threads"),
    s("trimcollapse",   1, NA,  4, "max", "trimcollapse_run()", "threads of the trimcollapse binary when trimming"),
    s("fastp",          1, NA,  2, "max", "anchored/barcode front trim, trim sample", "fastp / trimcollapse threads of helper trims"),
    s("fasterq_dump",   1, NA,  4, "max", "download_sra(): fasterq-dump --threads", "fasterq-dump threads"),
    s("io",             1, NA,  1, "max", "qs2 nthreads / fst::threads_fst", "compression threads of qs2 / fst writes"),
    s("demultiplex",    1, NA,  1, "max", "pipeline_demultiplex()", "files / collapse workers of the fast demultiplexer"),
    s("qc_table",       1, NA,  4, "max", "metadata_qc_table()", "QC table workers"),
    s("collection",     1, NA,  1, "max", "collection / megafst builders", "collection workers (memory checked by the builders)"),
    s("umap",           1, NA,  1, "max", "collection_umap", "UMAP threads"),
    s("blast",          1, NA,  1, "max", "annotation_build", "BLAST threads"),
    s("metadata",       0, NA,  1, "max", "metadata_bpparam()", "NCBI metadata download workers (rate limited, no CPU check)")
  )
}

resource_step <- function(step) {
  reg <- resource_steps()
  row <- reg[reg$step == step, , drop = FALSE]
  if (!nrow(row)) row <- reg[reg$step == "default", , drop = FALSE]
  as.list(row)
}

# ---- host state ----------------------------------------------------------------------------------------------------

#' Whether thread counts are load-based (config$dynamic_thread_count, else option, default TRUE)
#' @noRd
resource_dynamic_enabled <- function(config = NULL) {
  v <- if (is.list(config)) config[["dynamic_thread_count", exact = TRUE]] else NULL
  if (is.null(v)) v <- getOption("massiveNGSpipe.dynamic_thread_count", TRUE)
  !isFALSE(v)
}

#' Memory (GB) the thread decisions keep free (option massiveNGSpipe.mem_reserve_gb, default 80)
#' @noRd
resource_reserve_gb <- function() as.numeric(getOption("massiveNGSpipe.mem_reserve_gb", 80))

#' Host resources now: cores, load average, memory, swap and the resource mode
#'
#' Reads /proc (cheap); fields are NA where /proc is unavailable.
#' @param mode the resource mode (default: \code{resource_mode_read()}), NULL when none is active
#' @return list(cores, load1, load5, load15, mem_total_gb, mem_available_gb, swap_total_gb, swap_free_gb,
#'  reserve_gb, mode)
#' @export
resource_snapshot <- function(mode = resource_mode_read()) {
  cores <- tryCatch(max(1L, parallel::detectCores(logical = TRUE), na.rm = TRUE), error = function(e) NA_integer_)
  load <- tryCatch(as.numeric(strsplit(readLines("/proc/loadavg", n = 1L), " ")[[1]][1:3]),
                   error = function(e) rep(NA_real_, 3))
  mi <- tryCatch(readLines("/proc/meminfo"), error = function(e) character())
  kb <- function(key) {
    m <- grep(paste0("^", key, ":"), mi, value = TRUE)
    if (length(m) != 1) NA_real_ else as.numeric(gsub("[^0-9]", "", m)) / 1024^2
  }
  list(cores = cores, load1 = load[1], load5 = load[2], load15 = load[3],
       mem_total_gb = kb("MemTotal"), mem_available_gb = kb("MemAvailable"),
       swap_total_gb = kb("SwapTotal"), swap_free_gb = kb("SwapFree"),
       reserve_gb = resource_reserve_gb(), mode = mode)
}

# ---- the formula ---------------------------------------------------------------------------------------------------

#' The one thread formula (pure; see resource_threads())
#'
#' load-based = min(max, max(min, floor((cores - load1) / cores_per_worker)), floor((mem - reserve) / per_worker_gb), n);
#' the mode caps it further at floor(max * factor), max_threads and floor(max_mem_gb / per_worker_gb).
#' Never above \code{max}, never below 1.
#' @noRd
resource_threads_formula <- function(max, min = 1, cores = NA, load1 = NA, mem_gb = NA, reserve_gb = 80,
                                     cores_per_worker = 1, per_worker_gb = NA, n = Inf, mode = NULL,
                                     fallback = max) {
  min <- base::min(min, max)
  cpu_ok <- !is.na(cores_per_worker) && cores_per_worker > 0
  mem_ok <- !is.na(per_worker_gb) && per_worker_gb > 0
  if ((cpu_ok && (is.na(cores) || is.na(load1))) || (mem_ok && is.na(mem_gb))) {
    x <- fallback                                             # no /proc: the step's fallback
  } else {
    x <- max
    if (cpu_ok) x <- base::min(x, base::max(min, floor((cores - load1) / cores_per_worker)))
    if (mem_ok) x <- base::min(x, floor((mem_gb - reserve_gb) / per_worker_gb))
  }
  x <- base::min(x, n)
  if (!is.null(mode)) {                                       # master regulator: caps only
    x <- base::min(x, base::max(1, floor(max * mode$factor)))
    if (!is.na(mode$max_threads)) x <- base::min(x, mode$max_threads)
    if (!is.na(mode$max_mem_gb) && mem_ok) x <- base::min(x, floor(mode$max_mem_gb / per_worker_gb))
  }
  as.integer(base::max(1L, base::min(x, max)))
}

#' Thread / worker count for one step
#'
#' The single entry point for every thread and worker decision in the package. \code{max} is the configured value
#' (config$threads[[step]], a function default, an option). With dynamic thread counts on (default) the result is
#' load-based within 1..max (\code{resource_threads_formula()}), capped further by an active resource mode
#' (\code{resource_mode_read()}); with \code{config$dynamic_thread_count = FALSE} (or option
#' \code{massiveNGSpipe.dynamic_thread_count = FALSE}) it is \code{max}, unchanged.
#' @param step step name, see \code{resource_steps()} (unknown names use the "default" row)
#' @param max the configured maximum
#' @param per_worker_gb,min,cores_per_worker override the registry row
#' @param n another upper bound (e.g. number of files)
#' @param config pipeline config (for dynamic_thread_count), or NULL for the option
#' @param snapshot host state, default \code{resource_snapshot()} (argument for tests)
#' @return integer (or \code{max} itself when dynamic counts are off or max is not one positive number)
#' @export
resource_threads <- function(step, max, per_worker_gb = NULL, min = NULL, n = Inf, config = NULL,
                             cores_per_worker = NULL, snapshot = NULL) {
  if (!resource_dynamic_enabled(config)) return(max)
  if (is.null(max) || length(max) != 1 || !is.numeric(max) || is.na(max) || max < 1) return(max)
  r <- resource_step(step)
  if (is.null(snapshot)) snapshot <- resource_snapshot()
  fb <- switch(r$fallback, "1" = 1L, dt = base::min(16L, data.table::getDTthreads()), max)
  resource_threads_formula(max = max, min = min %||% r$min, cores = snapshot$cores, load1 = snapshot$load1,
                           mem_gb = snapshot$mem_available_gb, reserve_gb = snapshot$reserve_gb %||% resource_reserve_gb(),
                           cores_per_worker = cores_per_worker %||% r$cores_per_worker,
                           per_worker_gb = per_worker_gb %||% r$per_worker_gb, n = n, mode = snapshot$mode,
                           fallback = base::min(fb, max))
}

#' Memory budget (GB) of a step under the resource mode
#'
#' The configured budget, or less when a resource mode is active: budget * mode factor, capped at max_mem_gb.
#' @noRd
resource_mem_gb <- function(max_gb, config = NULL, mode = resource_mode_read()) {
  if (!resource_dynamic_enabled(config) || is.null(mode) || is.null(max_gb)) return(max_gb)
  x <- max_gb * mode$factor
  if (!is.na(mode$max_mem_gb)) x <- base::min(x, mode$max_mem_gb)
  x
}

#' Re-apply the BLAS/OpenMP thread cap of this process (stage start)
#'
#' No-op unless a resource mode is or was active, so without a mode file the value set by pipeline_config() stays.
#' @noRd
resource_apply_blas <- function(config) {
  cap <- config$threads_blas_cap
  if (is.null(cap) || !resource_dynamic_enabled(config)) return(invisible(NULL))
  mode <- resource_mode_read()
  if (is.null(mode) && is.null(.resource_env$blas)) return(invisible(NULL))
  th <- if (is.null(mode)) cap else resource_threads("blas", cap, config = config)
  if (!identical(.resource_env$blas, th)) {
    RhpcBLASctl::blas_set_num_threads(th); RhpcBLASctl::omp_set_num_threads(th)
    .resource_env$blas <- if (is.null(mode)) NULL else th
  }
  invisible(th)
}

# ---- STAR ----------------------------------------------------------------------------------------------------------

#' STAR --outBAMsortingThreadN for one alignment
#'
#' 16 by default (config$star_sort_threads / option massiveNGSpipe.star_sort_threads): -17% wall time on heavy
#' libraries, byte-identical BAM, ~30 GB instead of ~16 GB peak private RAM (star_speed/REVIEW.md, 2026-11-24).
#' Falls back to STAR's own 6 when the memory left after the planned BAM sort is below
#' \code{config$star_sort_extra_gb} (default 14), and never exceeds the mapping threads. NA when STAR does not sort.
#' @param sort_plan from bam_sort_plan()
#' @param run_threads STAR --runThreadN of this alignment
#' @noRd
star_sort_threads <- function(config, sort_plan, run_threads, snapshot = NULL) {
  if (!identical(sort_plan$sorter, "STAR")) return(NA_integer_)
  want <- config$star_sort_threads %||% getOption("massiveNGSpipe.star_sort_threads", 16L)
  th <- resource_threads("star_sort", want, config = config, snapshot = snapshot)
  headroom <- sort_plan$decision$free_gb - sort_plan$decision$needed_gb
  if (length(headroom) != 1 || is.na(headroom) || headroom < (config$star_sort_extra_gb %||% 14)) th <- base::min(th, 6L)
  as.integer(base::max(1L, base::min(th, run_threads)))
}

#' ORFIK_STAR_EXTRA_ARGS value carrying --outBAMsortingThreadN
#'
#' Keeps any user-set extra arguments; a user-set --outBAMsortingThreadN wins. Equal to STAR's default
#' (min(6, runThreadN)) or NA: the existing value unchanged, so the STAR command line stays as before.
#' @noRd
star_extra_args <- function(existing, sort_threads, run_threads) {
  existing <- trimws(existing %||% "")
  if (is.na(sort_threads) || grepl("outBAMsortingThreadN", existing, fixed = TRUE) ||
      sort_threads == base::min(6L, run_threads)) return(existing)
  trimws(paste(existing, "--outBAMsortingThreadN", sort_threads))
}

#' SSD scratch (GB) one STAR alignment needs
#'
#' input copy + STAR temp/BAM + margin. STAR's temp grows with the number of alignments, not the FASTQ size: ~1 GB
#' per million input reads for multi-mapper-heavy Ribo-seq (54 GB for 68 M reads; ~11x the final BAM). Reads are
#' estimated from the file size (massiveNGSpipe.scratch_bytes_per_read: 7 bytes gzipped, 30 plain; small values
#' over-estimate, which only sends a run to the main drive). Was 2 x input + margin (~51 GB for that library).
#' @noRd
align_scratch_need_gb <- function(file1, config, size = file.info(file1)$size) {
  size <- if (is.na(size)) 0 else size
  bpr <- getOption("massiveNGSpipe.scratch_bytes_per_read", c(gz = 7, plain = 30))
  per_read <- if (grepl("\\.gz$", file1)) bpr[["gz"]] else bpr[["plain"]]
  reads_m <- size / per_read / 1e6
  size / 1024^3 + reads_m * (config$scratch_gb_per_million_reads %||% 1) + (config$ssd_min_free_gb %||% 0)
}

# ---- the resource mode file ----------------------------------------------------------------------------------------

resource_mode_levels <- c(normal = 1, conserve = 0.5, minimal = 0.25)

#' Path of the server-wide resource mode file, or NULL when disabled
#'
#' Option \code{massiveNGSpipe.resource_mode_file}, else env \code{MNGSP_RESOURCE_MODE_FILE}, else
#' \code{~/livemount/Bio_data/NGS_pipeline/RESOURCE_MODE} (shared storage of all konsola containers; absent on
#' Valen = no mode). "", "none", "off" or FALSE disable reading.
#' @export
resource_mode_file <- function() {
  p <- getOption("massiveNGSpipe.resource_mode_file", NULL)
  if (is.null(p)) p <- Sys.getenv("MNGSP_RESOURCE_MODE_FILE", "~/livemount/Bio_data/NGS_pipeline/RESOURCE_MODE")
  if (isFALSE(p) || !nzchar(p) || tolower(p) %in% c("none", "off", "false")) return(NULL)
  path.expand(p)
}

resource_mode_max_hours <- function()
  as.numeric(getOption("massiveNGSpipe.resource_mode_max_hours", Sys.getenv("MNGSP_RESOURCE_MODE_MAX_HOURS", "6")))

resource_time_parse <- function(x) {
  if (is.null(x)) return(NA)
  x <- sub("Z$", "", trimws(x))
  t <- as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  if (is.na(t)) t <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  t
}
resource_time_format <- function(t) format(as.POSIXct(t), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

#' Parse and validate the lines of a resource mode file (pure)
#'
#' Valid: \code{mode} normal|conserve|minimal, \code{set_at} and \code{expires_at} (UTC, ISO 8601),
#' set_at <= now + skew (else future-dated: clock skew), expires_at > now, expires_at - set_at <= max_hours + skew;
#' optional \code{max_mem_gb} and \code{max_threads} positive numbers; \code{reason}, \code{set_by}, \code{set_host}
#' free text. Unknown keys are ignored.
#' @return list(ok, why, mode = list(...) or NULL)
#' @noRd
resource_mode_parse <- function(lines, now = Sys.time(), max_hours = resource_mode_max_hours(),
                                skew_min = getOption("massiveNGSpipe.resource_mode_skew_minutes", 10)) {
  bad <- function(why) list(ok = FALSE, why = why, mode = NULL)
  lines <- trimws(lines[!grepl("^\\s*(#|$)", lines)])
  if (!length(lines) || !all(grepl("^[A-Za-z_]+\\s*:", lines))) return(bad("malformed (not 'key: value' lines)"))
  kv <- stats::setNames(trimws(sub("^[^:]*:", "", lines)), tolower(trimws(sub(":.*$", "", lines))))
  get <- function(k) if (k %in% names(kv)) kv[[k]] else NULL
  mode <- get("mode")
  if (is.null(mode) || !(mode %in% names(resource_mode_levels))) return(bad("malformed (mode must be normal, conserve or minimal)"))
  num <- function(k) {
    v <- get(k); if (is.null(v) || !nzchar(v) || tolower(v) %in% c("na", "none")) return(NA_real_)
    x <- suppressWarnings(as.numeric(v)); if (is.na(x) || x <= 0) NaN else x
  }
  max_mem <- num("max_mem_gb"); max_thr <- num("max_threads")
  if (is.nan(max_mem) || is.nan(max_thr)) return(bad("malformed (max_mem_gb / max_threads must be positive numbers)"))
  set_at <- resource_time_parse(get("set_at")); expires <- resource_time_parse(get("expires_at"))
  if (is.na(set_at) || is.na(expires)) return(bad("malformed (set_at and expires_at are required, UTC YYYY-MM-DDTHH:MM:SSZ)"))
  skew <- skew_min * 60
  if (expires <= set_at) return(bad("malformed (expires_at not after set_at)"))
  if (as.numeric(set_at) > as.numeric(now) + skew) return(bad("future-dated set_at (clock skew?)"))
  if (as.numeric(expires) - as.numeric(set_at) > max_hours * 3600 + skew)
    return(bad(sprintf("expires_at more than %g h after set_at", max_hours)))
  if (as.numeric(expires) <= as.numeric(now)) return(bad("expired"))
  list(ok = TRUE, why = "active",
       mode = list(mode = mode, factor = resource_mode_levels[[mode]], max_mem_gb = max_mem,
                   max_threads = if (is.na(max_thr)) NA_real_ else floor(max_thr),
                   reason = get("reason") %||% "", set_by = get("set_by") %||% "", set_host = get("set_host") %||% "",
                   set_at = set_at, expires_at = expires))
}

#' Read the resource mode file: the active mode, or NULL (no file, disabled, invalid, expired)
#'
#' Cheap (one stat + one small read) and never fails: an unreadable, malformed, expired or future-dated file gives one
#' warning per file version and counts as no mode. Read again at every call, so a mode change applies at the next
#' stage start; running work is not touched. With option \code{massiveNGSpipe.resource_mode_rename_expired = TRUE} an
#' expired file is renamed to RESOURCE_MODE.expired (only if unchanged since it was read).
#' @param path file, default \code{resource_mode_file()}
#' @param now time (argument for tests)
#' @param details TRUE returns list(state, why, mode, path) instead
#' @export
resource_mode_read <- function(path = resource_mode_file(), now = Sys.time(), details = FALSE) {
  out <- function(state, why = state, mode = NULL) if (details) list(state = state, why = why, mode = mode, path = path) else mode
  if (is.null(path)) return(out("disabled"))
  info <- tryCatch(file.info(path, extra_cols = FALSE), error = function(e) NULL)
  if (is.null(info) || is.na(info$size)) return(out("absent"))
  key <- paste(path, as.numeric(info$mtime), info$size)
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) e)
  if (inherits(lines, "error")) {
    resource_warn_once(key, paste0("Resource mode file ", path, " unreadable (", conditionMessage(lines), "): ignored"))
    return(out("ignored", "unreadable"))
  }
  p <- resource_mode_parse(lines, now = now)
  if (!p$ok) {
    resource_warn_once(paste(key, p$why), paste0("Resource mode file ", path, " ignored: ", p$why))
    if (p$why == "expired" && isTRUE(getOption("massiveNGSpipe.resource_mode_rename_expired", FALSE))) {
      again <- tryCatch(file.info(path, extra_cols = FALSE), error = function(e) NULL)
      if (!is.null(again) && identical(as.numeric(again$mtime), as.numeric(info$mtime)))
        try(file.rename(path, paste0(path, ".expired")), silent = TRUE)
    }
    return(out("ignored", p$why))
  }
  if (!identical(.resource_env$announced, key)) {
    .resource_env$announced <- key
    message(sprintf("- Resource mode %s (max_threads %s, max_mem_gb %s) by %s until %s: %s", p$mode$mode,
                    p$mode$max_threads, p$mode$max_mem_gb, p$mode$set_by, resource_time_format(p$mode$expires_at),
                    p$mode$reason))
  }
  out("active", "active", p$mode)
}

resource_warn_once <- function(key, msg) {
  if (key %in% .resource_env$warned) return(invisible(FALSE))
  .resource_env$warned <- c(.resource_env$warned, key)
  warning(msg, call. = FALSE)
  invisible(TRUE)
}

#' Set, refresh, clear or show the server-wide resource mode
#'
#' For the master agent (claude-megafst) or an operator. \code{resource_mode_set()} writes the file atomically
#' (temporary file in the same directory + rename, so readers never see half a file; concurrent writers: the last
#' rename wins) and appends a line to RESOURCE_MODE.log. Every file expires (\code{hours} <= 6 by default, option
#' \code{massiveNGSpipe.resource_mode_max_hours}); to keep a mode, the writer must refresh it (heartbeat:
#' \code{resource_mode_refresh()} or set again) before it expires. A crashed writer therefore leaves no lasting
#' restriction.
#' @param mode "normal" (only the max_threads / max_mem_gb caps), "conserve" (x0.5 of each configured maximum) or
#'  "minimal" (x0.25)
#' @param hours hours until the file expires
#' @param reason why (one line)
#' @param max_mem_gb,max_threads optional absolute caps (memory per step, threads per step)
#' @param set_by who set it
#' @param path the mode file
#' @param now time (argument for tests)
#' @return the parsed mode, invisibly
#' @export
resource_mode_set <- function(mode = c("conserve", "minimal", "normal"), hours = 2, reason,
                              max_mem_gb = NA, max_threads = NA,
                              set_by = Sys.getenv("MNGSP_AGENT", Sys.info()[["user"]]),
                              path = resource_mode_file(), now = Sys.time()) {
  mode <- match.arg(mode)
  if (is.null(path)) stop("the resource mode file is disabled (option massiveNGSpipe.resource_mode_file / MNGSP_RESOURCE_MODE_FILE)")
  if (missing(reason) || !nzchar(trimws(reason))) stop("give a reason")
  stopifnot(is.numeric(hours), length(hours) == 1, hours > 0)
  if (hours > resource_mode_max_hours()) stop("hours must be <= ", resource_mode_max_hours(), " (refresh the file to keep the mode)")
  one_line <- function(x) gsub("[\r\n]+", " ", as.character(x))
  lines <- c("# massiveNGSpipe resource mode; see ?resource_mode_set. Expires automatically.",
             paste0("mode: ", mode),
             paste0("max_mem_gb: ", if (is.na(max_mem_gb)) "NA" else max_mem_gb),
             paste0("max_threads: ", if (is.na(max_threads)) "NA" else max_threads),
             paste0("reason: ", one_line(reason)),
             paste0("set_by: ", one_line(set_by)),
             paste0("set_host: ", Sys.info()[["nodename"]]),
             paste0("set_at: ", resource_time_format(now)),
             paste0("expires_at: ", resource_time_format(now + hours * 3600)))
  p <- resource_mode_parse(lines, now = now)
  if (!p$ok) stop("resource mode would be invalid: ", p$why)
  dir <- dirname(path)
  if (!dir.exists(dir)) stop("directory of the resource mode file does not exist: ", dir)
  tmp <- file.path(dir, sprintf(".%s.tmp.%s.%d.%s", basename(path), Sys.info()[["nodename"]], Sys.getpid(),
                                paste(sample(c(letters, 0:9), 8, TRUE), collapse = "")))
  ok <- tryCatch({ writeLines(lines, tmp); file.rename(tmp, path) }, error = function(e) { unlink(tmp); stop(e) })
  if (!isTRUE(ok)) { unlink(tmp); stop("could not write ", path) }
  resource_mode_log(path, sprintf("set %s %gh max_mem_gb=%s max_threads=%s by %s: %s", mode, hours, max_mem_gb,
                                  max_threads, set_by, one_line(reason)), now)
  invisible(p$mode)
}

#' @rdname resource_mode_set
#' @export
resource_mode_refresh <- function(hours = 2, path = resource_mode_file(), now = Sys.time()) {
  m <- resource_mode_read(path, now = now)
  if (is.null(m)) stop("no active resource mode to refresh: set one with resource_mode_set()")
  resource_mode_set(m$mode, hours = hours, reason = m$reason, max_mem_gb = m$max_mem_gb, max_threads = m$max_threads,
                    set_by = m$set_by, path = path, now = now)
}

#' @rdname resource_mode_set
#' @param only_set_by clear only when the file was set by this writer (NULL: any writer)
#' @export
resource_mode_clear <- function(only_set_by = NULL, path = resource_mode_file(), now = Sys.time()) {
  if (is.null(path) || !file.exists(path)) { message("No resource mode file"); return(invisible(FALSE)) }
  if (!is.null(only_set_by)) {
    d <- resource_mode_read(path, now = now, details = TRUE)
    if (!is.null(d$mode) && !identical(d$mode$set_by, only_set_by)) {
      message("Resource mode was set by ", d$mode$set_by, ", not ", only_set_by, ": left in place")
      return(invisible(FALSE))
    }
  }
  if (!file.rename(path, paste0(path, ".cleared"))) stop("could not clear ", path)
  resource_mode_log(path, paste("cleared by", Sys.getenv("MNGSP_AGENT", Sys.info()[["user"]])), now)
  message("Resource mode cleared (go wild: configured values with load-based thread counts)")
  invisible(TRUE)
}

#' @rdname resource_mode_set
#' @export
resource_mode_show <- function(path = resource_mode_file(), now = Sys.time()) {
  d <- resource_mode_read(path, now = now, details = TRUE)
  s <- resource_snapshot(mode = d$mode)
  cat("Resource mode file: ", if (is.null(path)) "(disabled)" else path, "\n", sep = "")
  cat("State:              ", d$state, if (d$state == "ignored") paste0(" (", d$why, ")"), "\n", sep = "")
  if (!is.null(d$mode)) {
    m <- d$mode
    cat(sprintf("Mode:               %s (x%g of configured maxima), max_threads %s, max_mem_gb %s\n",
                m$mode, m$factor, m$max_threads, m$max_mem_gb))
    cat(sprintf("Set by:             %s on %s at %s: %s\n", m$set_by, m$set_host, resource_time_format(m$set_at), m$reason))
    cat(sprintf("Expires:            %s (in %.0f min)\n", resource_time_format(m$expires_at),
                as.numeric(difftime(m$expires_at, now, units = "mins"))))
  } else cat("Mode:               none (configured maxima, load-based thread counts)\n")
  cat(sprintf("Host:               %s cores, load %.1f, MemAvailable %.0f GB (reserve %g), swap free %.0f/%.0f GB\n",
              s$cores, s$load1, s$mem_available_gb, s$reserve_gb, s$swap_free_gb, s$swap_total_gb))
  invisible(d)
}

resource_mode_log <- function(path, text, now = Sys.time())
  try(cat(resource_time_format(now), " ", text, "\n", sep = "", file = paste0(path, ".log"), append = TRUE), silent = TRUE)
