# Incremental, resource-checked megafst / megacovrle builders.
#
# The megafst of an organism collection (all_samples-<Organism>) stores every library's P-site
# coverage as one integer column per Run in fixed genomic pages (one fst file per page and strand,
# collection_tables_indexed/coverage_<chr>_part<N>_{forward,reverse}.fst) plus coverage_index.fst.
# Page layout, file names, column values and the index are the same as the original
# fst_megaformat_make() output; what changes is how pages are produced:
#
# * Libraries are streamed in batches straight from their per-library covRle files; the page
#   writer does not need the whole-collection megacovrle in memory.
# * A manifest (Run, covRle path, mtime, size) records what each page set contains. On rerun only
#   new libraries are appended (fst fork append_columns_fst, no rewrite of existing columns) and
#   libraries whose covRle changed are replaced (replace_existing_columns). Removed libraries,
#   changed layout or accumulated replacement garbage trigger a full rebuild.
# * Full rebuilds are written to a staging directory, validated, then each file is renamed into
#   place (index last), preserving unrelated files in the live directory.
# * Every load is preceded by a RAM check that waits for memory and finally errors.
# * Appends are idempotent: a page's existing column names are the source of truth, so an
#   interrupted build resumes where it stopped.
# Benchmarks and design notes: ~/shared_workspace/hakon_playground/megafst.

MEGAFST_FORMAT_VERSION <- 1L

#' Linux MemAvailable in GB (GiB), NA where unavailable
#' @return numeric
#' @noRd
megafst_memory_available_gb <- function() {
  if (!file.exists("/proc/meminfo")) return(NA_real_)
  line <- grep("^MemAvailable:", readLines("/proc/meminfo", warn = FALSE), value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line)) / 1024^2
}

#' Wait until a load of needed_gb fits, keeping keep_free_gb for everyone else; error on timeout
#'
#' Polls MemAvailable every poll_seconds. Errors after timeout_minutes instead of loading into a
#' machine that cannot hold it (an OOM kill would also hit unrelated jobs).
#' @param needed_gb numeric, estimated GB the next load will allocate
#' @param what character, label for messages
#' @param keep_free_gb numeric, memory that must remain available after the load
#' @param timeout_minutes numeric, give up after this long
#' @param poll_seconds numeric
#' @param available_fun function returning available GB (injectable for tests)
#' @return invisible(TRUE)
#' @noRd
wait_for_memory <- function(needed_gb, what, keep_free_gb = 32, timeout_minutes = 120,
                            poll_seconds = 60, available_fun = megafst_memory_available_gb) {
  start <- Sys.time()
  repeat {
    avail <- available_fun()
    if (is.na(avail)) {
      warning("Available memory unknown on this system; loading ", what, " without a RAM check")
      return(invisible(TRUE))
    }
    if (avail - needed_gb >= keep_free_gb) return(invisible(TRUE))
    waited <- as.numeric(difftime(Sys.time(), start, units = "mins"))
    if (waited >= timeout_minutes) {
      stop(sprintf(paste("Not enough memory to load %s: needs ~%.1f GB and %.0f GB must stay free,",
                         "but only %.1f GB is available after waiting %.0f minutes."),
                   what, needed_gb, keep_free_gb, avail, waited), call. = FALSE)
    }
    message(sprintf("-- Waiting for memory (%s): need %.1f GB + %.0f GB reserve, %.1f GB available (%.0f/%.0f min)",
                    what, needed_gb, keep_free_gb, avail, waited, timeout_minutes))
    Sys.sleep(poll_seconds)
  }
}

#' Per-library input signature: Run, covRle path, mtime and size
#' @param runs character
#' @param cov_paths character, same order as runs
#' @return data.table
#' @noRd
megafst_library_table <- function(runs, cov_paths) {
  stopifnot(length(runs) == length(cov_paths), !anyDuplicated(runs), all(file.exists(cov_paths)))
  info <- file.info(cov_paths)
  data.table(Run = runs, cov_path = normalizePath(cov_paths), cov_mtime = as.numeric(info$mtime),
             cov_size = as.numeric(info$size))
}

#' Page layout of a collection, identical to the original fst_megaformat_make() order and names
#' @param seqinfo Seqinfo of the collection genome
#' @param tile_size integer, page length in nt
#' @param output_dir character, directory the index will point to
#' @param bp2 integer, legacy chromosome interleave factor (kept for identical page order)
#' @return data.table chr, part, start, end, file_forward, file_reverse
#' @noRd
megafst_page_table <- function(seqinfo, tile_size, output_dir, bp2 = 5) {
  stopifnot(tile_size >= 1, tile_size <= .Machine$integer.max)
  seqinfo <- seqinfo[GenomeInfoDb::seqnames(seqinfo)[order(GenomeInfoDb::seqlengths(seqinfo), decreasing = TRUE)]]
  seqinfo <- seqinfo[GenomeInfoDb::seqnames(seqinfo)[order(rep(seq(bp2), length.out = length(GenomeInfoDb::seqlengths(seqinfo))))]]
  lens <- GenomeInfoDb::seqlengths(seqinfo)
  stopifnot(!anyNA(lens))
  rbindlist(lapply(names(lens), function(chr) {
    n_tiles <- ceiling(lens[[chr]] / tile_size)
    part <- seq_len(n_tiles)
    start <- as.integer((part - 1) * tile_size + 1)   # integer, as in the original IRanges-based index
    end <- as.integer(pmin(part * tile_size, lens[[chr]]))
    base <- file.path(output_dir, paste0("coverage_", chr, "_part", part))
    data.table(chr = chr, part = part, start = start, end = end,
               file_forward = paste0(base, "_forward.fst"), file_reverse = paste0(base, "_reverse.fst"))
  }))
}

#' Decide what an update needs
#' @param old data.table or NULL, previous manifest
#' @param new data.table, current megafst_library_table()
#' @param params_old list or NULL, previous layout parameters
#' @param params_new list, current layout parameters
#' @param compact_fraction numeric, rebuild when replaced columns since the last full build
#' would exceed this fraction of all columns (replacement payloads stay in the files)
#' @return list(mode = "none" | "incremental" | "full", reason, add, replace)
#' @noRd
megafst_plan <- function(old, new, params_old, params_new, compact_fraction = 0.25,
                         max_appends = 64, batch_size = 500) {
  full <- function(reason) list(mode = "full", reason = reason, add = new$Run, replace = character())
  if (is.null(old) || is.null(params_old)) return(full("no manifest"))
  layout_keys <- c("format_version", "tile_size", "compress", "seq_names", "seq_lengths")
  if (!identical(params_old[layout_keys], params_new[layout_keys])) return(full("layout parameters changed"))
  removed <- setdiff(old$Run, new$Run)
  if (length(removed)) return(full(paste(length(removed), "libraries removed")))
  m <- merge(new, old, by = "Run", all.x = TRUE, suffixes = c("", "_old"), sort = FALSE)
  add <- m[is.na(cov_path_old)]$Run
  replace <- m[!is.na(cov_path_old) & (cov_path != cov_path_old | cov_mtime != cov_mtime_old |
                                          cov_size != cov_size_old)]$Run
  replaced_total <- (params_old$replaced_since_full %||% 0) + length(replace)
  if (replaced_total > compact_fraction * nrow(new)) {
    return(full(sprintf("%d replaced columns since last full build exceed %.0f%% of %d libraries",
                        replaced_total, 100 * compact_fraction, nrow(new))))
  }
  # Each append/replace pass adds a column segment to every page; > ~80 segments slow cold reads
  # (benchmark B5), so compact with a staged full rebuild before that.
  passes <- ceiling(length(add) / batch_size) + ceiling(length(replace) / batch_size)
  appends_total <- (params_old$appends_since_full %||% 0) + passes
  if (passes > 0 && appends_total > max_appends) {
    return(full(sprintf("%d append passes since last full build exceed %d (segment compaction)",
                        appends_total, max_appends)))
  }
  if (!length(add) && !length(replace)) return(list(mode = "none", reason = "up to date",
                                                   add = character(), replace = character()))
  list(mode = "incremental", reason = sprintf("%d new, %d changed libraries", length(add), length(replace)),
       add = add, replace = replace, passes = passes)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Write one batch of libraries into every page (both strands)
#'
#' op "add": write the page if missing, else append only columns not yet present (idempotent
#' resume). op "replace": replace existing columns in place (fst fork). Each page's values are
#' as.integer(covRle_strand[[chr]][start:end]) exactly as in the original builder.
#' @param covs named list (names = Run) of covRle
#' @param pages data.table from megafst_page_table() with page files in the target directory
#' @param op "add" or "replace"
#' @param compress fst compression
#' @param BPPARAM BiocParallel param for pages of one chromosome (forked workers share covs)
#' @return invisible(NULL)
#' @noRd
megafst_write_batch <- function(covs, pages, op = c("add", "replace"), compress = 50,
                                BPPARAM = BiocParallel::SerialParam()) {
  op <- match.arg(op)
  runs <- names(covs)
  stopifnot(length(runs) == length(covs), !anyDuplicated(runs), !anyNA(runs))
  write_page <- function(i, pages, strand_list, runs, file_col, op, compress) {
    file <- pages[[file_col]][i]
    present <- if (file.exists(file)) fst::metadata_fst(file)$columnNames else character()
    todo <- if (op == "add") setdiff(runs, present) else intersect(runs, present)
    if (op == "replace" && length(todo) != length(runs)) {
      stop("Cannot replace libraries missing from page ", file, ": ",
           paste(setdiff(runs, present), collapse = ", "))
    }
    if (!length(todo)) return(invisible(NULL))
    chr <- pages$chr[i]
    window <- IRanges::IRanges(pages$start[i], pages$end[i])
    dt <- setDT(lapply(strand_list[todo], function(x) as.integer(x[[chr]][window])))
    setnames(dt, todo)
    if (op == "replace") {
      fst::replace_existing_columns(file, dt, compress = compress)
    } else if (!length(present)) {
      tmp <- paste0(file, ".tmp")
      fst::write_fst(dt, tmp, compress = compress)
      stopifnot(file.rename(tmp, file))
    } else {
      fst::append_columns_fst(dt, file, compress = compress)
    }
    invisible(NULL)
  }
  for (strand in c("forward", "reverse")) {
    strand_list <- lapply(covs, if (strand == "forward") f else r)
    file_col <- if (strand == "forward") "file_forward" else "file_reverse"
    # One parallel pass over all pages of this strand: forked workers share strand_list
    # (one fork set per strand per batch, not per chromosome).
    BiocParallel::bplapply(seq_len(nrow(pages)), write_page, pages = pages, strand_list = strand_list,
                           runs = runs, file_col = file_col, op = op, compress = compress,
                           BPPARAM = BPPARAM)
    rm(strand_list); gc(FALSE)
  }
  invisible(NULL)
}

#' Check every page has exactly the expected columns and row count
#' @param pages data.table from megafst_page_table()
#' @param runs character, expected column names
#' @return invisible(TRUE); errors listing offending files
#' @noRd
megafst_validate_pages <- function(pages, runs) {
  files <- c(pages$file_forward, pages$file_reverse)
  rows <- rep(pages$end - pages$start + 1, 2)
  bad <- character()
  for (k in seq_along(files)) {
    if (!file.exists(files[k])) { bad <- c(bad, paste(files[k], "missing")); next }
    m <- fst::metadata_fst(files[k])
    if (m$nrOfRows != rows[k]) bad <- c(bad, paste(files[k], "rows", m$nrOfRows, "!=", rows[k]))
    if (!setequal(m$columnNames, runs) || anyDuplicated(m$columnNames))
      bad <- c(bad, paste(files[k], "columns differ from manifest"))
  }
  if (length(bad)) stop("megafst validation failed:\n", paste(utils::head(bad, 20), collapse = "\n"),
                        if (length(bad) > 20) paste0("\n... and ", length(bad) - 20, " more"))
  invisible(TRUE)
}

#' Stream libraries in RAM-checked batches into the pages
#' @param lib data.table from megafst_library_table() restricted to the libraries to write
#' @param pages page table in the target directory
#' @param op "add" or "replace"
#' @param batch_size integer, libraries loaded at once
#' @param covrle_memory_ratio numeric, in-memory GB per GB of covRle file (measured ~10x)
#' @inheritParams wait_for_memory
#' @return invisible(NULL)
#' @noRd
megafst_stream_libraries <- function(lib, pages, op, batch_size, compress, BPPARAM,
                                     covrle_memory_ratio = 11, keep_free_gb = 32,
                                     timeout_minutes = 120, workers = BiocParallel::bpnworkers(BPPARAM)) {
  if (!nrow(lib)) return(invisible(NULL))
  batches <- split(seq_len(nrow(lib)), ceiling(seq_len(nrow(lib)) / batch_size))
  max_page_rows <- max(pages$end - pages$start + 1)
  for (b in seq_along(batches)) {
    rows <- batches[[b]]
    need <- sum(lib$cov_size[rows]) * covrle_memory_ratio / 1024^3 +   # covRle in memory
      workers * max_page_rows * length(rows) * 4 * 2 / 1024^3        # per-worker dense page (+copy)
    wait_for_memory(need, sprintf("megafst batch %d/%d (%d libraries)", b, length(batches), length(rows)),
                    keep_free_gb, timeout_minutes)
    message(sprintf("--- %s batch %d/%d: %d libraries", op, b, length(batches), length(rows)))
    covs <- lapply(lib$cov_path[rows], ORFik:::read_RDSQS, nthread = 5)
    names(covs) <- lib$Run[rows]
    megafst_write_batch(covs, pages, op, compress, BPPARAM)
    rm(covs); gc(FALSE)
  }
  invisible(NULL)
}

#' Build or incrementally update one collection's megafst
#'
#' @param runs character, Run IDs (column names), in collection order
#' @param cov_paths character, per-library covRle files, same order
#' @param seqinfo Seqinfo of the collection genome
#' @param output_dir character, live collection_tables_indexed directory
#' @param tile_size integer page length (default 2e6; must be < 2^31, keep small for cold reads)
#' @param batch_size integer libraries per loaded batch (large keeps segment count low)
#' @param max_appends integer, staged full rebuild after this many append/replace passes
#' @param compress fst compression
#' @param compact_fraction see megafst_plan()
#' @param force logical, force a full rebuild
#' @param adopt_existing logical, if no manifest exists but the live pages already contain exactly
#' these libraries, record a manifest instead of rebuilding (libraries whose covRle is newer than
#' the pages are then replaced). Avoids a full rebuild when first enabling caching.
#' @inheritParams megafst_stream_libraries
#' @return list(plan, output_dir), invisibly
#' @export
megafst_update <- function(runs, cov_paths, seqinfo, output_dir, tile_size = 2e6, batch_size = 500,
                           compress = 50, compact_fraction = 0.25, max_appends = 64, force = FALSE,
                           adopt_existing = TRUE, covrle_memory_ratio = 11, keep_free_gb = 32,
                           timeout_minutes = 120,
                           BPPARAM = BiocParallel::MulticoreParam(4, exportglobals = FALSE, log = FALSE),
                           bp2 = 5) {
  output_dir <- sub("/+$", "", output_dir)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  lock <- filelock::lock(file.path(dirname(output_dir), paste0(".", basename(output_dir), ".lock")),
                         exclusive = TRUE, timeout = 0)
  if (is.null(lock)) stop("Another megafst build holds the lock for ", output_dir)
  on.exit(filelock::unlock(lock), add = TRUE)

  manifest_path <- file.path(output_dir, "megafst_manifest.fst")
  params_path <- file.path(output_dir, "megafst_params.rds")
  index_path <- file.path(output_dir, "coverage_index.fst")
  lib <- megafst_library_table(runs, cov_paths)
  pages <- megafst_page_table(seqinfo, tile_size, output_dir, bp2)
  params <- list(format_version = MEGAFST_FORMAT_VERSION, tile_size = as.numeric(tile_size),
                 compress = as.numeric(compress), seq_names = unique(pages$chr),
                 seq_lengths = as.numeric(pages[, max(end), by = chr]$V1), replaced_since_full = 0,
                 appends_since_full = 0)
  old <- if (file.exists(manifest_path)) fst::read_fst(manifest_path, as.data.table = TRUE) else NULL
  params_old <- if (file.exists(params_path)) readRDS(params_path) else NULL

  if (is.null(old) && adopt_existing && !force && file.exists(index_path)) {
    adopted <- megafst_try_adopt(lib, pages, index_path, params)
    if (!is.null(adopted)) { old <- adopted$old; params_old <- adopted$params
      message("-- Adopted existing megafst pages as the cache baseline (", nrow(old), " libraries)") }
  }
  plan <- if (force) list(mode = "full", reason = "forced", add = lib$Run, replace = character()) else
    megafst_plan(old, lib, params_old, params, compact_fraction, max_appends, batch_size)
  message("-- megafst plan: ", plan$mode, " (", plan$reason, ")")
  if (plan$mode == "none") return(invisible(list(plan = plan, output_dir = output_dir)))

  if (plan$mode == "incremental") {
    params$replaced_since_full <- (params_old$replaced_since_full %||% 0) + length(plan$replace)
    params$appends_since_full <- (params_old$appends_since_full %||% 0) + plan$passes
    megafst_stream_libraries(lib[Run %in% plan$replace], pages, "replace", batch_size, compress, BPPARAM,
                             covrle_memory_ratio, keep_free_gb, timeout_minutes)
    megafst_stream_libraries(lib[Run %in% plan$add], pages, "add", batch_size, compress, BPPARAM,
                             covrle_memory_ratio, keep_free_gb, timeout_minutes)
    megafst_validate_pages(pages, lib$Run)
    megafst_save_state(lib, params, manifest_path, params_path)
    if (!file.exists(index_path)) fst::write_fst(pages, index_path)
    return(invisible(list(plan = plan, output_dir = output_dir)))
  }

  # Full rebuild in staging (resumable: existing staging pages keep their columns)
  staging <- paste0(output_dir, ".staging")
  dir.create(staging, showWarnings = FALSE)
  staging_pages <- copy(pages)
  staging_pages[, file_forward := file.path(staging, basename(file_forward))]
  staging_pages[, file_reverse := file.path(staging, basename(file_reverse))]
  megafst_stream_libraries(lib, staging_pages, "add", batch_size, compress, BPPARAM,
                           covrle_memory_ratio, keep_free_gb, timeout_minutes)
  megafst_validate_pages(staging_pages, lib$Run)
  # Publish: rename each page into place (atomic per file), drop stale pages, index + state last
  new_files <- c(staging_pages$file_forward, staging_pages$file_reverse)
  stopifnot(all(file.rename(new_files, file.path(output_dir, basename(new_files)))))
  stale <- setdiff(list.files(output_dir, "^coverage_.*_part[0-9]+_(forward|reverse)\\.fst$"),
                   basename(new_files))
  if (length(stale)) unlink(file.path(output_dir, stale))
  tmp_index <- paste0(index_path, ".tmp"); fst::write_fst(pages, tmp_index)
  stopifnot(file.rename(tmp_index, index_path))
  megafst_save_state(lib, params, manifest_path, params_path)
  unlink(staging, recursive = TRUE)
  invisible(list(plan = plan, output_dir = output_dir))
}

megafst_save_state <- function(lib, params, manifest_path, params_path) {
  tmp <- paste0(manifest_path, ".tmp"); fst::write_fst(lib, tmp); stopifnot(file.rename(tmp, manifest_path))
  tmp <- paste0(params_path, ".tmp"); saveRDS(params, tmp); stopifnot(file.rename(tmp, params_path))
  invisible(NULL)
}

#' Use pre-caching pages as the baseline if they match the current libraries and layout
#' @return NULL (cannot adopt) or list(old, params)
#' @noRd
megafst_try_adopt <- function(lib, pages, index_path, params) {
  idx <- tryCatch(fst::read_fst(index_path, as.data.table = TRUE), error = function(e) NULL)
  if (is.null(idx) || nrow(idx) != nrow(pages) ||
      !identical(idx[, .(chr = as.character(chr), start = as.numeric(start), end = as.numeric(end))],
                 pages[, .(chr = as.character(chr), start = as.numeric(start), end = as.numeric(end))])) return(NULL)
  files <- file.path(dirname(index_path), basename(c(pages$file_forward, pages$file_reverse)))
  if (!all(file.exists(files))) return(NULL)
  cols <- fst::metadata_fst(files[1])$columnNames
  if (!setequal(cols, lib$Run)) return(NULL)
  if (!all(vapply(files, function(f) setequal(fst::metadata_fst(f)$columnNames, cols), logical(1)))) return(NULL)
  oldest_page <- min(as.numeric(file.info(files)$mtime))
  old <- copy(lib)
  # A library whose covRle is newer than the oldest page may not be reflected: mark it changed
  old[cov_mtime > oldest_page, cov_mtime := -1]
  list(old = old, params = params)
}

#' Build or incrementally update one collection's per-strand megacovrle
#'
#' Product: megalist_covrle/all_covrles_<Org>_{forward,reverse}.qs, lists (names = Run) of
#' RleList per library, in collection order. Cached with its own manifest: unchanged runs are kept,
#' changed runs re-read, new runs added, removed runs dropped; nothing is rewritten when up to date.
#' Every load is RAM-checked (in-memory size ~10x the qs file size).
#' @inheritParams megafst_update
#' @param file_paths character(2), forward and reverse output files
#' @param collection_csv character, collection experiment CSV (for adopting pre-caching files)
#' @return invisibly the plan mode
#' @export
megacovrle_update <- function(runs, cov_paths, file_paths, force = FALSE, covrle_memory_ratio = 11,
                              keep_free_gb = 32, timeout_minutes = 120, collection_csv = NULL) {
  stopifnot(length(file_paths) == 2)
  dir.create(dirname(file_paths[1]), showWarnings = FALSE, recursive = TRUE)
  manifest_path <- file.path(dirname(file_paths[1]), "megacovrle_manifest.fst")
  lib <- megafst_library_table(runs, cov_paths)
  old <- if (file.exists(manifest_path) && all(file.exists(file_paths)) && !force)
    fst::read_fst(manifest_path, as.data.table = TRUE) else NULL
  # Adopt pre-caching files (unnamed, collection order) without loading them, if they are newer
  # than the collection definition and every covRle: names come from the manifest on next load.
  if (is.null(old) && !force && all(file.exists(file_paths)) && !is.null(collection_csv) &&
      file.exists(collection_csv)) {
    built <- min(as.numeric(file.info(file_paths)$mtime))
    if (built > as.numeric(file.info(collection_csv)$mtime) && built > max(lib$cov_mtime)) {
      tmp <- paste0(manifest_path, ".tmp"); fst::write_fst(lib, tmp); stopifnot(file.rename(tmp, manifest_path))
      message("-- Adopted existing megacovrle as cache baseline (", nrow(lib), " libraries)")
      return(invisible("adopted"))
    }
  }
  if (is.null(old)) {
    mode <- "full"; read_runs <- lib$Run; removed <- character()
  } else {
    m <- merge(lib, old, by = "Run", all.x = TRUE, suffixes = c("", "_old"), sort = FALSE)
    read_runs <- m[is.na(cov_path_old) | cov_path != cov_path_old | cov_mtime != cov_mtime_old |
                     cov_size != cov_size_old]$Run
    removed <- setdiff(old$Run, lib$Run)
    mode <- if (length(read_runs) || length(removed) || !identical(old$Run, lib$Run)) "incremental" else "none"
  }
  message(sprintf("-- megacovrle plan: %s (%d to read, %d removed)", mode, length(read_runs), length(removed)))
  if (mode == "none") return(invisible("none"))
  new_gb <- sum(lib[Run %in% read_runs]$cov_size) * covrle_memory_ratio / 1024^3
  for (s in 1:2) {
    keep_gb <- if (mode == "full") 0 else file.size(file_paths[s]) * covrle_memory_ratio / 1024^3
    wait_for_memory(keep_gb + new_gb, paste("megacovrle", c("forward", "reverse")[s]), keep_free_gb, timeout_minutes)
    strand_list <- if (mode == "full") list() else ORFik:::read_RDSQS(file_paths[s], nthread = 5)
    if (length(strand_list) && is.null(names(strand_list)) && length(strand_list) == nrow(old))
      names(strand_list) <- old$Run   # adopted pre-caching file: manifest order = list order
    if (length(strand_list) && (is.null(names(strand_list)) || !all(names(strand_list) %in% old$Run)))
      stop("Existing megacovrle ", file_paths[s], " does not match its manifest; rerun with force = TRUE")
    strand_list <- strand_list[intersect(names(strand_list), setdiff(lib$Run, read_runs))]
    for (k in which(lib$Run %in% read_runs)) {
      cov <- ORFik:::read_RDSQS(lib$cov_path[k], nthread = 5)
      strand_list[[lib$Run[k]]] <- if (s == 1) f(cov) else r(cov)
    }
    strand_list <- strand_list[lib$Run]
    stopifnot(identical(names(strand_list), lib$Run))
    attr(strand_list, "strand") <- c("+", "-")[s]
    tmp <- paste0(file_paths[s], ".tmp")
    qs2::qs_save(strand_list, tmp, nthreads = 5)
    stopifnot(file.rename(tmp, file_paths[s]))
    rm(strand_list); gc(FALSE)
  }
  tmp <- paste0(manifest_path, ".tmp"); fst::write_fst(lib, tmp); stopifnot(file.rename(tmp, manifest_path))
  invisible(mode)
}
