#' Build the frame counts by Run and read length for an organism collection
#'
#' Writes QC_STATS/frames_per_library_length.fst under the collection result
#' folder. Rows are sorted by Run, mapper and read_length. The three frame
#' columns count P-sites in coding exons (score weighted; transcripts counted
#' separately, as in the study frame QC). Offset is the signed ORFik
#' offsets_start value.\cr
#' No alignment file is read when the study frame QC exists. Post-shift
#' counts come, in order, from the per-sample shift QC cache
#' (QC_STATS/per_sample_qc/v2_complete_TIS_windows/<mapper>/<Run>_frames.csv,
#' keyed by Run), then from the study table QC_STATS/Ribo_frames_all.csv
#' (labels mapped to Runs through the experiment CSV), and only when neither
#' is usable from the P-shifted ofst. Both QC files must be at least as new
#' as the P-shifted library.\cr
#' The three pre columns (frames of the read 5' ends before shifting) are the
#' post-shift counts rotated per read length: a P-site in CDS frame f has
#' its 5' end in frame (f + offset) mod 3. Frames are in CDS (spliced)
#' coordinates, so this is exact for reads whose 5' end and P-site both lie
#' in the same annotated CDS, including reads across annotated junctions. It
#' is an approximation for reads whose 5' end lies outside the CDS (they are
#' counted although the unshifted 5' end would not hit the CDS, and vice
#' versa) and for reads spliced differently from the annotation.
#' @param collection Collection experiment name (e.g. all_samples-Homo_sapiens)
#'   or an ORFik experiment object for that collection.
#' @param force Recompute every library with available inputs.
#' @param BPPARAM BiocParallel parameter for processing studies.
#' @return The output data.table, invisibly.
#' @export
collection_frames_update <- function(collection, force = FALSE,
                                     BPPARAM = BiocParallel::SerialParam()) {
  df <- if (is.character(collection)) ORFik::read.experiment(collection, validate = FALSE,
                                                           output.env = new.env()) else collection
  stopifnot(all(df$libtype == "RFP"))
  runs <- as.character(ORFik::runIDs(df))
  if (anyNA(runs) || any(!nzchar(runs) | basename(runs) != runs) || anyDuplicated(runs))
    stop("Collection frame QC requires one row per nonempty Run")
  qc_dir <- file.path(ORFik::resFolder(df), "QC_STATS")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(qc_dir, ".frames_per_library_length.lock"), timeout = 0)
  if (is.null(lock)) stop("Another frame build holds the lock for ", qc_dir)
  on.exit(filelock::unlock(lock), add = TRUE)

  build_started <- Sys.time()
  aligned <- ORFik::libFolder(df, mode = "all", unique_mappers = FALSE)
  libraries <- data.table::data.table(Run = runs,
    library_key = if ("library_key" %in% colnames(df)) as.character(df$library_key) else NA_character_,
    experiment = basename(dirname(aligned)), aligned = aligned,
    stem = tools::file_path_sans_ext(basename(df$filepath)))
  groups <- split(libraries, libraries$aligned)
  groups <- groups[order(vapply(groups, nrow, integer(1)), decreasing = TRUE)]
  # Set the fst thread count in the parent. Calling threads_fst() after
  # forking can deadlock if the parent has already used OpenMP.
  fst_threads <- fst::threads_fst(1L)
  on.exit(fst::threads_fst(fst_threads), add = TRUE)
  exons <- qc_cds_exon_table(ORFik::loadRegion(df, part = "cds"))
  annotation <- digest::digest(as.list(exons))
  cache_dir <- file.path(qc_dir, "frames_per_library_length_cache")
  dir.create(cache_dir, showWarnings = FALSE)
  exp_dir <- ORFik::config()["exp"]
  txdb <- df@txdb
  # One failing study is reported and skipped; it never stops the collection build.
  results <- BiocParallel::bplapply(groups, function(libs) {
    tryCatch(collection_frames_study(libs, exons, annotation, txdb, cache_dir, exp_dir, force),
             error = function(e) structure(list(), failed = paste0(libs$experiment[1L], ": ", conditionMessage(e))))
  }, BPPARAM = BPPARAM)
  failed <- unlist(lapply(results, attr, "failed"))
  if (length(failed))
    warning(length(failed), " study(ies) skipped after errors: ", paste(failed, collapse = "; "), call. = FALSE)

  out <- data.table::rbindlist(c(list(cbind(collection_frames_empty(), source = character(),
                                              cached = logical())),
                                 unlist(results, recursive = FALSE)))
  sources <- unique(out[, .(Run, mapper, source, cached)])
  if (nrow(sources)) {
    tab <- sources[, .N, by = .(mapper, source)]
    message("- Frame sources (Run x mapper): ",
            paste0(tab$mapper, "/", tab$source, " = ", tab$N, collapse = ", "),
            "; ", sum(sources$cached), " from the collection cache")
    fallback <- sources[source == "pshifted_ofst", uniqueN(Run)]
    if (fallback) message("- ", fallback, " Run(s) had no usable stored frame QC and were counted from the P-shifted ofst")
  }
  out[, c("source", "cached") := NULL]
  data.table::setkey(out, Run, mapper, read_length)
  # The collection lock does not coordinate with the study P-shift writer.
  # Keep the previous report if a study published new shift inputs mid-build.
  branches <- c(libraries$aligned, file.path(libraries$aligned, "unique_mappers"))
  shifted <- file.path(branches, "pshifted", paste0(libraries$stem, "_pshifted.ofst"))
  shift_tables <- file.path(unique(branches), "pshifted", "shifting_table.rds")
  inputs <- c(shifted, shift_tables)
  changed <- which(file.mtime(inputs) >= build_started)
  if (length(changed))
    stop("P-shift inputs changed during the collection frame build; retry after P-shifting finishes: ",
         paste(head(inputs[changed], 5L), collapse = ", "), call. = FALSE)

  missing <- setdiff(runs, out$Run)
  if (length(missing)) {
    warning(length(missing), " Run(s) have no frame rows (missing shift inputs or no accepted lengths): ",
            paste(head(missing, 10L), collapse = ", "))
  }
  path <- file.path(qc_dir, "frames_per_library_length.fst")
  tmp <- tempfile(tmpdir = qc_dir)
  on.exit(unlink(tmp), add = TRUE)
  fst::write_fst(out, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish collection frame QC: ", path)
  message("- Wrote ", path, " (", nrow(out), " rows, ", data.table::uniqueN(out$Run), " Runs)")
  invisible(out)
}

collection_frames_empty <- function() {
  data.table::data.table(Run = character(), library_key = character(), experiment = character(),
    read_length = integer(), frame0 = numeric(), frame1 = numeric(), frame2 = numeric(),
    pre0 = numeric(), pre1 = numeric(), pre2 = numeric(), offset = integer(), mapper = character())
}

# Version 2: post-shift counts from stored QC, pre-shift counts by rotation (no raw ofst read).
COLLECTION_FRAMES_VERSION <- 2L

collection_frames_study <- function(libs, exons, annotation, txdb, cache_dir, exp_dir, force) {
  experiment <- libs$experiment[1L]
  exp_file <- file.path(exp_dir, paste0(experiment, ".csv"))
  folder <- libs$aligned[1L]
  same_txdb <- NULL
  study_txdb_ok <- function() {
    if (is.null(same_txdb)) same_txdb <<- collection_frames_same_txdb(exp_file, txdb)
    same_txdb
  }
  lapply(c("all", "unique"), function(mapper) {
    branch <- if (mapper == "all") folder else file.path(folder, "unique_mappers")
    raw <- file.path(branch, "ofst", paste0(libs$stem, ".ofst"))
    shifted <- file.path(branch, "pshifted", paste0(libs$stem, "_pshifted.ofst"))
    available <- which(file.exists(shifted))
    if (!length(available)) return(NULL)
    table_dir <- file.path(branch, "pshifted")
    table_path <- file.path(table_dir, "shifting_table.rds")
    table_sig <- qc_file_sig(table_path)
    index <- pshift_table_index(table_path)
    if (is.null(index)) {
      warning("Missing shifting table for ", experiment, " / ", mapper, " (", length(available), " Runs)")
      return(NULL)
    }
    # Some legacy tables contain pre-migration absolute paths. Match only an
    # unambiguous OFST basename within this study and mapper's own table.
    keys <- basename(index$keys)
    exact <- match(normalizePath(raw, mustWork = FALSE), index$keys)
    unique_keys <- which(!duplicated(keys) & !duplicated(keys, fromLast = TRUE))
    legacy <- unique_keys[match(basename(raw), keys[unique_keys])]
    hits <- ifelse(is.na(exact), legacy, exact)

    qc_path <- file.path(branch, "QC_STATS", "Ribo_frames_all.csv")
    sample_dir <- file.path(branch, "QC_STATS", "per_sample_qc", "v2_complete_TIS_windows", mapper)
    missing <- available[is.na(hits[available])]
    if (length(missing))
      warning("Missing offset entries for ", experiment, " / ", mapper, ": ",
              paste(libs$Run[missing], collapse = ", "))
    available <- setdiff(available, missing)
    if (!length(available)) return(NULL)

    qc <- NULL
    qc_loaded <- FALSE
    shifted_info <- file.info(shifted)
    raw_info <- file.info(raw)
    mapper_cache <- file.path(cache_dir, experiment, mapper)
    dir.create(mapper_cache, recursive = TRUE, showWarnings = FALSE)
    tables <- lapply(available, function(i) {
      run <- libs$Run[i]
      offsets <- pshift_canonical_offsets(index$tab[[hits[i]]])
      record_paths <- pshift_offset_record_paths(raw[i], shifted[i])
      if (file.exists(record_paths$pending)) {
        warning("P-site offsets are pending for ", experiment, " / ", mapper, " / ", run,
                call. = FALSE)
        return(NULL)
      }
      if (file.exists(record_paths$applied)) {
        record <- readRDS(record_paths$applied)
        if (!identical(record$offsets, offsets) ||
            !pshift_same_file(record$shifted, pshift_file_identity(shifted[i])) ||
            (!is.na(raw_info$size[i]) && !pshift_same_file(record$raw, pshift_file_identity(raw[i])))) {
          warning("P-site offset provenance differs for ", experiment, " / ", mapper, " / ", run,
                  call. = FALSE)
          return(NULL)
        }
      }
      # Every source is derived from the P-shifted file and the offsets, so these key the cache.
      sig <- list(version = COLLECTION_FRAMES_VERSION, annotation = annotation, offsets = offsets,
                  files = shifted[i], size = shifted_info$size[i],
                  mtime = as.numeric(shifted_info$mtime[i]))
      cache <- file.path(mapper_cache, paste0(run, ".rds"))
      old <- if (!force && file.exists(cache)) tryCatch(readRDS(cache), error = function(e) NULL)
      cached <- identical(old$sig, sig)
      if (cached) {
        counts <- old$counts
        source <- old$source
      } else {
        usable <- function(post) !is.null(post) && nrow(post) > 0
        fresh <- function(path) isTRUE(file.mtime(path) >= shifted_info$mtime[i])
        post <- NULL
        source <- "sample_qc"
        sample_qc <- file.path(sample_dir, paste0(run, "_frames.csv"))
        if (fresh(sample_qc) && study_txdb_ok())
          post <- collection_frames_sample_qc(sample_qc, txdb)
        if (!usable(post)) {
          source <- "study_qc"
          if (!qc_loaded) {
            qc <<- collection_frames_study_qc(qc_path, exp_file, txdb)
            qc_loaded <<- TRUE
          }
          post <- if (!is.null(qc) && fresh(qc_path)) qc[[run]]
        }
        if (!usable(post)) {
          source <- "pshifted_ofst"
          post <- qc_frame_sums(qc_read_ofst(shifted[i]), exons)
        }
        counts <- collection_frames_count(post, offsets)
        dropped <- attr(counts, "dropped_lengths")
        if (length(dropped))
          warning(experiment, " / ", mapper, " / ", run, ": ", length(dropped),
                  " read length(s) without a P-site offset dropped (", paste(dropped, collapse = ","),
                  "); the shifting table may not match the P-shifted file", call. = FALSE)
        if (!identical(sig$size, file.info(sig$files)$size) ||
            !identical(sig$mtime, as.numeric(file.mtime(sig$files))))
          stop("P-shifted OFST changed while counting ", experiment, " / ", mapper, " / ", run)
        pshift_atomic_record(list(sig = sig, counts = counts, source = source), cache)
      }
      if (!nrow(counts)) return(NULL)
      counts <- data.table::copy(counts)
      counts[, `:=`(Run = run, library_key = libs$library_key[i],
                    experiment = libs$experiment[i], mapper = mapper)]
      data.table::setcolorder(counts, names(collection_frames_empty()))
      counts[, `:=`(source = source, cached = cached)]
      counts
    })
    if (!identical(table_sig, qc_file_sig(table_path)))
      stop("Shifting table changed while counting ", experiment, " / ", mapper)
    data.table::rbindlist(tables)
  })
}

# The study QC was counted on the study annotation, so it is only reused when
# the study uses the collection's txdb.
collection_frames_same_txdb <- function(exp_file, txdb) {
  tryCatch({
    study <- ORFik::read.experiment(exp_file, validate = FALSE, output.env = new.env())
    identical(normalizePath(study@txdb), normalizePath(txdb))
  }, error = function(e) FALSE)
}

collection_frames_read_qc <- function(path) {
  qc <- data.table::fread(path, select = c("fraction", "length", "frame", "score"))
  if (!nrow(qc) || anyNA(qc) || any(!qc$frame %in% 0:2) ||
      any(!is.finite(qc$score) | qc$score < 0) ||
      anyDuplicated(qc, by = c("fraction", "length", "frame"))) return(NULL)
  qc
}

# Per-sample shift QC cache: one library per file, named by Run, so no label mapping.
collection_frames_sample_qc <- function(path, txdb) {
  if (isTRUE(file.mtime(path) < file.mtime(txdb))) return(NULL)
  tryCatch({
    qc <- collection_frames_read_qc(path)
    if (is.null(qc) || data.table::uniqueN(qc$fraction) != 1L) return(NULL)
    qc[, .(fraction = length, frame, score)]
  }, error = function(e) NULL)
}

# Label variants a study frame table can carry for each library: the current
# shift QC label (bamVarName without experiment and libtype), with libtype,
# and both with "_NA" tokens dropped (older label scheme).
collection_frames_label_keys <- function(study) {
  full <- shift_qc_sample_names(study)
  short <- sub("^RFP_", "", full)
  strip <- function(x) gsub("_NA(?=_|$)", "", x, perl = TRUE)
  list(short, full, strip(short), strip(full))
}

# Reuse study QC only if its labels map uniquely across the entire report.
# A newer experiment CSV can change the meaning of a previously valid label.
collection_frames_study_qc <- function(path, exp_file, txdb) {
  mtimes <- file.mtime(c(path, exp_file, txdb))
  if (anyNA(mtimes) || any(mtimes[1L] < mtimes[-1L])) return(NULL)
  tryCatch({
    qc <- collection_frames_read_qc(path)
    if (is.null(qc)) return(NULL)
    study <- ORFik::read.experiment(exp_file, validate = FALSE, output.env = new.env())
    if (!identical(normalizePath(study@txdb), normalizePath(txdb))) return(NULL)
    runs <- as.character(ORFik::runIDs(study))
    if (anyDuplicated(runs)) return(NULL)
    # Tables with a Run column (shift QC since 2026-11-26) need no label mapping
    run_col <- tryCatch(data.table::fread(path, select = "Run")$Run, error = function(e) NULL)
    if (length(run_col) == nrow(qc) && !anyNA(run_col) && all(nzchar(run_col)) && all(run_col %in% runs)) {
      if (anyDuplicated(data.frame(run_col, qc$length, qc$frame))) return(NULL)
      return(split(qc[, .(fraction = length, frame, score)], run_col))
    }
    labels <- unique(as.character(qc$fraction))
    mapped <- setNames(rep(NA_character_, length(labels)), labels)
    for (key in collection_frames_label_keys(study)) {
      unique_key <- !duplicated(key) & !duplicated(key, fromLast = TRUE)
      todo <- is.na(mapped)
      mapped[todo] <- runs[unique_key][match(labels[todo], key[unique_key])]
    }
    if (anyNA(mapped) || any(!nzchar(mapped)) || anyDuplicated(mapped)) return(NULL)
    split(qc[, .(fraction = length, frame, score)], unname(mapped[as.character(qc$fraction)]))
  }, error = function(e) NULL)
}

# post: data.table(fraction = read length, frame, score) of P-sites in CDS frames.
# Pre-shift frames by rotation: a P-site in frame f has its 5' end in frame
# (f + offset) mod 3 (offset = offsets_start <= 0). See collection_frames_update().
collection_frames_count <- function(post, offsets) {
  # Lengths without an offset cannot be rotated (and have no output row); they are dropped and returned
  # as attribute dropped_lengths for the caller's per-library warning.
  dropped <- sort(unique(as.integer(post$fraction[!post$fraction %in% offsets$fraction])))
  post <- data.table::as.data.table(post)[fraction %in% offsets$fraction]
  post <- post[, .(score = sum(as.numeric(score))), by = .(read_length = as.integer(fraction), frame = as.integer(frame))]
  post[, offset := offsets$offsets_start[match(read_length, offsets$fraction)]]
  post[, pre := (frame + offset) %% 3L]
  # Retain every length with an estimated shift, including those with no CDS hits.
  counts <- data.table::as.data.table(offsets)
  data.table::setnames(counts, c("read_length", "offset"))
  for (f in 0:2) {
    after <- post[frame == f]
    before <- post[pre == f, .(score = sum(score)), by = read_length]
    counts[, (paste0("frame", f)) := data.table::fcoalesce(after$score[match(read_length, after$read_length)], 0)]
    counts[, (paste0("pre", f)) := data.table::fcoalesce(before$score[match(read_length, before$read_length)], 0)]
  }
  if (length(dropped)) data.table::setattr(counts, "dropped_lengths", dropped)
  counts
}
