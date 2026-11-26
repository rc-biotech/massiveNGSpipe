#' Build the frame counts by Run and read length for an organism collection
#'
#' Writes QC_STATS/frames_per_library_length.fst under the collection result
#' folder. Rows are sorted by Run, mapper and read_length. The three frame
#' columns count P-sites in coding exons, and the three pre columns count
#' aligned 5' ends in coding exons, weighted by read score in both cases.
#' Transcripts are counted separately, as in the study frame QC. Offset is
#' the signed ORFik offsets_start value.
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
  results <- BiocParallel::bplapply(groups, function(libs) {
    collection_frames_study(libs, exons, annotation, txdb, cache_dir, exp_dir, force)
  }, BPPARAM = BPPARAM)

  out <- data.table::rbindlist(c(list(collection_frames_empty()), unlist(results, recursive = FALSE)))
  data.table::setkey(out, Run, mapper, read_length)
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

collection_frames_study <- function(libs, exons, annotation, txdb, cache_dir, exp_dir, force) {
  experiment <- libs$experiment[1L]
  exp_file <- file.path(exp_dir, paste0(experiment, ".csv"))
  folder <- libs$aligned[1L]
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
    missing <- available[is.na(hits[available])]
    if (length(missing))
      warning("Missing offset entries for ", experiment, " / ", mapper, ": ",
              paste(libs$Run[missing], collapse = ", "))
    available <- setdiff(available, missing)
    if (!length(available)) return(NULL)

    qc <- NULL
    qc_loaded <- FALSE
    qc_mtime <- file.mtime(qc_path)
    shifted_info <- file.info(shifted)
    raw_info <- file.info(raw)
    mapper_cache <- file.path(cache_dir, experiment, mapper)
    dir.create(mapper_cache, recursive = TRUE, showWarnings = FALSE)
    tables <- lapply(available, function(i) {
      run <- libs$Run[i]
      offsets <- pshift_canonical_offsets(index$tab[[hits[i]]])
      record_paths <- pshift_offset_record_paths(raw[i], shifted[i])
      if (file.exists(record_paths$pending))
        stop("P-site offsets are pending for ", experiment, " / ", mapper, " / ", run)
      if (file.exists(record_paths$applied)) {
        record <- readRDS(record_paths$applied)
        if (!identical(record$offsets, offsets) ||
            !pshift_same_file(record$shifted, pshift_file_identity(shifted[i])) ||
            (!is.na(raw_info$size[i]) && !pshift_same_file(record$raw, pshift_file_identity(raw[i]))))
          stop("P-site offset provenance differs for ", experiment, " / ", mapper, " / ", run)
      }
      sig <- list(version = 1L, annotation = annotation, offsets = offsets,
                  files = c(raw[i], shifted[i]),
                  size = c(raw_info$size[i], shifted_info$size[i]),
                  mtime = as.numeric(c(raw_info$mtime[i], shifted_info$mtime[i])))
      cache <- file.path(mapper_cache, paste0(run, ".rds"))
      old <- if (!force && file.exists(cache)) tryCatch(readRDS(cache), error = function(e) NULL)
      if (identical(old$sig, sig)) {
        counts <- old$counts
      } else {
        if (!qc_loaded) {
          qc <<- collection_frames_study_qc(qc_path, exp_file, txdb)
          qc_loaded <<- TRUE
        }
        post <- if (!is.null(qc) && isTRUE(qc_mtime >= shifted_info$mtime[i]))
          qc[[run]]
        if (!is.null(post) && any(!post$fraction %in% offsets$fraction)) post <- NULL
        counts <- collection_frames_count(raw[i], shifted[i], offsets, exons, post)
        if (!identical(sig$size, file.info(sig$files)$size) ||
            !identical(sig$mtime, as.numeric(file.mtime(sig$files))))
          stop("OFST sources changed while counting ", experiment, " / ", mapper, " / ", run)
        pshift_atomic_record(list(sig = sig, counts = counts), cache)
      }
      if (!nrow(counts)) return(NULL)
      counts[, `:=`(Run = run, library_key = libs$library_key[i],
                    experiment = libs$experiment[i], mapper = mapper)]
      data.table::setcolorder(counts, names(collection_frames_empty()))
      counts
    })
    if (!identical(table_sig, qc_file_sig(table_path)))
      stop("Shifting table changed while counting ", experiment, " / ", mapper)
    data.table::rbindlist(tables)
  })
}

# Reuse study QC only if its labels map uniquely across the entire report.
# A newer experiment CSV can change the meaning of a previously valid label.
collection_frames_study_qc <- function(path, exp_file, txdb) {
  mtimes <- file.mtime(c(path, exp_file, txdb))
  if (anyNA(mtimes) || any(mtimes[1L] < mtimes[-1L])) return(NULL)
  tryCatch({
    qc <- data.table::fread(path, select = c("fraction", "length", "frame", "score"))
    if (!nrow(qc) || anyNA(qc) || any(!qc$frame %in% 0:2) ||
        any(!is.finite(qc$score) | qc$score < 0) ||
        anyDuplicated(qc, by = c("fraction", "length", "frame"))) return(NULL)
    study <- ORFik::read.experiment(exp_file, validate = FALSE, output.env = new.env())
    if (!identical(normalizePath(study@txdb), normalizePath(txdb))) return(NULL)
    labels <- sub("^RFP_", "", shift_qc_sample_names(study))
    runs <- as.character(ORFik::runIDs(study))
    unique_label <- !duplicated(labels) & !duplicated(labels, fromLast = TRUE)
    mapped <- runs[unique_label][match(qc$fraction, labels[unique_label])]
    if (anyNA(mapped) || any(!nzchar(mapped)) || anyDuplicated(runs)) return(NULL)
    split(qc[, .(fraction = length, frame, score)], mapped)
  }, error = function(e) NULL)
}

collection_frames_count <- function(raw, shifted, offsets, exons, post) {
  points <- NULL
  if (is.null(post) || !file.exists(raw)) points <- qc_read_ofst(shifted)
  if (is.null(post)) post <- qc_frame_sums(points, exons)
  if (any(!post$fraction %in% offsets$fraction)) stop("Shifted read lengths lack offsets in ", shifted)
  if (file.exists(raw)) {
    points <- NULL
    reads <- ORFik::import.ofst(raw)
    sizes <- ORFik::readWidths(reads, along.reference = TRUE)
    keep <- sizes %in% offsets$fraction
    reads <- reads[keep]
    strand <- as.character(BiocGenerics::strand(reads))
    score <- S4Vectors::mcols(reads)$score
    pre_points <- data.table::data.table(
      seqnames = as.character(GenomicRanges::seqnames(reads)),
      start = ifelse(strand == "-", BiocGenerics::end(reads), BiocGenerics::start(reads)),
      strand = strand, size = sizes[keep], score = if (is.null(score)) rep(1, length(reads)) else score)
  } else {
    warning("Raw OFST missing; using genomic offset to recover 5' ends: ", raw)
    pre_points <- points
    offsets_by_read <- offsets$offsets_start[match(pre_points$size, offsets$fraction)]
    if (anyNA(offsets_by_read)) stop("Shifted read lengths lack offsets in ", shifted)
    pre_points[, start := start + ifelse(strand == "-", -offsets_by_read, offsets_by_read)]
  }
  pre <- qc_frame_sums(pre_points, exons)
  # Retain every length with an estimated shift, including those with no CDS hits.
  counts <- data.table::as.data.table(offsets)
  data.table::setnames(counts, c("read_length", "offset"))
  for (f in 0:2) {
    after <- post[frame == f]
    before <- pre[frame == f]
    counts[, (paste0("frame", f)) := data.table::fcoalesce(as.numeric(after$score[
      match(read_length, after$fraction)]), 0)]
    counts[, (paste0("pre", f)) := data.table::fcoalesce(as.numeric(before$score[
      match(read_length, before$fraction)]), 0)]
  }
  counts
}
