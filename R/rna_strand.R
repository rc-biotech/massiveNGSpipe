# Final-BAM RNA orientation, distinct from pre-alignment mate selection.
# Cache annotations in the parent process so forked converters share them.
.rna_strand_annotation_cache <- new.env(parent = emptyenv())

rna_strand_signature <- function(path) {
  st <- file.info(path)
  if (!file.exists(path)) stop("RNA strand input missing: ", path)
  c(size = st$size, mtime = as.numeric(st$mtime))
}

rna_strand_regions <- function(df) {
  annotation <- df@txdb
  key <- digest::digest(list(annotation, rna_strand_signature(annotation)), algo = "xxhash64")
  if (exists(key, .rna_strand_annotation_cache, inherits = FALSE))
    return(get(key, .rna_strand_annotation_cache, inherits = FALSE))
  ex <- GenomicRanges::reduce(unlist(ORFik::loadRegion(df, "mrna"), use.names = FALSE))
  ex <- ex[!grepl("^(chr)?(M|MT)$", as.character(seqnames(ex)), ignore.case = TRUE)]
  plus <- ex[as.character(strand(ex)) == "+"]
  minus <- ex[as.character(strand(ex)) == "-"]
  ex <- c(IRanges::subsetByOverlaps(plus, minus, ignore.strand = TRUE, invert = TRUE),
          IRanges::subsetByOverlaps(minus, plus, ignore.strand = TRUE, invert = TRUE))
  value <- list(regions = ex, annotation_signature = key)
  assign(key, value, .rna_strand_annotation_cache)
  value
}

# One vote per informative read, even if its spliced blocks hit several exons.
rna_strand_classify <- function(reads, regions, min_informative = 1000L,
                                stranded_fraction = .8, unstranded_tolerance = .1) {
  same <- GenomicRanges::countOverlaps(reads, regions, ignore.strand = FALSE) > 0L
  any <- GenomicRanges::countOverlaps(reads, regions, ignore.strand = TRUE) > 0L
  # Regions have no opposite-strand overlaps, but a spliced read can bridge
  # different genes. Exclude those conflicting assignments too.
  flipped <- reads
  strand(flipped) <- ifelse(as.character(strand(reads)) == "+", "-", "+")
  opposite <- GenomicRanges::countOverlaps(flipped, regions, ignore.strand = FALSE) > 0L
  keep <- any & xor(same, opposite)
  rna_strand_decide(sum(same & keep), sum(opposite & keep), min_informative,
                    stranded_fraction, unstranded_tolerance)
}

rna_strand_decide <- function(sense, antisense, min_informative,
                              stranded_fraction, unstranded_tolerance) {
  n <- sense + antisense
  fraction <- if (n) sense / n else NA_real_
  status <- "unknown"; mode <- NA_integer_; read_strand <- NA_character_
  if (n >= min_informative) {
    if (fraction >= stranded_fraction) {status <- "stranded"; mode <- 1L; read_strand <- "alignment"}
    else if (fraction <= 1 - stranded_fraction) {status <- "stranded"; mode <- 2L; read_strand <- "opposite"}
    else if (abs(fraction - .5) <= unstranded_tolerance) {status <- "unstranded"; mode <- 0L; read_strand <- "unstranded"}
  }
  list(status = status, strandMode = mode, read_strand = read_strand,
       informative_reads = n, sense_reads = sense, antisense_reads = antisense,
       sense_fraction = fraction)
}

rna_strand_detect <- function(bam, annotation, max_records = 1000000L,
                              min_informative = 1000L, stranded_fraction = .8,
                              unstranded_tolerance = .1, sampling = NULL,
                              sampling_fraction = .1, sampling_seed = 17L) {
  stopifnot(length(max_records) == 1L, is.finite(max_records), max_records >= 1,
    length(min_informative) == 1L, min_informative >= 1,
    stranded_fraction > .5, stranded_fraction <= 1,
    unstranded_tolerance >= 0, unstranded_tolerance < stranded_fraction - .5)
  start <- proc.time()[[3]]
  random <- identical(sampling, "whole_bam_random")
  stopifnot(is.null(sampling) || random)
  input <- bam
  parameters <- list(max_records = max_records, min_informative = min_informative,
    stranded_fraction = stranded_fraction, unstranded_tolerance = unstranded_tolerance)
  if (random) {
    stopifnot(length(sampling_fraction) == 1L, is.finite(sampling_fraction),
      sampling_fraction > 0, sampling_fraction <= 1,
      length(sampling_seed) == 1L, is.finite(sampling_seed), sampling_seed >= 0,
      sampling_seed == as.integer(sampling_seed))
    parameters <- c(parameters, list(sampling = sampling,
      sampling_fraction = sampling_fraction, sampling_seed = sampling_seed))
    if (sampling_fraction < 1) {
      samtools <- Sys.which("samtools")
      if (!nzchar(samtools)) stop("Whole-BAM strand sampling requires samtools")
      input <- tempfile(fileext = ".bam")
      on.exit(unlink(input), add = TRUE)
      fraction <- format(sampling_fraction, scientific = FALSE, digits = 15, trim = TRUE)
      seed_fraction <- paste0(sampling_seed, sub("^0", "", fraction))
      status <- system2(samtools, c("view", "-s", seed_fraction, "-F", "2308",
        "-b", "-o", shQuote(input), shQuote(path.expand(bam))))
      if (status != 0L || !file.exists(input)) stop("Whole-BAM strand subsampling failed")
    }
  }
  bf <- Rsamtools::BamFile(input, yieldSize = as.integer(min(100000L, max_records)))
  open(bf); on.exit(close(bf))
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE), tag = "NH", what = "flag")
  scanned <- unique <- 0L; sense <- antisense <- 0L; paired_seen <- single_seen <- FALSE
  repeat {
    reads <- GenomicAlignments::readGAlignments(bf, param = param)
    if (!length(reads)) break
    flags <- mcols(reads)$flag
    is_pair <- bitwAnd(flags, 1L) != 0L
    paired_seen <- paired_seen || any(is_pair); single_seen <- single_seen || any(!is_pair)
    if (paired_seen && single_seen) stop("Mixed single/paired records in RNA BAM; split the library before orientation inference")
    scanned <- scanned + length(reads)
    if (paired_seen) reads <- reads[bitwAnd(flags, 64L) != 0L]
    if (!"NH" %in% names(mcols(reads))) stop("RNA strand inference requires NH unique-mapping tags")
    reads <- reads[!is.na(mcols(reads)$NH) & mcols(reads)$NH == 1L]
    if (any(!as.character(seqnames(reads)) %in% seqlevels(annotation$regions))) {
      # Harmonise naming style (1 vs chr1); reads on sequences absent from the annotation cannot overlap a
      # region, so they are dropped. (ORFik::matchSeqStyle is not exported and refuses to drop used levels.)
      try(suppressWarnings(GenomeInfoDb::seqlevelsStyle(reads) <- GenomeInfoDb::seqlevelsStyle(annotation$regions)[1]),
          silent = TRUE)
      reads <- reads[as.character(seqnames(reads)) %in% seqlevels(annotation$regions)]
    }
    unique <- unique + length(reads)
    chunk <- rna_strand_classify(reads, annotation$regions, min_informative,
                                stranded_fraction, unstranded_tolerance)
    sense <- sense + chunk$sense_reads; antisense <- antisense + chunk$antisense_reads
    result <- rna_strand_decide(sense, antisense, min_informative,
                               stranded_fraction, unstranded_tolerance)
    # Coordinate-sorted prefixes can look strongly stranded even when the
    # whole library is balanced. Consume the ENTIRE random subset; never
    # accept an early classified prefix of that subset.
    if (scanned >= max_records) {
      if (random && length(GenomicAlignments::readGAlignments(bf, param = param)))
        stop("Whole-BAM strand sample exceeds max_records; lower sampling_fraction or increase max_records")
      break
    }
    if (!random && result$status != "unknown") break
    Rsamtools::yieldSize(bf) <- as.integer(min(100000L, max_records - scanned))
  }
  if (!exists("result", inherits = FALSE)) result <- rna_strand_decide(0L, 0L, min_informative,
    stranded_fraction, unstranded_tolerance)
  c(result, list(version = 1L, source_bam = normalizePath(bam),
    product_signature = rna_strand_signature(bam),
    source_bam_signature = rna_strand_signature(bam),
    annotation_signature = annotation$annotation_signature,
    parameters = parameters,
    scanned_reads = scanned, unique_reads = unique,
    elapsed_seconds = proc.time()[[3]] - start, timestamp = as.character(Sys.time()),
    bam_layout = if (paired_seen) "paired_first_mate" else "single_mate",
    selected_mate = if (paired_seen) 1L else NA_integer_,
    normalized = FALSE))
}

#' Retry an under-informed coordinate-prefix strand scan on a whole-BAM random sample
#'
#' The default scan reads the first \code{max_records} reads of the
#' coordinate-sorted BAM. For deep libraries those can all fall in a gene-poor
#' region (confirmed 2026-10-20, PRJNA765941: 1e6 reads all in chr1:10-632 kb,
#' 0-4 informative reads, while 81 percent of reads map uniquely). Only when
#' that prefix scan is "unknown" because it saw fewer than
#' \code{min_informative} informative reads, the run is re-scanned once on a
#' reproducible whole-BAM random sample sized to ~90 percent of
#' \code{max_records}. Resolved prefix decisions are never changed, and the
#' record keeps the requested parameters (so cached decisions stay fresh) with
#' the sample used in \code{record$sampling_fallback}.
#' @noRd
rna_strand_prefix_fallback <- function(record, bam, annotation, parameters) {
  if (!identical(record$status, "unknown") || !is.null(parameters$sampling) ||
      isTRUE(record$informative_reads >= parameters$min_informative)) return(record)
  total <- rna_strand_mapped_reads(bam)
  if (is.na(total) || total <= parameters$max_records) return(record)
  fraction <- signif(0.9 * parameters$max_records / total, 3)
  message("- RNA strand: prefix scan of ", basename(bam), " saw ", record$informative_reads,
          " informative reads; retrying on a whole-BAM random sample (fraction ", fraction, ")")
  retry <- tryCatch(do.call(rna_strand_detect, c(list(bam = bam, annotation = annotation),
    modifyList(parameters, list(sampling = "whole_bam_random", sampling_fraction = fraction,
                                sampling_seed = 17L)))), error = function(e) {
    message("- RNA strand fallback failed: ", conditionMessage(e)); NULL })
  if (is.null(retry)) return(record)
  retry$sampling_fallback <- retry$parameters[c("sampling", "sampling_fraction", "sampling_seed")]
  retry$parameters <- parameters
  retry
}

#' Mapped primary reads in an indexed BAM (samtools idxstats), NA if unavailable
#' @noRd
rna_strand_mapped_reads <- function(bam) {
  samtools <- Sys.which("samtools")
  if (!nzchar(samtools)) return(NA_real_)
  out <- tryCatch(system2(samtools, c("idxstats", shQuote(path.expand(bam))), stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
  if (!length(out)) return(NA_real_)
  sum(as.numeric(vapply(strsplit(out, "\t"), `[`, character(1), 3)), na.rm = TRUE)
}

rna_strand_ensure <- function(df, annotation, options = list(), require_resolved = TRUE) {
  bam <- unlist(ORFik::filepath(df, "default"), use.names = FALSE)[1]
  if (length(bam) != 1L || !grepl("\\.bam$", bam))
    stop("RNA normalization expects the final single-mate BAM, not paired/derived inputs")
  path <- ORFik:::rna_strand_sidecar(bam)
  parameters <- modifyList(list(max_records = 1000000L, min_informative = 1000L,
    stranded_fraction = .8, unstranded_tolerance = .1), options)
  if (identical(parameters$sampling, "whole_bam_random"))
    parameters <- modifyList(list(sampling_fraction = .1, sampling_seed = 17L), parameters)
  record <- if (file.exists(path)) readRDS(path) else NULL
  fresh <- !is.null(record) && identical(record$version, 1L) &&
    identical(record$product_signature, rna_strand_signature(bam)) &&
    identical(record$annotation_signature, annotation$annotation_signature) &&
    identical(record$parameters[sort(names(record$parameters))], parameters[sort(names(parameters))])
  # A cached under-informed prefix result from before the whole-BAM fallback existed (PRJNA765941, 2026-10-20) is
  # re-evaluated once; records that already went through the fallback carry sampling_fallback and stay cached.
  if (fresh && identical(record$status, "unknown") && is.null(record$sampling_fallback) &&
      is.null(parameters$sampling) && isTRUE(record$informative_reads < parameters$min_informative))
    fresh <- FALSE
  if (!fresh) {
    record <- do.call(rna_strand_detect, c(list(bam = bam, annotation = annotation), parameters))
    record <- rna_strand_prefix_fallback(record, bam, annotation, parameters)
    record$decision_id <- digest::digest(record[c("version", "product_signature",
      "annotation_signature", "parameters", "status", "read_strand")], algo = "xxhash64")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(path)); saveRDS(record, tmp)
    if (!file.rename(tmp, path)) stop("Cannot save RNA orientation: ", path)
  }
  if (require_resolved && record$status == "unknown") stop("RNA orientation unresolved for ", runIDs(df),
    ": ", record$informative_reads, " informative reads, sense fraction ", record$sense_fraction)
  record
}

rna_strand_ofst_current <- function(df, record) {
  path <- unlist(ORFik::filepath(df, "ofst"), use.names = FALSE)[1]
  sidecar <- ORFik:::rna_strand_sidecar(path)
  if (!file.exists(path) || !file.exists(sidecar)) return(FALSE)
  product <- readRDS(sidecar)
  isTRUE(product$normalized) && identical(product$decision_id, record$decision_id) &&
    identical(product$product_signature, rna_strand_signature(path))
}

rna_strand_write_table <- function(df, records) {
  rows <- lapply(seq_len(nrow(df)), function(i) {
    run <- runIDs(df[i, ]); r <- records[[run]]
    if (is.null(r)) return(data.table::data.table(Run = run,
      assay = as.character(df$libtype[i]), status = "not_applicable"))
    data.table::data.table(Run = run, assay = as.character(df$libtype[i]),
      status = r$status, strandMode = r$strandMode, read_strand = r$read_strand,
      informative_reads = r$informative_reads, sense_reads = r$sense_reads,
      antisense_reads = r$antisense_reads, sense_fraction = r$sense_fraction,
      scanned_reads = r$scanned_reads, unique_reads = r$unique_reads,
      detection_seconds = r$elapsed_seconds, decision_id = r$decision_id,
      annotation_signature = r$annotation_signature, timestamp = r$timestamp,
      source_bam = r$source_bam)
  })
  root <- if (length(records)) dirname(records[[1]]$source_bam) else dirname(unlist(filepath(df[1, ], "default"), use.names = FALSE)[1])
  dir <- file.path(root, "QC_STATS"); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(tmpdir = dir)
  data.table::fwrite(data.table::rbindlist(rows, fill = TRUE), tmp)
  if (!file.rename(tmp, file.path(dir, "rna_strand.csv"))) stop("Cannot save study RNA strand table")
}

rna_strand_backup <- function(df, runs) {
  root <- dirname(unlist(filepath(df[match(runs[1], runIDs(df)), ], "default"), use.names = FALSE)[1])
  # A fresh dataset has no derived products to back up. Repairs preserve
  # historical data before flags are reset or conversion overwrites it.
  paths <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  paths <- paths[!grepl("/\\.rna_strand/backups/", paths)]
  selected <- paths[vapply(basename(paths), function(p)
    any(startsWith(p, paste0(runs, ".")) | startsWith(p, paste0(runs, "_"))), logical(1)) &
    !grepl("\\.bam(\\.bai)?$", paths) &
    !grepl("/LOGS/|/raw_data/", paths)]
  selected <- unique(c(selected, paths[grepl("/QC_STATS/countTable_|/pshifted_merged/", paths)]))
  if (!length(selected)) return(invisible(NULL))
  backup_root <- file.path(root, ".rna_strand", "backups")
  dir.create(backup_root, recursive = TRUE, showWarnings = FALSE)
  destination <- tempfile(pattern = paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid(), "_"),
    tmpdir = backup_root)
  for (path in selected) {
    target <- file.path(destination, substring(path, nchar(root) + 2L))
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(path, target, copy.date = TRUE)) stop("Cannot preserve prior strand product: ", path)
  }
  invisible(destination)
}

# Existing siblings remain skipped; regional count matrices are regenerated
# together from fresh and unaffected coverage, preserving every column.
rna_strand_invalidate <- function(df, config, runs) {
  rna_strand_backup(df, runs)
  rna_strand_retire_products(df, runs)
  experiment <- name(df); siblings <- setdiff(runIDs(df), runs)
  for (parent in intersect(c("ofst", "covrle", "bigwig"), names(config$flag))) {
    for (step in c(parent, paste0(parent, "_unique"))) {
      if (step_is_done(config, parent, experiment))
        for (run in setdiff(siblings, samples_done(config, step, experiment)))
          set_sample_flag(config, step, experiment, run)
      for (run in runs) unlink(file.path(sample_flag_dir(config, step, experiment), paste0(run, ".rds")))
    }
    remove_flag(config, parent, experiment)
  }
  for (parent in intersect(c("cigar_collapse", "merged_lib", "pcounts", "QC"), names(config$flag))) remove_flag(config, parent, experiment)
}

# After preservation, remove obsolete active formats so filepath/countTable
# cannot select an old .covrds, split bigWig or count matrix alongside its repair.
rna_strand_retire_products <- function(df, runs) {
  root <- dirname(unlist(filepath(df[match(runs[1], runIDs(df)), ], "default"), use.names = FALSE)[1])
  stems <- unique(c(runs, tools::file_path_sans_ext(basename(df$filepath[runIDs(df) %in% runs]))))
  for (base in c(root, file.path(root, "unique_mappers"))) {
    for (folder in c("ofst", "pshifted", "cov_RLE", "cov_RLE_List", "bigwig")) {
      paths <- list.files(file.path(base, folder), recursive = TRUE, full.names = TRUE)
      remove <- paths[vapply(basename(paths), function(p)
        any(startsWith(p, paste0(stems, ".")) | startsWith(p, paste0(stems, "_"))), logical(1))]
      if (length(remove) && !all(file.remove(remove))) stop("Cannot retire old RNA products")
    }
    matrices <- list.files(file.path(base, "QC_STATS"), pattern = "^countTable_", full.names = TRUE)
    if (length(matrices) && !all(file.remove(matrices))) stop("Cannot retire old RNA count matrices")
  }
  invisible(NULL)
}

# RNA cigar-collapse rewrites OFST in place; carry the decision forward and
# refresh its product signature rather than invalidating normalization.
rna_strand_collapse_ofst <- function(ofst, path) {
  record <- ORFik:::rna_strand_record(path)
  ORFik::export.ofst(collapseDuplicatedReads(ofst), path)
  if (!is.null(record)) {
    record$product_signature <- rna_strand_signature(path)
    saveRDS(record, ORFik:::rna_strand_sidecar(path))
  }
  invisible(NULL)
}

# A merged browser library containing unstranded RNA must be treated as total
# coverage, not a directional transcript track. Native TE counts remain per Run.
rna_strand_tag_merged <- function(df, directory) {
  rna <- toupper(as.character(df$libtype)) %in% c("RNA", "RNA-SEQ", "RNASEQ")
  if (!any(rna)) return(invisible(NULL))
  groups <- bamVarName(df, TRUE, TRUE, TRUE, TRUE, TRUE)
  inputs <- filepath(df, "pshifted")
  for (group in unique(groups)) {
    rows <- which(groups == group & rna)
    if (!length(rows)) next
    records <- lapply(inputs[rows], ORFik:::rna_strand_record)
    if (all(vapply(records, is.null, logical(1)))) next
    if (any(vapply(records, is.null, logical(1))))
      stop("Missing RNA orientation provenance in merged library ", group)
    policies <- vapply(records, function(r) r$read_strand, character(1))
    path <- file.path(directory, paste0(group, ".ofst"))
    record <- list(version = 1L, normalized = TRUE, status = "derived",
      read_strand = if (any(policies == "unstranded")) "unstranded" else "alignment",
      source_decisions = vapply(records, function(r) if (is.null(r$decision_id))
        digest::digest(r, algo = "xxhash64") else r$decision_id, character(1)),
      mixed_strandedness = any(policies == "unstranded") && any(policies != "unstranded"),
      product_signature = rna_strand_signature(path))
    sidecar <- ORFik:::rna_strand_sidecar(path)
    dir.create(dirname(sidecar), recursive = TRUE, showWarnings = FALSE)
    saveRDS(record, sidecar)
  }
  invisible(NULL)
}
