#' Demultiplex pooled Ribo-seq lanes to collapsed per-sample reads
#'
#' Fast replacement for the R implementation (\code{demultiplex_r()}): a small C program
#' (\code{inst/demultiplexer/demultiplexer.c}, compiled on first use) streams each pool file once and
#' assigns every read by the barcode right before the 3' adapter. Read 1 layout:
#' \code{[UMI 5'][insert][UMI 3'][BARCODE][ADAPTER]}. Per read: the first adapter hit (first 12 nt of the
#' adapter, at most \code{adapter_mismatches} mismatches) that is preceded by a valid barcode is used;
#' the barcode must match exactly or, if unambiguous, within \code{barcode_mismatches}; both UMIs are
#' removed and inserts shorter than \code{min_insert} are dropped. Reads of the same \code{final_sample}
#' (e.g. one library sequenced in several lanes or runs) are merged and collapsed into one ORFik
#' collapsed fasta (\code{>seq<i>_x<count>}). Only read 1 is used, also for paired-end lanes.
#'
#' Outputs, ready for the normal pipeline (alignment onwards):
#' \itemize{
#'  \item \code{<out_dir>/SINGLE/collapsed_trimmed_<final_sample>.fasta.gz}
#'  \item \code{<out_dir>/report_<final_sample>.json}: fastp-style read counts (raw = reads carrying the
#'   sample's barcode; trimmed = kept inserts), read by ORFik's trimming QC
#'  \item \code{<path>/demultiplex_stats.csv}: reads per sample and pool file (exact / rescued barcode,
#'   too short), and \code{<path>/demultiplex_pool_stats.csv}: per pool file totals, reads without adapter,
#'   unknown and ambiguous barcodes, and the most frequent unknown barcodes.
#' }
#' The run is resumable: pool files whose per-file statistics exist in \code{tmp_dir} are not demultiplexed
#' again, and collapsed outputs are written atomically.
#' @param path directory of the multiplexed pool files (by convention \code{raw_data/<study>/multiplexed/}).
#' @param controlfile path to the control csv. Required columns: \code{sample}, \code{BARCODE} and either
#'  \code{file} (pool file name in \code{path}, read 1) or the legacy \code{LANE} (file \code{L<LANE>*}).
#'  Optional: \code{final_sample} (merge target; default \code{sample}), \code{umi_5p}, \code{umi_3p},
#'  \code{adapter} (all constant within one pool file; they override the arguments below).
#' @param out_dir output directory (trim directory of the study), default
#'  \code{processed_data/<study>/trim} derived from \code{path}.
#' @param umi_5p_size,umi_3p_size integer, UMI lengths (default 0; control file columns take precedence).
#' @param adapter 3' adapter, default the Illumina TruSeq read 1 adapter.
#' @param min_insert integer, minimum insert length after UMI removal (default 15). The aligner's own
#'  minimum read length applies later.
#' @param barcode_mismatches 0 or 1 (default 1; a 1-mismatch match is used only if unique).
#' @param adapter_mismatches integer, default 2 (over the first 12 adapter nt).
#' @param parallel_files integer, pool files demultiplexed at the same time (each uses pigz with 4 threads
#'  plus one demultiplexing core).
#' @param collapse_workers integer, samples collapsed in parallel.
#' @param tmp_dir directory for per-sample insert files (fast local disk; about 25 bytes per kept read).
#'  Deleted after a successful run unless \code{keep_tmp = TRUE}.
#' @param keep_tmp logical, default FALSE.
#' @return data.table of per-sample statistics (invisibly); files written as described.
#' @export
#' @examples
#' \dontrun{
#' path <- "~/livemount/Bio_data/raw_data/japan_yeast_stress/multiplexed/"
#' demultiplex(path, umi_5p_size = 2, umi_3p_size = 5)
#' }
demultiplex <- function(path, controlfile = file.path(path, "controlfile.csv"),
                        out_dir = file.path(sub("/raw_data/", "/processed_data/", dirname(normalizePath(path))), "trim"),
                        umi_5p_size = 0, umi_3p_size = 0, adapter = "AGATCGGAAGAGC", min_insert = 15,
                        barcode_mismatches = 1, adapter_mismatches = 2, parallel_files = 4, collapse_workers = 8,
                        tmp_dir = file.path(if (dir.exists(path.expand("~/livemount/ssd/tmp"))) "~/livemount/ssd/tmp" else tempdir(),
                                            paste0("demultiplex_", basename(dirname(normalizePath(path))))),
                        keep_tmp = FALSE) {
  path <- normalizePath(path, mustWork = TRUE)
  ctrl <- demultiplex_read_controlfile(controlfile, path, umi_5p_size, umi_3p_size, adapter)
  bin <- demultiplexer_binary()
  tmp_dir <- path.expand(tmp_dir)
  out_dir <- path.expand(out_dir)
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "SINGLE"), showWarnings = FALSE, recursive = TRUE)

  # 1. one streaming pass per pool file (resumable per file)
  files <- unique(ctrl[, .(file, tag, umi_5p, umi_3p, adapter)])
  message("Demultiplexing ", nrow(files), " pool file(s), ", uniqueN(ctrl$final_sample), " final sample(s) -> ", out_dir)
  run_one <- function(i) {
    f <- files[i]
    stats_file <- file.path(tmp_dir, paste0("stats.", f$tag, ".tsv"))
    if (file.exists(stats_file)) return(paste("reused", f$tag))
    bc <- ctrl[file == f$file, .(sample_tag = paste0(sample, "__", f$tag), BARCODE)]
    bc_file <- file.path(tmp_dir, paste0("barcodes.", f$tag, ".tsv"))
    fwrite(bc, bc_file, sep = "\t", col.names = FALSE)
    log <- file.path(tmp_dir, paste0("demux.", f$tag, ".log"))
    args <- c("demux", "-i", shQuote(file.path(path, f$file)), "-b", shQuote(bc_file), "-o", shQuote(tmp_dir), "-t", paste0(f$tag, ".part"),
              "-a", f$adapter, "-u5", f$umi_5p, "-u3", f$umi_3p, "-m", min_insert,
              "-bm", barcode_mismatches, "-am", adapter_mismatches)
    status <- system2(bin, args, stdout = log, stderr = log)
    if (status != 0) stop("demultiplexer failed for ", f$file, " (see ", log, ")")
    # publish this file's outputs only when it completed: rename part files, then the stats file last
    for (s in bc$sample_tag) file.rename(file.path(tmp_dir, paste0(s, ".", f$tag, ".part.seq")),
                                         file.path(tmp_dir, paste0(s, ".seq")))
    file.rename(file.path(tmp_dir, paste0("stats.", f$tag, ".part.tsv")), stats_file)
    paste("done", f$tag)
  }
  res <- parallel::mclapply(seq_len(nrow(files)), run_one, mc.cores = max(1, min(parallel_files, nrow(files))),
                            mc.preschedule = FALSE)
  errs <- vapply(res, function(x) inherits(x, "try-error"), logical(1))
  if (any(errs)) stop("Demultiplexing failed: ", paste(unlist(res[errs]), collapse = "; "))
  message(paste(unlist(res), collapse = ", "))

  # 2. statistics
  st <- rbindlist(lapply(files$tag, function(tg) {
    x <- fread(file.path(tmp_dir, paste0("stats.", tg, ".tsv")), sep = "\t", colClasses = list(character = c("sample", "barcode")))
    x[, tag := tg]
  }))
  per_sample <- st[key == "sample"]
  per_sample[, c("sample", "tag2") := tstrsplit(sample, "__", fixed = TRUE)]
  per_sample <- merge(per_sample, ctrl[, .(sample, tag, file, final_sample)], by.x = c("sample", "tag"), by.y = c("sample", "tag"))
  per_sample[, `:=`(reads_kept = reads_exact + reads_rescued, reads_barcode = reads_exact + reads_rescued + reads_short)]
  fwrite(per_sample[, .(final_sample, sample, file, barcode, reads_barcode, reads_kept, reads_exact, reads_rescued, reads_short,
                        mean_insert = round(insert_len_sum / pmax(reads_kept, 1), 2))],
         file.path(path, "demultiplex_stats.csv"))
  pool_stats <- st[key != "sample", .(file = files$file[match(tag, files$tag)], key, barcode, reads = reads_exact)]
  fwrite(pool_stats, file.path(path, "demultiplex_pool_stats.csv"))

  # 3. merge + collapse per final sample (atomic outputs), fastp-style reports
  targets <- split(per_sample, per_sample$final_sample)
  collapse_one <- function(fs) {
    x <- targets[[fs]]
    out <- file.path(out_dir, "SINGLE", paste0("collapsed_trimmed_", fs, ".fasta.gz"))
    inputs <- file.path(tmp_dir, paste0(x$sample, "__", x$tag, ".seq"))
    status <- system2(bin, c("collapse", "-o", shQuote(out), "-m", min_insert, shQuote(inputs)), stdout = FALSE, stderr = FALSE)
    if (status != 0) stop("collapse failed for ", fs)
    demultiplex_write_report(file.path(out_dir, paste0("report_", fs, ".json")), x, ctrl[final_sample == fs]$adapter[1])
    fs
  }
  res <- parallel::mclapply(names(targets), collapse_one, mc.cores = max(1, collapse_workers), mc.preschedule = FALSE)
  errs <- vapply(res, function(x) inherits(x, "try-error"), logical(1))
  if (any(errs)) stop("Collapsing failed: ", paste(unlist(res[errs]), collapse = "; "))

  tot <- pool_stats[key %in% c("total", "kept"), .(reads = sum(reads)), by = key]
  message(sprintf("Demultiplexing done: %s of %s reads kept (%.1f%%), %d samples. Stats: %s",
                  format(tot[key == "kept"]$reads, big.mark = ","), format(tot[key == "total"]$reads, big.mark = ","),
                  100 * tot[key == "kept"]$reads / max(tot[key == "total"]$reads, 1), length(targets),
                  file.path(path, "demultiplex_stats.csv")))
  if (!keep_tmp) unlink(tmp_dir, recursive = TRUE)
  invisible(per_sample)
}

#' Read and validate a demultiplex control file
#' @return data.table with columns sample, BARCODE, file, tag, final_sample, umi_5p, umi_3p, adapter
#' @noRd
demultiplex_read_controlfile <- function(controlfile, path, umi_5p_size = 0, umi_3p_size = 0, adapter = "AGATCGGAAGAGC") {
  stopifnot(file.exists(controlfile))
  ctrl <- fread(controlfile, colClasses = "character")
  if (!all(c("sample", "BARCODE") %in% colnames(ctrl))) stop("controlfile needs columns sample and BARCODE")
  if (!"file" %in% colnames(ctrl)) {
    if (!"LANE" %in% colnames(ctrl)) stop("controlfile needs a 'file' column (pool file in ", path, ") or 'LANE'")
    lane_files <- list.files(path, "\\.(fastq|fq)(\\.gz)?$")
    ctrl[, file := vapply(LANE, function(l) {
      hit <- grep(paste0("^L", l, "([._]R?1)?[._]"), lane_files, value = TRUE)
      hit <- hit[!grepl("[._]R?2[._]", hit)]
      if (length(hit) != 1) stop("Lane ", l, ": expected one read 1 file named L", l, "*, found: ", paste(hit, collapse = ", "))
      hit
    }, character(1))]
  }
  missing <- setdiff(unique(ctrl$file), list.files(path))
  if (length(missing)) stop("Pool files not found in ", path, ": ", paste(missing, collapse = ", "))
  if (!"final_sample" %in% colnames(ctrl)) ctrl[, final_sample := sample]
  ctrl[is.na(final_sample) | final_sample == "", final_sample := sample]
  if (!"umi_5p" %in% colnames(ctrl)) ctrl[, umi_5p := as.character(umi_5p_size)]
  if (!"umi_3p" %in% colnames(ctrl)) ctrl[, umi_3p := as.character(umi_3p_size)]
  default_adapter <- adapter
  if (!"adapter" %in% colnames(ctrl)) ctrl[, adapter := default_adapter]
  ctrl[is.na(adapter) | adapter == "", adapter := default_adapter]
  ctrl[, BARCODE := toupper(trimws(BARCODE))]
  if (any(!grepl("^[ACGTN]+$", ctrl$BARCODE))) stop("Barcodes must be A/C/G/T/N")
  if (uniqueN(nchar(ctrl$BARCODE)) != 1) stop("All barcodes must have the same length")
  if (any(grepl("__|[/\t ]", ctrl$sample))) stop("sample names must not contain '__', '/', spaces or tabs")
  if (any(grepl("[/\t ]", ctrl$final_sample))) stop("final_sample names must not contain '/', spaces or tabs")
  dup <- ctrl[, .N, by = .(file, BARCODE)][N > 1]
  if (nrow(dup)) stop("Duplicated barcode within a pool file: ", paste(dup$file, dup$BARCODE, collapse = "; "))
  dup <- ctrl[, .N, by = .(file, sample)][N > 1]
  if (nrow(dup)) stop("Duplicated sample within a pool file: ", paste(dup$file, dup$sample, collapse = "; "))
  per_file <- ctrl[, .(n = uniqueN(paste(umi_5p, umi_3p, adapter))), by = file][n > 1]
  if (nrow(per_file)) stop("umi_5p/umi_3p/adapter must be constant within a pool file: ", paste(per_file$file, collapse = ", "))
  ctrl[, tag := make.names(sub("\\.(fastq|fq)(\\.gz)?$", "", file), unique = FALSE)]
  if (uniqueN(ctrl$tag) != uniqueN(ctrl$file)) stop("Pool file names must stay distinct without their fastq extension")
  ctrl[]
}

#' fastp-style report for a demultiplexed sample (read by ORFik's trimming QC)
#' @noRd
demultiplex_write_report <- function(file, x, adapter) {
  raw <- sum(x$reads_barcode); kept <- sum(x$reads_kept)
  rep <- list(summary = list(
    before_filtering = list(total_reads = raw, read1_mean_length = round(sum(x$raw_len_sum) / max(raw, 1), 2)),
    after_filtering = list(total_reads = kept, read1_mean_length = round(sum(x$insert_len_sum) / max(kept, 1), 2))),
    adapter_cutting = list(read1_adapter_sequence = adapter),
    demultiplex = list(source = "massiveNGSpipe::demultiplex", pool_files = x$file, samples = x$sample,
                       reads_exact = sum(x$reads_exact), reads_rescued_1_mismatch = sum(x$reads_rescued),
                       reads_short_insert = sum(x$reads_short)))
  tmp <- paste0(file, ".tmp")
  jsonlite::write_json(rep, tmp, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  file.rename(tmp, file)
}

#' Path to the compiled demultiplexer, compiling it on first use
#'
#' The C source ships in \code{inst/demultiplexer/}; the binary is cached per source version in the user
#' cache directory. Needs a C compiler and zlib headers.
#' @return character, path to the executable
#' @export
demultiplexer_binary <- function() {
  src <- system.file("demultiplexer", "demultiplexer.c", package = "massiveNGSpipe")
  if (!nzchar(src)) {  # devtools::load_all() from a checkout
    src <- file.path(find.package("massiveNGSpipe"), "inst", "demultiplexer", "demultiplexer.c")
  }
  if (!file.exists(src)) stop("demultiplexer.c not found in the massiveNGSpipe installation")
  key <- substr(tools::md5sum(src), 1, 10)
  dir <- tools::R_user_dir("massiveNGSpipe", "cache")
  bin <- file.path(dir, paste0("demultiplexer_", key))
  if (!file.exists(bin)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    cc <- Sys.getenv("CC", "cc")
    tmp <- paste0(bin, ".tmp", Sys.getpid())
    out <- system2(cc, c("-O3", "-o", tmp, shQuote(src), "-lz"), stdout = TRUE, stderr = TRUE)
    if (!file.exists(tmp)) stop("Compiling the demultiplexer failed:\n", paste(out, collapse = "\n"))
    file.rename(tmp, bin)
  }
  bin
}
