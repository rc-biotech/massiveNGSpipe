#' Fast single-pass trimming and collapsing (fastp 1.3.7 compatible)
#'
#' Replaces the fastp trimming pass plus \code{ORFik::collapse.fastq()} for single-end reads with one pass of a
#' small C program (\code{inst/trimcollapse/trimcollapse.c}, compiled on first use). It reproduces fastp 1.3.7's
#' per-read processing exactly (front/tail trimming, polyG auto-trimming for two-colour instruments, adapter
#' trimming including fastp's one-insertion/one-deletion fallback, adapter-dimer and length filters, optional
#' quality/N filter) and writes a fastp-compatible JSON report with the same keys and numbers that massiveNGSpipe
#' and ORFik read. Output files keep the pipeline's names and locations:
#' \itemize{
#'  \item \code{<output.dir>/trim/report_<id>.json}
#'  \item \code{<output.dir>/trim/SINGLE/collapsed_trimmed_<id>.fasta.gz} (ORFik collapsed fasta)
#'  \item or, with \code{fastq_out}, a FASTQ of trimmed reads (barcode detection sample)
#' }
#' \code{id} is the raw file name without its fastq extension, as in fastp/ORFik. Decompression uses ISA-L when
#' the library is available (else zlib); results do not depend on the decompressor or thread count.
#' @param file raw single-end fastq(.gz)
#' @param output.dir directory holding (or receiving) the \code{trim} folder
#' @param adapter adapter sequence or \code{"disable"} (fastp \code{"auto"} detection is not reimplemented)
#' @param trim.front,trim.tail fastp \code{--trim_front1} / \code{--trim_tail1}
#' @param tail_after_adapter integer, 3' cut after adapter removal (as \code{fastp_post_adapter_tail()})
#' @param min.length fastp \code{--length_required}
#' @param quality_filtering logical, fastp's default quality/N filter (the pipeline's main pass disables it)
#' @param reads_to_process numeric or NULL, process only the first reads (detection sample)
#' @param fastq_out character or NULL, write trimmed reads here instead of collapsing
#' @param json character, report path
#' @param threads integer, worker threads
#' @return path of the collapsed fasta (or of \code{fastq_out}), invisibly
#' @export
trimcollapse_run <- function(file, output.dir, adapter, trim.front = 0, trim.tail = 0, tail_after_adapter = 0,
                             min.length = 20, quality_filtering = FALSE, reads_to_process = NULL, fastq_out = NULL,
                             json = NULL, threads = min(16L, parallel::detectCores())) {
  file <- path.expand(file); output.dir <- path.expand(output.dir)
  stopifnot(length(file) == 1, file.exists(file), is.character(adapter), length(adapter) == 1, adapter != "auto")
  id <- trimcollapse_id(file)
  trim_dir <- file.path(output.dir, "trim")
  dir.create(file.path(trim_dir, "SINGLE"), showWarnings = FALSE, recursive = TRUE)
  if (is.null(json)) json <- file.path(trim_dir, paste0("report_", id, ".json"))
  out <- if (is.null(fastq_out)) file.path(trim_dir, "SINGLE", paste0("collapsed_trimmed_", id, ".fasta.gz")) else path.expand(fastq_out)
  args <- c("-i", shQuote(file), if (is.null(fastq_out)) c("-o", shQuote(out)) else c("--fastq-out", shQuote(out)),
            "-j", shQuote(json), "-a", adapter, "-f", trim.front, "-t", trim.tail, "--tail-after-adapter", tail_after_adapter,
            "-l", min.length, "-w", max(1L, as.integer(threads)),
            if (!quality_filtering) "-Q",
            if (!is.null(reads_to_process)) c("-r", format(reads_to_process, scientific = FALSE)))
  message("trimcollapse (fastp-compatible trim + collapse): ", id)
  log <- system2(trimcollapse_binary(), args, stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status")
  if (!is.null(status) && status != 0) stop("trimcollapse failed for ", file, ":\n", paste(log, collapse = "\n"))
  if (!file.exists(out) || !file.exists(json)) stop("trimcollapse produced no output for ", file)
  if (is.null(fastq_out)) fast_trim_state$collapsed <- out
  invisible(out)
}

#' Can the fast trim + collapse path replace fastp for this input?
#'
#' Single-end fastq with an explicit adapter (or \code{"disable"}), when \code{options(massiveNGSpipe.fast_trim)}
#' is not FALSE (set from \code{pipeline_config(fast_trim)} inside the trim stage) and the binary builds.
#' Paired-end data, fasta input and fastp's \code{"auto"} adapter detection keep using fastp.
#' @noRd
fast_trim_enabled <- function(file, file2 = NULL, adapter = "disable") {
  if (isFALSE(getOption("massiveNGSpipe.fast_trim", TRUE))) return(FALSE)
  paired <- !is.null(file2) && length(file2) > 0 && !all(is.na(file2)) && any(nzchar(file2))
  if (paired || length(file) != 1 || !file.exists(file)) return(FALSE)   # missing input: let fastp report it
  if (!grepl("\\.(fastq|fq)(\\.gz)?$", file)) return(FALSE)
  if (identical(adapter, "auto") || !nzchar(adapter) || !grepl("^([ACGT]+|disable)$", adapter)) return(FALSE)
  ok <- tryCatch({ trimcollapse_binary(); TRUE }, error = function(e) { message("fast trim unavailable: ", conditionMessage(e)); FALSE })
  ok
}

fast_trim_state <- new.env(parent = emptyenv())

#' Run id used in trim output names (raw file name without its fastq extension)
#' @noRd
trimcollapse_id <- function(file) tools::file_path_sans_ext(sub("\\.gz$", "", basename(file)))

#' Path to the compiled trimcollapse program, compiling it on first use
#'
#' Tries ISA-L (\code{libisal}, header shipped in \code{inst/trimcollapse/isal}) first and falls back to zlib.
#' Cached per source version in the user cache directory.
#' @return character, path to the executable
#' @export
trimcollapse_binary <- function() {
  src <- system.file("trimcollapse", "trimcollapse.c", package = "massiveNGSpipe")
  if (!nzchar(src)) src <- file.path(find.package("massiveNGSpipe"), "inst", "trimcollapse", "trimcollapse.c")
  if (!file.exists(src)) stop("trimcollapse.c not found in the massiveNGSpipe installation")
  key <- substr(tools::md5sum(src), 1, 10)
  dir <- tools::R_user_dir("massiveNGSpipe", "cache")
  bin <- file.path(dir, paste0("trimcollapse_", key))
  if (file.exists(bin)) return(bin)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  cc <- Sys.getenv("CC", "cc")
  tmp <- paste0(bin, ".tmp", Sys.getpid())
  isal_inc <- file.path(dirname(src), "isal")
  out <- suppressWarnings(system2(cc, c("-O3", "-DHAVE_ISAL", paste0("-I", shQuote(isal_inc)), "-o", shQuote(tmp), shQuote(src),
                                        "-lz", "-lpthread", "-l:libisal.so.2"), stdout = TRUE, stderr = TRUE))
  if (!file.exists(tmp)) {
    out2 <- suppressWarnings(system2(cc, c("-O3", "-o", shQuote(tmp), shQuote(src), "-lz", "-lpthread"), stdout = TRUE, stderr = TRUE))
    if (!file.exists(tmp)) stop("Compiling trimcollapse failed:\n", paste(c(out, out2), collapse = "\n"))
    message("trimcollapse built with zlib (ISA-L not available; install libisal for ~1.5x faster decompression)")
  }
  file.rename(tmp, bin)
  bin
}
