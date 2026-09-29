# Shared fixture builders for massiveNGSpipe unit tests.
#
# These build the minimal-but-structurally-real config/pipelines shapes
# that most massiveNGSpipe orchestration functions expect, backed by real
# tempdir() directories (not mocks) so flag/file-based functions behave
# exactly as they do in production, just against throwaway paths.

#' Minimal real mNGSp config for tests, with real flag directories created
#' under a tempdir() project.
fake_config <- function(project = tempfile("mNGSp_test_"), preset = "RNA-seq",
                        mode = "local", contam = FALSE, session_dir = NULL,
                        extra = list()) {
  flags <- pipeline_flags(project, mode = mode, preset = preset,
                          contam = contam, create_dirs = TRUE)
  flag_steps <- if (preset == "empty") list() else flag_grouping(flags)
  config <- list(
    project = project,
    preset = preset,
    mode = mode,
    flag = flags,
    flag_steps = flag_steps,
    config = c(ref = file.path(project, "ref"), bam = file.path(project, "bam"),
              fastq = file.path(project, "fastq"), exp = file.path(project, "exp")),
    session_dir = session_dir,
    error_dir = NULL,
    discord_webhook = NULL,
    stop_downloading_new_data_at_drive_usage = 92,
    compress_raw_data = TRUE,
    min_raw_reads_pshift = 1e5,
    min_alignment_rate_pshift = 10,
    max_no_adapter_removed_pct = 80,
    skip_pshift_on_low_reads = FALSE,
    skip_pshift_on_low_alignment = FALSE
  )
  # utils::modifyList (not c()) so that names present in `extra` REPLACE
  # the default above rather than being appended as a duplicate --
  # config$x on a list with two same-named elements silently returns
  # only the first one, which would make `extra` overrides no-ops.
  utils::modifyList(config, extra)
}

#' Minimal fake `pipelines` list matching the shape used across the
#' package: named list of accession -> list(accession, organisms, study).
#' @param bam_dir character, default a fresh tempdir(). conf["bam"] for
#' the one organism, used by e.g. qc_verdict_report().
fake_pipelines <- function(accession = "PRJNA000001", organism = "Homo sapiens",
                           runs = data.table::data.table(
                             Run = c("SRR001", "SRR002"),
                             LibraryLayout = c("SINGLE", "SINGLE"),
                             LIBRARYTYPE = c("RFP", "RFP"),
                             ScientificName = organism
                           ),
                           exp_name = paste0(accession, "-", gsub(" ", "_", tolower(organism))),
                           bam_dir = tempfile("fake_bam_")) {
  stats::setNames(list(list(
    accession = accession,
    organisms = stats::setNames(list(list(conf = c(exp = exp_name, bam = bam_dir))), organism),
    study = runs
  )), accession)
}

#' Mark every step of `steps` done for `experiment` in `config` (helper
#' for tests that need a "this experiment is already fully done" state).
fake_mark_all_done <- function(config, steps, experiment) {
  for (s in steps) set_flag(config, s, experiment)
}

#' Synthetic Ribo-seq-like fastq data with a fixed 5'/3' barcode and a
#' REALISTIC, VARYING insert (footprint) length -- adapted from
#' fastqc_adapters_info()'s own roxygen @examples pattern (R/fastq_helpers.R),
#' which already demonstrates building synthetic reads with
#' sample(string_length, ...) varying insert length.
#'
#' Each read is barcode5p + <random insert> + barcode3p + adapter,
#' truncated to raw_len (must be long enough for the full adapter to
#' actually appear post-barcode -- see barcode_detector_single()'s own
#' `max_size_before`/`max_size_after` naming: 75/~30-45 is the user's own
#' worked example), then adapter-trimmed via the package's own
#' remove_adapter_ORFik() to simulate what fastp's own trim step produces.
#' @return list(raw_file, trimmed_file, barcode5p, barcode3p,
#' barcode5p_size, barcode3p_size, insert_lengths, adapter, raw_len, n)
fake_barcode_fixture <- function(n = 2000, insert_lengths = 26:34,
                                 barcode5p = "AATT", barcode3p = "GGCC",
                                 adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                                 raw_len = 75) {
  inserts_len <- sample(insert_lengths, n, replace = TRUE)
  bases <- c("A", "C", "G", "T")
  inserts <- vapply(inserts_len, function(l)
    paste(sample(bases, l, replace = TRUE), collapse = ""), character(1))
  full_reads <- substr(paste0(barcode5p, inserts, barcode3p, adapter), 1, raw_len)
  quals <- strrep("I", nchar(full_reads))

  raw <- Biostrings::DNAStringSet(full_reads)
  names(raw) <- paste0("read", seq_len(n))
  qual_bstrings <- Biostrings::BStringSet(quals)

  raw_file <- tempfile(fileext = ".fastq")
  ShortRead::writeFastq(
    ShortRead::ShortReadQ(sread = raw, quality = ShortRead::FastqQuality(qual_bstrings),
                          id = Biostrings::BStringSet(names(raw))),
    raw_file, compress = FALSE)

  invisible(capture.output(suppressMessages(
    trimmed <- remove_adapter_ORFik(raw, adapter, max.mismatch = 2, add_statistics = FALSE)
  )))
  trimmed_quals <- Biostrings::subseq(qual_bstrings, start = 1, width = Biostrings::width(trimmed))

  trimmed_file <- tempfile(fileext = ".fastq")
  ShortRead::writeFastq(
    ShortRead::ShortReadQ(sread = trimmed, quality = ShortRead::FastqQuality(trimmed_quals),
                          id = Biostrings::BStringSet(names(raw))),
    trimmed_file, compress = FALSE)

  list(raw_file = raw_file, trimmed_file = trimmed_file,
       barcode5p = barcode5p, barcode3p = barcode3p,
       barcode5p_size = nchar(barcode5p), barcode3p_size = nchar(barcode3p),
       insert_lengths = insert_lengths, adapter = adapter, raw_len = raw_len, n = n)
}
