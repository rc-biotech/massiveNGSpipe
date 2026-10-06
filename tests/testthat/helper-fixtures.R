# Shared fixture builders for massiveNGSpipe unit tests.
#
# These build the minimal-but-structurally-real config/pipelines shapes
# that most massiveNGSpipe orchestration functions expect, backed by real
# tempdir() directories (not mocks) so flag/file-based functions behave
# exactly as they do in production, just against throwaway paths.

#' A tiny, real .ofst file with known per-length score-weighted read counts
#'
#' Builds a real GAlignments + ORFik::export.ofst() round trip (not a
#' mock) -- cheap enough for a fast unit test, but exercises the real
#' ofst read/write path read_length_distribution() depends on.
#' @param widths integer vector, one width per row
#' @param scores numeric vector, same length as widths (each row's true
#' read-multiplicity weight)
#' @return character, path to the written .ofst file
fake_ofst <- function(widths = c(28, 28, 30, 30), scores = c(100, 1, 5, 5)) {
  stopifnot(length(widths) == length(scores))
  # GAlignments() itself rejects zero-length seqnames/cigar/strand inputs
  # ("not parallel to x") -- always build >= 1 row, then subset down to
  # the requested (possibly empty) length, which GAlignments handles fine.
  n <- max(length(widths), 1)
  w <- if (length(widths) == 0) 1 else widths
  aln <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(rep("chr1", n)),
    pos = rep(1L, n),
    cigar = paste0(w, "M"),
    strand = S4Vectors::Rle(BiocGenerics::strand(rep("+", n)))
  )
  S4Vectors::mcols(aln)$score <- if (length(scores) == 0) 1 else scores
  if (length(widths) == 0) aln <- aln[0]
  path <- tempfile(fileext = ".ofst")
  ORFik::export.ofst(aln, path)
  path
}

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
    # Matches pipeline_config()'s own default (file.path(project_dir, "FINAL_LIST.csv")).
    complete_metadata = file.path(project, "FINAL_LIST.csv"),
    session_dir = session_dir,
    error_dir = NULL,
    discord_webhook = NULL,
    # FALSE by default -- safe/conservative for a test (never deletes
    # whatever fixture files a test created), matching real
    # pipeline_config()'s own default of mode == "online" rather than
    # hardcoding TRUE.
    delete_raw_files = FALSE,
    delete_trimmed_files = FALSE,
    delete_collapsed_files = FALSE,
    stop_downloading_new_data_at_drive_usage = 92,
    compress_raw_data = TRUE,
    min_raw_reads_pshift = 1e5,
    min_alignment_rate_pshift = 10,
    max_no_adapter_removed_pct = 80,
    skip_pshift_on_low_reads = FALSE,
    skip_pshift_on_low_alignment = FALSE,
    # SerialParam by default -- fast and deterministic for a unit test;
    # bpparam_from_config() needs config$threads to have an entry for
    # every step it's asked about, config$thread_type to be a function,
    # and config$parallel_conf to have exactly these 4 names (see its
    # own stopifnot()s, R/pipeline_init_helpers.R).
    threads = list(main = 1, default = 1, trim = 1, collapse = 1,
                  pshifted = 1, valid_pshift = 1, pcounts = 1),
    thread_type = BiocParallel::SerialParam,
    parallel_conf = BiocParallel::bpoptions(log = FALSE, logdir = NA_character_,
                                            jobname = "test", stop.on.error = TRUE)
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

#' Minimal S4 stand-in for an ORFik experiment, supporting only what
#' convert_per_sample() and its pipeline_convert_*() callers
#' (R/pipeline_preset_steps_sub.R) actually call: name(), runIDs(),
#' single-row `[` subsetting, and uniqueMappers<-. Avoids constructing a
#' real ORFik experiment (needs real bam/reference files) just to test
#' massiveNGSpipe's own per-sample resume orchestration.
methods::setClass("fake_exp_stub", representation(name = "character", run_ids = "character",
                                                   unique_mappers = "logical"),
                  prototype(unique_mappers = FALSE))
methods::setMethod("name", "fake_exp_stub", function(x) x@name)
methods::setMethod("runIDs", "fake_exp_stub", function(x) x@run_ids)
methods::setMethod("[", "fake_exp_stub", function(x, i, ...) {
  methods::new("fake_exp_stub", name = x@name, run_ids = x@run_ids[i],
              unique_mappers = x@unique_mappers)
})
methods::setMethod("uniqueMappers<-", "fake_exp_stub", function(x, value) {
  x@unique_mappers <- value
  x
})

fake_experiment_stub <- function(run_ids = c("SRR001", "SRR002", "SRR003"),
                                 exp_name = "PRJNA000001-homo_sapiens") {
  methods::new("fake_exp_stub", name = exp_name, run_ids = run_ids)
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
