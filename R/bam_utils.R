bam_flags_to_dt <- function(bam = "~/livemount/Bio_data/processed_data/PRJEB23398-homo_sapiens/aligned/ERR2193146.bam") {
  # Quick sanity:
  # - STAR unique mappers typically have NH==1 and MAPQ==255
  # - multimappers have NH>1 and MAPQ==0 (unless you changed STAR defaults)
  head(df)

  # Helper for flags
  has_bit <- function(flag, bit) bitwAnd(flag, bit) != 0L

  param <- ScanBamParam(
    what = c("qname", "flag", "mapq", "cigar"),  # add cigar here
    tag  = c("NH", "HI", "AS", "nM")             # NH=#loci, HI=hit index, AS=score, NM=#edits/mismatches
  )

  x <- scanBam(bam, param = param)[[1]]

  dt <- data.table(
    qname  = x$qname,
    flag   = x$flag,
    mapq   = x$mapq,
    Alignment_score     = x$tag$AS,
    cigar  = x$cigar,        # CIGAR string
    NH     = x$tag$NH,
    HI     = x$tag$HI,
    NM     = x$tag$nM,       # edit distance (≈ mismatches/indels)
    is_unmapped      = has_bit(x$flag, 0x4),
    is_secondary     = has_bit(x$flag, 0x100),
    is_supplementary = has_bit(x$flag, 0x800)
  )

  dt$is_primary <- with(dt, !is_unmapped & !is_secondary & !is_supplementary)
  return(dt)
}

STAR_junctions_to_dt <- function(file, motif_in_mrna_sense = TRUE,
                                 verbose = TRUE) {
  stopifnot(file.exists(file))
  stopifnot(identical(tools::file_ext(file), "tab"))
  dt <- fread(file, header = FALSE)
  setnames(dt, c(
    "chromosome",               # V1
    "intron_start",             # V2
    "intron_end",               # V3
    "strand",                   # V4
    "motif",                    # V5
    "annotated",               # V6
    "unique_reads",             # V7
    "multi_reads",              # V8
    "max_overhang"              # V9
  ))
  # Strand: 0 = unstranded, 1 = +, 2 = -
  dt[, strand := factor(
    strand,
    levels = c(0, 1, 2),
    labels = c("*", "+", "-"))
  ]
  if (motif_in_mrna_sense) {
    dt[strand == "-", motif := fifelse(motif %in% c(1L, 2L), 1L,   # GT-AG
                                       fifelse(motif %in% c(3L, 4L), 3L,   # GC-AG
                                               fifelse(motif %in% c(5L, 6L), 5L,   # AT-AC
                                                       motif)))                    # keep 0 as non-canonical
    ]
  }


  # Motif: 0 = non-canonical, 1..6 are splice motifs
  dt[, motif := factor(
    motif,
    levels = 0:6,
    labels = c("non-canonical",
               "GT/AG",
               "CT/AC",
               "GC/AG",
               "CT/GC",
               "AT/AC",
               "GT/AT"))
  ]

  if (verbose) message("- There are ", nrow(dt[annotated == 0]), " novel junctions")
  return(dt)
}

STAR_junctions_to_dt_all <- function(log_dir, verbose = TRUE) {
  junction_files <- list.files(log_dir, "_SJ\\.out\\.tab$", full.names = TRUE)
  list <- lapply(junction_files, function(file) {
    if (verbose) message(basename(file))
    STAR_junctions_to_dt(file, verbose = verbose)
  })
  dt <- rbindlist(list, idcol = TRUE)
  if (verbose) {
    message("Motifs total:")
    print(table(dt$motif))
    message("Novel (0) / annotated (1) junctions total:")
    print(table(dt$annotated))
  }
  return(dt)
}

get_qname <- function(bam, yieldSize = 1e6) {
  bf <- BamFile(bam, yieldSize = yieldSize)
  open(bf)
  on.exit(close(bf))

  qnames <- character()

  repeat {
    sb <- scanBam(
      bf,
      param = ScanBamParam(what = "qname")
    )[[1]]

    if (length(sb$qname) == 0L)
      break

    qnames <- c(qnames, sb$qname)
  }

  return(qnames)
}

#' Recreate alignment report for the decollapsed reads
get_expanded_alignment_metrics <- function(fasta_file, bam_file) {
  message(basename(bam_file))
  fasta_headers <- names(readDNAStringSet(fasta_file, use.names = TRUE))
  system.time(bam_headers <- get_qname(bam_file, 1e7))
  stopifnot(all(bam_headers %in% fasta_headers)) # If false, not from the same file
  stopifnot(length(unique(fasta_headers)) == length(fasta_headers))

  dt_bam <- data.table(seq_id = bam_headers)
  dt_bam <- dt_bam[, .N, by = seq_id]

  dt <- data.table::merge.data.table(data.table(seq_id = fasta_headers),
                                     dt_bam, by = "seq_id", sort = FALSE, all = TRUE)
  dt[, scores := as.integer(gsub(".*_x", "", seq_id))] # TODO: support all formats
  dt[is.na(N), N := 0]
  table(dt$N)

  # Note: the two "unique_aligned_reads" columns used to share the exact
  # same name (a copy-paste slip missing "_collapsed" on the second one)
  # -- data.table allows duplicate column names silently, and $-access
  # picked the first (correct, true-read-count) one, so this wasn't
  # causing a wrong ratio below, but it broke CSV export (two identically
  # -headered columns) and any column-name-based access to the collapsed
  # count. Fixed to match the naming convention every other pair here
  # already uses.
  dt_summary <- data.table(total_input_seqs = sum(dt$scores),
                           total_input_seqs_collapsed = nrow(dt),
                           total_aligned_reads = sum(dt[N > 0]$scores),
                           total_aligned_reads_collapsed = nrow(dt[N > 0]),
                           total_alignments = sum(dt$scores*dt$N),
                           total_alignments_collapsed = sum(dt$N),
                           unique_aligned_reads = sum(dt[N == 1]$scores),
                           unique_aligned_reads_collapsed = nrow(dt[N == 1]),
                           multimapping_aligned_reads = sum(dt[N > 1]$scores),
                           multimapping_aligned_reads_collapsed = nrow(dt[N > 1]))

  dt_summary_per_million <- round(dt_summary / 1e6, 2)
  colnames(dt_summary_per_million) <- paste0(colnames(dt_summary_per_million), "(million reads)")
  dt_summary_per_million

  dt_summary_relative <- data.table(total_aligned_reads = dt_summary$total_aligned_reads / dt_summary$total_input_seqs,
                                    unique_alignments = dt_summary$unique_aligned_reads / dt_summary$total_input_seqs,
                                    multimapping_aligned_reads = dt_summary$multimapping_aligned_reads / dt_summary$total_input_seqs)

  dt_summary_relative_percentage <- 100*round(dt_summary_relative, 2)
  colnames(dt_summary_relative_percentage) <- paste0(colnames(dt_summary_relative_percentage), "%")

  # dt_summary_relative's own (pre-"%") column names collide with
  # dt_summary's raw counts -- both have e.g. "total_aligned_reads", one
  # a count and one a 0-1 ratio. Same silent-duplicate-column issue as
  # unique_aligned_reads above. Renamed only after being used to build
  # dt_summary_relative_percentage's names above, so the final "%"
  # columns keep their existing names unchanged.
  colnames(dt_summary_relative) <- paste0(colnames(dt_summary_relative), "_ratio")

  res <- cbind(dt_summary, dt_summary_per_million, dt_summary_relative, dt_summary_relative_percentage)
  return(res)
}

#' Recreate alignment report for the decollapsed reads for whole exp
#' @param df ORFik experiment
#' @param fasta_dir path of collapsed files
#' @param BPPARAM a BPPARAM object
#' @return a data.table, with uncollapsed alignment statistics
#' @export
get_expanded_alignment_metrics_exp <- function(df, fasta_dir = file.path(dirname(libFolder(df)), "trim", "SINGLE"),
                                               BPPARAM = MulticoreParam(8, exportglobals = FALSE)) {
  stopifnot(is(df, "experiment"))
  stopifnot(dir.exists(fasta_dir))
  bam_files <- filepath(df, "default")
  fasta_files <- file.path(fasta_dir, paste0(ORFik:::remove.file_ext(bam_files, TRUE), ".fasta"))
  stopifnot(all(tools::file_ext(bam_files) == "bam"))
  stopifnot(all(file.exists(bam_files)))
  if (!all(file.exists(fasta_files))) {
    fasta_files[!file.exists(fasta_files)] <- paste0(fasta_files[!file.exists(fasta_files)], ".gz")
    stopifnot(all(file.exists(fasta_files)))
  }

  identical(ORFik:::remove.file_ext(bam_files, TRUE), ORFik:::remove.file_ext(fasta_files, TRUE))
  list <- bpmapply(get_expanded_alignment_metrics, fasta_files, bam_files, SIMPLIFY = FALSE,
                   BPPARAM = BPPARAM)
  dt_final <- rbindlist(list, idcol = TRUE)
  dt_final[, `.id` := ORFik:::remove.file_ext(`.id`, TRUE)]
  dt_final[]
  return(dt_final)
}

#' Stacked bar plot of true (uncollapsed) per-sample alignment rate
#'
#' Reuses the already-computed relative/percentage columns from
#' \code{\link{get_expanded_alignment_metrics}} (\code{unique_alignments%},
#' \code{multimapping_aligned_reads%}) rather than recomputing anything,
#' plus the complement of \code{total_aligned_reads%} for the unmapped
#' share.
#' @param dt data.table, as returned by \code{\link{save_expanded_alignment_metrics}}'s
#' internal table build (one row per sample, with a \code{Run} column
#' and the standard \code{get_expanded_alignment_metrics()} columns)
#' @return a ggplot object
#' @noRd
expanded_alignment_metrics_plot <- function(dt) {
  plot_dt <- data.table::data.table(
    Run = rep(dt$Run, 3),
    category = factor(rep(c("Unique", "Multimapping", "Unmapped"), each = nrow(dt)),
                      levels = c("Unmapped", "Multimapping", "Unique")),
    percent = c(dt[["unique_alignments%"]],
               dt[["multimapping_aligned_reads%"]],
               100 - dt[["total_aligned_reads%"]])
  )
  ggplot2::ggplot(plot_dt, ggplot2::aes(x = Run, y = percent, fill = category)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(Unique = "#2c7bb6", Multimapping = "#fdae61",
                                          Unmapped = "#d7191c")) +
    ggplot2::labs(title = "True (uncollapsed) per-read alignment rate",
                 subtitle = "Collapsed STAR reports weight every unique sequence equally, regardless of how many original reads it represents -- this reweights by each sequence's own collapse multiplicity.",
                 x = NULL, y = "% of true input reads", fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Save the uncollapsed (true per-read) STAR alignment metrics table + plot
#'
#' Collapsed-fasta-based alignment (the norm for Ribo-seq in this
#' pipeline) means STAR's own Log.final.out/full_process.csv report
#' alignment rate per UNIQUE collapsed sequence, not per true read: a
#' library with one sequence collapsing 1,000,000 reads and a second,
#' unrelated singleton sequence that fails to map would show as a 50%
#' alignment rate, when in true-read terms it's over 99.9%.
#' \code{\link{get_expanded_alignment_metrics}} reconstructs the true
#' rate from each sequence's own collapse weight (its fasta header's
#' \code{_x<N>} suffix); this saves that reconstruction (table + a
#' stacked-percentage plot) to \code{out_dir}, following the same
#' \code{QCfolder(df)}-style \code{<bam_dir>/QC_STATS/} convention every
#' other QC output in this package already uses.
#'
#' Called from \code{\link{pipeline_cleanup}} as its final step, once
#' BAMs already have their final \code{<Run>.bam} names, so no ORFik
#' experiment object needs to exist yet (\code{pipeline_cleanup()} runs
#' well before the "exp" step creates one). SINGLE-end samples only:
#' collapsed PAIRED-end fasta live in a different location/shape
#' (\code{trim/PAIRED/}, read-1-only -- see \code{pipeline_collapse()}'s
#' own comment) and aren't covered by
#' \code{get_expanded_alignment_metrics()}'s current design; PAIRED
#' samples are skipped here with a message, not silently mishandled.
#'
#' Never throws: a failure anywhere in this (e.g. one sample's fasta/bam
#' mismatch) must never fail \code{pipeline_cleanup()} itself -- every
#' \code{pipe_*()} wrapper already wraps its whole per-study body in its
#' own \code{try()}, and a hard error here would permanently skip the
#' WHOLE study for the rest of the session (see
#' \code{report_failed_pipe()}) over what is only a diagnostics gap.
#' @param bam_dir character, this organism's final \code{<bam>/aligned}
#' dir (already has the renamed \code{<Run>.bam} files by the time this
#' is called)
#' @param collapsed_dir character, where SINGLE-end collapsed fasta live
#' (\code{<bam>/trim/SINGLE}, matching \code{pipeline_collapse()}'s own
#' output location). Resolved per sample via \code{\link{run_files_organizer}}
#' (same prefix-matching logic fastq/fasta discovery already uses
#' elsewhere in this package), not a hardcoded \code{<Run>.fasta.gz} --
#' real collapsed fasta are written as \code{collapsed_trimmed_<Run>.fasta.gz}
#' (verified live against production data; a plain \code{<Run>.fasta.gz}
#' guess never matches anything).
#' @param study_org data.table, this organism's metadata subset (needs
#' \code{Run}, \code{LibraryLayout})
#' @param out_dir character, where to save the table/plot
#' @param BPPARAM a BPPARAM object
#' @return invisible(NULL)
#' @noRd
save_expanded_alignment_metrics <- function(bam_dir, collapsed_dir, study_org,
                                            out_dir = file.path(bam_dir, "QC_STATS"),
                                            BPPARAM = BiocParallel::SerialParam()) {
  result <- try({
    single_runs <- study_org[LibraryLayout != "PAIRED"]
    if (nrow(single_runs) == 0) {
      message("-- No SINGLE-end samples, skipping expanded alignment metrics")
      return(invisible(NULL))
    }
    paired_runs <- study_org[LibraryLayout == "PAIRED"]
    if (nrow(paired_runs) > 0) {
      message("-- Skipping expanded alignment metrics for ", nrow(paired_runs),
             " PAIRED-end sample(s) (not supported by get_expanded_alignment_metrics())")
    }

    bam_files <- file.path(bam_dir, paste0(single_runs$Run, ".bam"))
    # Resolved one sample at a time (not in one run_files_organizer() call
    # across the whole study) so one sample's unmatched/ambiguous fasta
    # can't abort the rest -- run_files_organizer() itself hard-stops via
    # stopifnot(!anyNA(file_vec)) the moment any row fails to resolve.
    fasta_files <- vapply(seq_len(nrow(single_runs)), function(i) {
      resolved <- tryCatch(
        run_files_organizer(single_runs[i], collapsed_dir, format = ".fasta",
                            extra_files_warning = FALSE)[[1]],
        error = function(e) NA_character_)
      if (length(resolved) != 1) NA_character_ else resolved
    }, character(1))
    have_both <- file.exists(bam_files) & !is.na(fasta_files) & file.exists(fasta_files)
    if (!any(have_both)) {
      message("-- No matching bam/collapsed-fasta pairs found, skipping expanded alignment metrics")
      return(invisible(NULL))
    }
    if (!all(have_both)) {
      message("-- Missing bam or collapsed fasta for ", sum(!have_both),
             " sample(s), computing expanded alignment metrics for the rest")
    }

    dt_list <- BiocParallel::bpmapply(get_expanded_alignment_metrics,
                                      fasta_files[have_both], bam_files[have_both],
                                      SIMPLIFY = FALSE, BPPARAM = BPPARAM)
    names(dt_list) <- single_runs$Run[have_both]
    dt <- data.table::rbindlist(dt_list, idcol = "Run")

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(dt, file.path(out_dir, "expanded_alignment_metrics.csv"))

    plot <- expanded_alignment_metrics_plot(dt)
    ggplot2::ggsave(file.path(out_dir, "expanded_alignment_metrics.png"), plot,
                    width = max(7, 0.3 * nrow(dt) + 2), height = 6, dpi = 150,
                    limitsize = FALSE)
  }, silent = TRUE)

  if (is(result, "try-error")) {
    warning("save_expanded_alignment_metrics() failed (QC-only, not fatal): ",
           conditionMessage(attr(result, "condition")))
  }
  invisible(NULL)
}

#' Detect if R1 or R2 of read pair is primary read direction
#' @import GenomicAlignments GenomicRanges
#' @param R1 path to R1 fasta/fastq file
#' @param R2 path to R2 fasta/fastq file
#' @param genomeDir path to STAR index of genome, default:
#' "~/livemount/Bio_data/references/homo_sapiens/STAR_index/genomeDir/"
#' @param nreads numeric, default 1e6 (number of reads to use)
#' @param tx GRangesList, the transcripts to count overlaps on
#' @param td path to tempdir to use, default tempdir()
#' @param threads numeric, default 20
#' @param star path to STAR, default STAR.install()
#' @param keepGenomeLoaded character, default c("LoadAndRemove", "LoadAndKeep", "NoSharedMemory")[1]
#' @return a named numeric (1 or 2), 1 means R1 is primary. Name gives ratio of overlaps of R1 / R2.
#' So 0.07 means R2 is bigger and 7% overlaps in R1 compared to R2.
#' @examples
#' genomeDir <- "~/livemount/Bio_data/references/mus_musculus/STAR_index/genomeDir/"
#' files <- c("~/livemount/Bio_data/raw_data/RNA-seq/PRJNA985729-mus_musculus/SRR24972841_1.fastq",
#'  "~/livemount/Bio_data/raw_data/RNA-seq/PRJNA985729-mus_musculus/SRR24972841_2.fastq")
#' R1 <- files[1]
#' R2 <- files[2]
#' #detect_strand_mode(R1, R2, genomeDir)
detect_strand_mode <- function(R1, R2,
                               genomeDir = "~/livemount/Bio_data/references/homo_sapiens/STAR_index/genomeDir/",
                               nreads = 1e6,
                               tx = loadRegion(list.files(dirname(dirname(genomeDir)), ".db$", full.names = TRUE)[1], "tx"),
                               td = tempdir(), threads = 20, star = STAR.install(),
                               keepGenomeLoaded = c("LoadAndRemove", "LoadAndKeep", "NoSharedMemory")[1]) {
  stopifnot(file.exists(R1) & length(R1) == 1)
  stopifnot(file.exists(R2) & length(R2) == 1)
  stopifnot(file.exists(genomeDir) & length(genomeDir) == 1)
  stopifnot(is.character(td) & length(td) == 1)
  stopifnot(is(tx, "GRangesList"))
  star_load_options <- c("LoadAndRemove", "LoadAndKeep", "NoSharedMemory")
  stopifnot(is.character(keepGenomeLoaded) & keepGenomeLoaded %in% star_load_options)
  if (!dir.exists(td)) dir.create(td)
  bam <- file.path(td, "Aligned.out.bam")
  if (file.exists(bam)) file.remove(bam)
  is_compressed <- grepl(".gz$", R1)

  # Compression check
  cmd_read <- ifelse(is_compressed, "zcat", "-")

  bash_script <- sprintf(
    '%s \\
  --genomeDir %s \\
  --readFilesIn %s %s \\
  --readFilesCommand %s \\
  --readMapNumber %d \\
  --runThreadN %d \\
  --outSAMtype BAM Unsorted \\
  --genomeLoad %s \\
  --outFileNamePrefix %s',
    shQuote(star),
    shQuote(normalizePath(genomeDir)),
    shQuote(normalizePath(R1)),
    shQuote(normalizePath(R2)),
    cmd_read,
    as.integer(nreads),
    threads,
    keepGenomeLoaded,
    shQuote(paste0(normalizePath(td), "/"))
  )


  star_status <- system(bash_script)
  if (star_status != 0) stop("STAR failed to finish without errors!")

  if (!file.exists(bam)) {
    stop("STAR did not produce a BAM file at: ", bam, "\n\nSTAR output:\n", paste(star_status, collapse = "\n"))
  }

  if (file.info(bam)$size == 0) {
    stop("STAR produced empty BAM at: ", bam, "\n\nSTAR output:\n", paste(star_status, collapse = "\n"))
  }

  aln <- readGAlignmentPairs(bam)
  counts <- c(sum(countOverlaps(tx, aln@first, ignore.strand = FALSE)),
              sum(countOverlaps(tx, aln@last, ignore.strand = FALSE)))

  if (min(counts) == 0) stop("No overlaps to transcripts in either direction, either it is the wrong
                             genome/annotation or try to increase nreads")

  max <- which.max(counts)
  names(max) <- round(counts[1] / counts[2], 2)
  return(max)
}
