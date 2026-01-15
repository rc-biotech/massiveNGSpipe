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

  dt_summary <- data.table(total_input_seqs = sum(dt$scores),
                           total_input_seqs_collapsed = nrow(dt),
                           total_aligned_reads = sum(dt[N > 0]$scores),
                           total_aligned_reads_collapsed = nrow(dt[N > 0]),
                           total_alignments = sum(dt$scores*dt$N),
                           total_alignments_collapsed = sum(dt$N),
                           unique_aligned_reads = sum(dt[N == 1]$scores),
                           unique_aligned_reads = nrow(dt[N == 1]),
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

  res <- cbind(dt_summary, dt_summary_per_million, dt_summary_relative, dt_summary_relative_percentage)
  return(res)
}

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

#' Detect if R1 or R2 of read pair is primary read direction
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
  stopifnot(dir.exists(td) & length(td) == 1)
  stopifnot(is(tx, "GRangesList"))
  star_load_options <- c("LoadAndRemove", "LoadAndKeep", "NoSharedMemory")
  stopifnot(is.character(keepGenomeLoaded) & keepGenomeLoaded %in% star_load_options)
  if (!dir.exists(td)) dir.create(td)
  nlines <- nreads*4




  bash_script <- sprintf(
    '%s \\
  --genomeDir %s \\
  --readFilesIn %s %s \\
  --readFilesCommand head -n %d \\
  --runThreadN %d \\
  --outSAMtype BAM Unsorted \\
  --genomeLoad %s \\
  --outFileNamePrefix %s/',
    shQuote(star),
    shQuote(normalizePath(genomeDir)),
    shQuote(normalizePath(R1)),
    shQuote(normalizePath(R2)),
    nlines,
    threads,
    keepGenomeLoaded,
    shQuote(normalizePath(td))
  )
  star_status <- system(bash_script)
  if (star_status != 0) stop("STAR failed to finish without errors!")
  bam <- file.path(td, "Aligned.out.bam")
  if (!file.exists(bam) || file.info(bam)$size == 0) {
    stop("STAR did not produce BAM at: ", bam, "\n\nSTAR output:\n", paste(star_out, collapse = "\n"))
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
