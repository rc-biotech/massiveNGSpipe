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
