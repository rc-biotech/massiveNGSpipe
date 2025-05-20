#' Helper function to find the main adapter in a .fastq.gz file
#' @param file full path to fastq file to check with qc. Supports both
#' ".gz" and non compressed.
#' @param nreads integer, default 6000000 (6e6),
#' number of reads to read in for check
#' @param adapters_file = fs::path(tempdir(), "adapter_candidates.txt"),
#'  will use a internal set of known adapters, saved to a tempfile.
#' @return adapter candidates
#' @export
#' @examples
#' # Define parameters
#' n_strings <- 500000  # Number of strings
#' fixed_total_size <- 51
#' string_length <- seq(20, 33)  # Length of each string
#' string_lengths <- sample(string_length, n_strings, replace = TRUE)  # Length of each string
#' letters <- c("A", "T", "C", "G")  # DNA letters
#' # Generate random DNA strings
#' all_chars <- DNAString(paste(sample(letters, sum(string_lengths), replace = TRUE), collapse = ""))
#' dna_strings <- DNAStringSet(successiveViews(all_chars, width = string_lengths))
#' # View the first few strings
#' head(dna_strings)
#'
#' # Make barcode, dna, adapter
#' illumina_true_seq <- DNAStringSet(DNAString("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"))
#' barcode5p <- DNAStringSet(DNAString("AATTGG"))
#' barcode3p <- DNAStringSet(DNAString("AATTGG"))
#' dna <- DNAStringSet(paste0(barcode5p, dna_strings, barcode3p, illumina_true_seq))
#' dna <- subseq(dna, 1, fixed_total_size)
#' names(dna) <- paste0("seq_", seq(length(dna)))
#' qualities <- substr("#<GGGGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFGGGGGGGGG", 1, fixed_total_size)
#' mcols(dna)$qualities <- rep(BStringSet(qualities), length(dna))
#'
#' tempfile <- tempfile(fileext = ".fastq")
#' writeXStringSet(dna, tempfile, format = "fastq")
#'
#' fastqc_adapters_info(tempfile, 1e6)
fastqc_adapters_info <- function(file, nreads = 2e6,
                                 adapters_file = fs::path(tempdir(), "adapter_candidates.txt")) {
  message("- Auto detecting adapter with fastQC candidate list:")
  # Create and save the candidate adapter table.
  candidates <- adapter_list(adapters_file)

  adapters <- fastqc_parse_adapters(file, nreads, adapters_file)

  status <- attributes(adapters)$status
  no_adapter_found <- status == "pass" | nrow(adapters) == 0
  if (no_adapter_found) {
    return("disable")
  }
  adapter_name <- names(which.max(colSums(adapters[, -c("#Position")])))
  found_adapter <- candidates[name %in% adapter_name]$value
  names(found_adapter) <- candidates[name %in% adapter_name]$name
  return(found_adapter)
}

fastqc_parse_adapters <- function(file, nreads, adapters_file) {

  qc_report_path <- run_fastqc(file, nreads, adapters_file)

  # Parse the report and return the adapter sequence with most hits
  # (or "disable" if none found).
  report <- fread(text = readLines(qc_report_path), fill = Inf)

  adapter_section_index <- grep("^>>Adapter Content", report$V1)

  adapters <- report[-c((seq(adapter_section_index + 1)), nrow(report)),]
  colnames(adapters) <- unlist(report[adapter_section_index + 1,])
  adapters[, (names(adapters)) := lapply(.SD, as.numeric)]
  status <- report[adapter_section_index, 2][[1]]
  attributes(adapters) <- c(attributes(adapters), status = status)
  return(adapters)
}

run_fastqc <- function(file, nreads, adapters_file) {
  stopifnot(file.exists(file))
  stopifnot(file.exists(adapters_file))

  file_extension <- tools::file_ext(file)
  cat <- ifelse(file_extension == "gz", "zcat ", "cat ")
  temp_dir_unique <- tempfile()
  dir.create(temp_dir_unique)
  nreads <- nreads*4 # 4 lines per read

  res <- system(paste0(
    cat, file,
    " 2>/dev/null | head -n ", format(nreads, scientific = FALSE),
    " | fastqc stdin --extract",
    " --adapters ", adapters_file,
    " -o ", temp_dir_unique
  ))
  qc_report_path <- fs::path(temp_dir_unique, "stdin_fastqc", "fastqc_data.txt")
  return(qc_report_path)
}

adapter_list <- function(candidates_file = fs::path(tempdir(), "adapter_candidates.txt")) {

  if (!file.exists(candidates_file)) {
    candidates <- data.table::rbindlist(list(
      list(name = "Illumina Universal Adapter",    value = "AGATCGGAAGAG"), #AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
      list(name = "Illumina Small RNA 3' Adapter", value = "TGGAATTCTCGG"),
      list(name = "Illumina Small RNA 5' Adapter", value = "GATCGTCGGACT"),
      list(name = "Nextera Transposase Sequence",  value = "CTGTCTCTTATA"),
      list(name = "SOLID Small RNA Adapter",       value = "CGCCTTGGCCGT"),
      list(name = "Ingolia 2012 adapter",          value = "CTGTAGGCACCA"),#CTGTAGGCACCATCAAT
      list(name = "Illumina Uni. Adapter var2",    value = "ATCGTAGATCGG"), #ATCGTAGATCGGAAG
      list(name = "polyA",   value = "AAAAAAAAAA"),
      list(name = "polyT",   value = "TTTTTTTTTT"),
      list(name = "New CAC", value = "CACTCGGGCACC"), #CACTCGGGCACCAAGGA
      list(name = "New GTG", value = "GTGTCAGTCACT"),#GTGTCAGTCACTTCCAGCGG
      list(name = "New TGT", value = "TGTAGGCACCAT"), #TGTAGGCACCATC
      list(name = "New TCG", value = "TCGTATGCCGTC") #TCGTATGCCGTCTTCTGCTTG
    ))
    data.table::fwrite(
      candidates,
      file = candidates_file, sep = "\t", col.names = FALSE
    )
  } else {
    candidates <- fread(candidates_file, header = FALSE)
    stopifnot(ncol(candidates_file) == 2)
    colnames(candidates) <- c("name", "value")
  }
  return(candidates)
}

#' Fastq/fastq file organizer
run_files_organizer <- function(runs, source_dir, exclude = c("json", "html"),
                                format = c(".fastq", ".fq", ".fa", ".fasta"),
                                compressions = c("", ".gz")) {
  all_files <-
    lapply(seq_len(nrow(runs)), function(i) {
      # browser()
      filenames <-
        if (runs[i]$LibraryLayout == "PAIRED") {
          paste0(runs[i]$Run, c("_1", "_2"))
        } else {
          paste0(runs[i]$Run)
        }


      prefix1 <- c("", "trimmed_")
      prefix2 <- c("", "trimmed2_")

      format_used <- filenames[1]
      file <- list.files(source_dir, filenames[1], full.names = TRUE)
      file <- file[!(tools::file_ext(file) %in% exclude)]
      file2 <- NULL
      if(!is.na(filenames[2])) {
        file2 <- list.files(source_dir, filenames[2], full.names = TRUE)
        file2 <- file2[!(tools::file_ext(file2) %in% exclude)]
      }
      if (length(file) != 1) {
        if (length(file) == 0) {
          stop("File does not exist to trim (both .gz and unzipped): ",
               runs[i]$Run)
        } else {
          read_1_search <- paste0(filenames[1], format,
                                  rep(compressions, each = length(format)))
          read_1_search <- paste0(rep(prefix1, each = length(read_1_search)), read_1_search)

          temp_file <- file[basename(file) %in% read_1_search]
          if (length(temp_file) != 1) {
            if (length(temp_file) == 0) {
              file <- runs[i]$Run
            }
            stop("File format could not be detected",
                 file[1])
          } else file <- temp_file
        }
      }
      if (!is.null(file2) && length(file2) != 1) {


        if (length(file2) == 0) {
          stop("File does not exist to trim (both .gz and unzipped): ",
               runs[i]$Run)
        } else {
          read_2_search <- paste0(filenames[1], format,
                                  rep(compressions, each = length(format)))
          read_2_search <- paste0(rep(prefix2, each = length(read_2_search)), read_2_search)
          temp_file <- file2[basename(file2) %in% read_2_search]
          if (length(temp_file) != 1) {
            if (length(temp_file) == 0) {
              file <- rep(runs[i]$Run, 2)
            }
            stop("File format could not be detected",
                 file[2])
          } else file2 <- temp_file
        }
      }
      return(c(file, file2))
    })

  file_vec <- unlist(all_files, use.names = FALSE)
  no_duplicates <- length(unique(file_vec)) == length(file_vec)
  stopifnot(no_duplicates)
  return(all_files)
}

#' NGS barcode detector
#' Detects both 5' and 3' barcodes
#'
#' @param pipeline a pipeline object
#' @param redownload_raw_if_needed logical, default TRUE.
#' For additional stats redownload fastq if removed.
#' @return data.table with stats for all runs in experiment
#' @examples
#' config <- pipeline_config()
#' pipelines <- pipeline_init_all(config, only_complete_genomes = TRUE, gene_symbols = FALSE, simple_progress_report = FALSE)
#' pipeline <- pipelines[["PRJNA634994"]]
#' barcode_detector(pipeline, FALSE)
barcode_detector_pipeline <- function(pipeline, redownload_raw_if_needed = TRUE) {
  message("Barcode detection for study (", pipeline$accession, ")")
  dt_stats <- data.table()
  study_all <- pipeline$study

  organism <- names(pipeline$organisms)[1]
  for (organism in names(pipeline$organisms)) {
    message("- ", organism)
    study <- study[ScientificName == organism]
    if (!all(study$LibraryLayout == "SINGLE")) {
      message("Paired end data, can only detect barcode for single end, next experiment..")
      next
    }
    index <- pipeline$organisms[[organism]]$index
    conf <- pipeline$organisms[[organism]]$conf
    fastq_dir <- conf["fastq"]
    process_dir <- conf["bam"]
    trimmed_dir <- fs::path(process_dir, "trim")

    i <- 1
    dt_stats_org <- data.table()

    dt_stats_org <- rbindlist(lapply(seq(nrow(study)), function(i)
      barcode_detector_single(study[i], fastq_dir, process_dir, trimmed_dir,
                              redownload_raw_if_needed)))

    dt_stats <- rbindlist(list(dt_stats, dt_stats_org))
  }
  return(dt_stats)
}

barcode_detector_single <- function(study_sample, fastq_dir, process_dir, trimmed_dir,
                                    redownload_raw_if_needed = TRUE, check_at_mean_size = 34,
                                    minimum_size = 28, max_barcode_left_size = 18) {
  sample <- study_sample$Run
  message("-- ", sample)
  stopifnot(is(study_sample, "data.table") && nrow(study_sample) == 1)
  barcode5p_size <- barcode3p_size <- cut_left_rel_pos <- cut_right_rel_pos <- 0
  consensus_string_5p <- consensus_string_3p <- ""

  json_file <- grep(sample, dir(trimmed_dir, "\\.json$", full.names = TRUE), value = TRUE)
  if (length(json_file) == 1) {
    trim_stats <- ORFik::trimming.table(trimmed_dir, json_file, TRUE)
    if (trim_stats$trim_mean_length <= check_at_mean_size) {
      return(data.table(id = sample, adapter = trim_stats$adapter,
                        barcode_detected = FALSE,
                        max_length_raw = trim_stats$raw_mean_length,
                        mean_length_raw = trim_stats$raw_mean_length,
                        mean_length_adapter_filtered = trim_stats$trim_mean_length,
                        mean_length_adapter_barcode_filtered = NA,
                        barcode5p_size, barcode3p_size,
                        `reads_no_adapter_removed fastp(%)` = NA,
                        `reads_no_adapter_removed ORFik(%)` = NA,
                        consensus_string_5p, consensus_string_3p))
    }
  }



  file_trim <- try(massiveNGSpipe:::run_files_organizer(study_sample, trimmed_dir)[[1]], silent = TRUE)
  trimmed_file_exists <- !is(file_trim, "try-error")
  file <- file.path(fastq_dir, paste0(sample, ".fastq"))
  raw_file_exists <- file.exists(file)
  # Download
  if (!raw_file_exists & redownload_raw_if_needed) {
    massiveNGSpipe:::download_sra(study_sample, fastq_dir, compress = FALSE)
  }
  # Trim
  if (!trimmed_file_exists) {
    adapter <- massiveNGSpipe:::fastqc_adapters_info(file)

    ORFik::STAR.align.single(
      file, NULL,
      output.dir = process_dir,
      adapter.sequence = adapter,
      index.dir = "", steps = "tr"
    )
  }

  file_trim <- massiveNGSpipe:::run_files_organizer(study_sample, trimmed_dir)[[1]]

  a <- jsonlite::fromJSON(sub("/trimmed_", "/report_", sub("\\.fastq$", ".json", file_trim)))
  adapter <- a$adapter_cutting$read1_adapter_sequence
  if (is.null(adapter)) adapter <- "passed"
  adapter_trimmed_reads <- a$adapter_cutting$adapter_trimmed_reads
  if (is.null(adapter_trimmed_reads)) adapter_trimmed_reads <- NA
  reads_no_adapter_removed <- round(100 - (100* (adapter_trimmed_reads / a$summary$before_filtering$total_reads)), 1)
  reads_no_adapter_removed_ORFik <- NA
  max_size_before <- a$read1_before_filtering$total_cycles
  mean_size_before <- a$summary$before_filtering$read1_mean_length
  max_size_after <- a$summary$after_filtering$read1_mean_length
  max_size_after_all <- max_size_after
  fastq_cut <- NULL
  if (max_size_after > check_at_mean_size) {

    curves <- a$read1_after_filtering[c("quality_curves", "content_curves")]
    list <- lapply(curves, function(curve_type) {
      lapply(curve_type, function(curve_vec) {
        data.table(changepoint::cpt.mean(curve_vec, Q = 3)@cpts)
      })
    })

    set <- unlist(list)
    tab <- table(set)
    tab <- tab[-length(tab)]
    tab <- tab[as.numeric(names(tab)) <= max_size_after]
    pos <- as.numeric(names(tab))

    cut_left <- as.numeric(names(which.max(rev(tab[pos < max_barcode_left_size]))))
    pos_high <- as.numeric(names(which.max(tab[pos > max_size_before-cut_left])))

    #max(pos[pos > cut_left*2])

    barcode5p_size <- max(0, cut_left, na.rm = TRUE)
    barcode3p_size <- max(0, max_size_before - pos_high - 2, na.rm = TRUE)
    barcodes_detected_too_big <- (max_size_after - (barcode5p_size + barcode3p_size)) < minimum_size
    if (barcodes_detected_too_big) { # Reduce 3p barcode size
      amount_too_big <- minimum_size - (max_size_after - (barcode5p_size + barcode3p_size))
      barcode3p_size <- max(barcode3p_size - amount_too_big, 0)
      barcodes_detected_too_big <- (max_size_after - (barcode5p_size + barcode3p_size)) < minimum_size
      if (barcodes_detected_too_big) { # Reduce 5p barcode siz
        amount_too_big <- minimum_size - (max_size_after - (barcode5p_size + barcode3p_size))
        barcode5p_size <- max(barcode5p_size - amount_too_big, 0)
      }
    }

    cut_right_rel_pos <- barcode3p_size + 1
    cut_left_rel_pos <- barcode5p_size + 1


    # a$summary$after_filtering$read1_mean_length
    # x <- rbindlist(list[[1]], fill = TRUE)
    if (raw_file_exists & adapter != "passed") {
      fastq_raw <- readDNAStringSet(file, format = "fastq", nrec = 100000)
      adapter_ext <- paste0(adapter, paste0(rep("N", max_size_before - nchar(adapter)), collapse = ""))
      fastq_raw_trimmed <- trimLRPatterns(Lpattern = "", Rpattern = adapter_ext, subject = fastq_raw, max.Rmismatch = 3)
      fastq_raw_trimmed_len <- fastq_raw_trimmed[width(fastq_raw_trimmed) > 20]
      untrimmed_reads_raw_ORFik <- sum(width(fastq_raw_trimmed_len) == max(width(fastq_raw_trimmed_len)))
      reads_no_adapter_removed_ORFik <- round(100 - (100* ( untrimmed_reads_raw_ORFik / length(fastq_raw_trimmed_len))), 1)
    }


    fastq <- readDNAStringSet(file_trim, format = "fastq", nrec = 100000)
    fastq_cut <- subseq(fastq, cut_left_rel_pos, width(fastq) - cut_right_rel_pos)
    max_size_after_all <- mean(width(fastq_cut))
    if (cut_left_rel_pos > 1) {
      consensus_string_5p <- consensusString(subseq(fastq, 1, cut_left_rel_pos - 1))
    }
    if (cut_right_rel_pos > 1) {
      consensus_string_3p <- consensusString(subseq(fastq, width(fastq) - cut_right_rel_pos, width(fastq)))
    }
  }
  # subseq(dna, cut_left_rel_pos, width(dna) - cut_right_rel_pos)
  dt_stats_this <- data.table(id = sample, adapter = adapter,
                              barcode_detected = (cut_left_rel_pos > 1) | (cut_right_rel_pos > 1),
                              max_length_raw = max_size_before,
                              mean_length_raw = mean_size_before,
                              mean_length_adapter_filtered = max_size_after,
                              mean_length_adapter_barcode_filtered = round(max_size_after_all, 0),
                              barcode5p_size, barcode3p_size,
                              `reads_no_adapter_removed fastp(%)` = reads_no_adapter_removed,
                              `reads_no_adapter_removed ORFik(%)` = reads_no_adapter_removed_ORFik,
                              consensus_string_5p, consensus_string_3p)
  if (nrow(dt_stats_this) == 0) stop("Malformed values for barcode table, some are NULL")
  if (is.na(dt_stats_this$barcode_detected)) dt_stats_this[, barcode_detected := FALSE]

  print(dt_stats_this)

  return(dt_stats_this)
}
