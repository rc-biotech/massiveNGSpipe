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
                                 adapters_file = fs::path(tempdir(), "adapter_candidates.txt"),
                                 adapter_manual_dir = dirname(file[1])) {
  manual_defined_adapters <- file.path(adapter_manual_dir, "adapters_manual.csv")
  if (file.exists(manual_defined_adapters)) {
    message("Using manually defined adapters from file")
    adapters <- fread(manual_defined_adapters, header = TRUE)
    adapters[, LibraryLayout := "SINGLE"]
    run_to_path <- unlist(run_files_organizer(adapters, dirname(file)))
    index <- run_to_path == ORFik:::pasteDir(path.expand(file))
    if (sum(index) != 1) stop("Could not find Run id matching file: ", file)
    return(adapters[index]$adapter)
  }
  message("- Auto detecting adapter with fastQC candidate list:")
  # Create and save the candidate adapter table.
  candidates <- adapter_list(adapters_file)

  adapters <- fastqc_parse_adapters(file, nreads, adapters_file)

  status <- attributes(adapters)$status
  no_adapter_found <- status == "pass" | nrow(adapters) == 0
  if (no_adapter_found) {
    message("- No adapter found, searching with 5x reads sampled")
    adapters <- fastqc_parse_adapters(file, nreads*5, adapters_file)

    status <- attributes(adapters)$status
    no_adapter_found <- status == "pass" | nrow(adapters) == 0
    if (no_adapter_found) {
      message("- No adapter found, returning disable")
      return("disable")
    }
  }
  adapter_name <- names(which.max(colSums(adapters[, -c("#Position")])))
  found_adapter <- candidates[name %in% adapter_name]$value
  names(found_adapter) <- candidates[name %in% adapter_name]$name
  return(found_adapter)
}

fastqc_parse_adapters <- function(file, nreads, adapters_file =
                                  fs::path(tempdir(), "adapter_candidates.txt")) {

  qc_report_path <- run_fastqc(file, nreads, adapters_file)
  return(read_fastq_adapter_output(qc_report_path))
}

#' Read fastqc output (subset adapter)
#'
#' Parse the report and return the adapter sequence with most hits
#' (or "disable" if none found). Threshold defined by fastqc and
#' "fail" means adapter was found.
read_fastq_adapter_output <- function(qc_report_path) {
  report <- fread(text = readLines(qc_report_path), fill = Inf)

  adapter_section_index <- grep("^>>Adapter Content", report$V1)
  adapters <- report[-c((seq(adapter_section_index + 1)), nrow(report)),]
  colnames(adapters) <- unlist(report[adapter_section_index + 1,])
  score_columns <- adapters[, -1]
  adapters <- cbind(`#Position` = adapters[, `#Position`], score_columns[, (names(score_columns)) := lapply(.SD, as.numeric)])
  status <- report[adapter_section_index, 2][[1]]
  attributes(adapters) <- c(attributes(adapters), status = status, qc_report_path = qc_report_path)
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
      list(name = "RiboLace adapter",    value = "TCTCCTTGCATA"), #TCTCCTTGCATAATCACCAACC
      list(name = "polyA",   value = "AAAAAAAAAA"),
      list(name = "polyT",   value = "TTTTTTTTTT"),
      list(name = "polyN",   value = "NNNNNNNNNN"),
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

#' Get file path from study id
#' @param runs a data.table of study, with columns Run and LibraryLayout
#' @param source_dir character path, directory to match study ids
#' @param exclude formats to exclude, default: c("","json", "html", "md5")
#' @param format formats to search for, default: c(".fastq", ".fq", ".fa", ".fasta")
#' @param compressions character, suffixes to allow for formats,
#' default: c("", ".gz")
#' @param prefixes a list of length 2: prefixes allowed for files, default
#' list(c("", "trimmed_", "collapsed_trimmed_"), c("", "trimmed2_", "collapsed_trimmed2_"))
#' @param paired_end_suffixes a list of length 2: suffixes allowed for files, default
#' list(c("_1", "_2"), c("_R1_001", "_R2_001"))
#' @return a list of character vectors (paths)
#' @export
run_files_organizer <- function(runs, source_dir, exclude = c("","json", "html", "md5"),
                                format = c(".fastq", ".fq", ".fa", ".fasta"),
                                compressions = c("", ".gz"),
                                prefixes = list(c("", "trimmed_", "collapsed_trimmed_"),
                                                c("", "trimmed2_", "collapsed_trimmed2_")),
                                paired_end_suffixes = list(c("_1", "_2"), c("_R1_001", "_R2_001")),
                                extra_files_warning = TRUE) {
  stopifnot(is(runs, "data.table"))
  stopifnot(all(c("Run", "LibraryLayout") %in% colnames(runs)))

  if (!all(dir.exists(source_dir))) stop("Input directory of runs does not exist!")
  files <- list.files(source_dir, full.names = TRUE)
  files <- files[!(tools::file_ext(files) %in% exclude)]
  if (length(files) == 0) stop("There are no files of specified formats + compressions
                               in this directory!")
  if (length(files) > sum(ifelse(runs$LibraryLayout == "SINGLE", 1, 2)) & extra_files_warning)
    warning("You have additional files of wanted format in the folder, please move them!")

  all_files <-
    lapply(seq_len(nrow(runs)), function(i)
      run_files_organizer_internal(i, runs, files, format, compressions,
                                   prefixes, paired_end_suffixes))

  file_vec <- unlist(all_files, use.names = FALSE)
  no_duplicates <- length(unique(file_vec)) == length(file_vec)
  stopifnot(no_duplicates)
  stopifnot(!anyNA(file_vec))
  return(all_files)
}

run_files_organizer_internal <- function(i, runs, files,
                                         format = c(".fastq", ".fq", ".fa", ".fasta"),
                                         compressions = c("", ".gz"),
                                         prefixes = list(c("", "trimmed_", "collapsed_trimmed_"),
                                                         c("", "trimmed2_", "collapsed_trimmed2_")),
                                         paired_end_suffixes = list(c("_1", "_2"), c("_R1_001", "_R2_001"))) {
  for (paired_end_suffix in paired_end_suffixes) {
    filenames <-
      if (runs[i]$LibraryLayout == "PAIRED") {
        paste0(runs[i]$Run, paired_end_suffix)
      } else {
        runs[i]$Run
      }

    file <- grep(filenames[1], files, value = TRUE)
    file2 <- NULL
    if (length(file) == 2 & runs[i]$LibraryLayout == "PAIRED") {
      file2 <- file[2]
      file <- file[1]
    } else if(!is.na(filenames[2])) {
      file2 <- grep(filenames[2], files, value = TRUE)
    }
    if (length(file) != 0) break
  }

  file_1_single_match <- length(file) == 1
  if (!file_1_single_match) {
    if (length(file) == 0) {
      stop("Input File does not exist (Compressions checked ", paste0(compressions, collapse = ", "), ": ",
           runs[i]$Run)
    } else {
      read_1_search <- paste0(filenames[1], format,
                              rep(compressions, each = length(format)))
      read_1_search <- paste0(rep(prefixes[[1]], each = length(read_1_search)), read_1_search)

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

  file_2_multi_match <- length(file2) > 1
  if (file_2_multi_match) {
    if (length(file2) == 0) {
      stop("Input File does not exist (both .gz and unzipped): ",
           runs[i]$Run)
    } else {
      read_2_search <- paste0(filenames[1], format,
                              rep(compressions, each = length(format)))
      read_2_search <- paste0(rep(prefixes[[2]], each = length(read_2_search)), read_2_search)
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
#' pipelines <- pipeline_init_all(config, only_complete_genomes = TRUE, gene_symbols = FALSE, progress_report = FALSE)
#' pipeline <- pipelines[["PRJNA634994"]]
#' # barcode_detector_pipeline(pipeline, FALSE)
#' # debug(detect_adapter_and_trim)
#' # pipeline$study <- pipeline$study[Run == "SRR10846527"]
#' # barcode_detector_pipeline(pipeline, TRUE)
barcode_detector_pipeline <- function(pipeline, redownload_raw_if_needed = TRUE) {
  message("Barcode detection for study (", pipeline$accession, ")")
  dt_stats <- data.table()
  study_all <- pipeline$study

  organism <- names(pipeline$organisms)[1]
  for (organism in names(pipeline$organisms)) {
    message("- ", organism)
    study <- study_all[ScientificName == organism]
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

#' Internal barcode detector
#' @noRd
#' @examples
#' library(massiveNGSpipe)
#' config <- ORFik::config()
#' fastq_dir <- "directory_with_fastq_run" # Change this to correct
#' stopifnot(dir.exists(fastq_dir))
#' # Specify run and paired end status. This is filename in fastq_dir, can have suffix .fastq or .fastq.gz
#' run <- data.table(Run = "SRR12031232", LibraryLayout = "SINGLE")
#' stopifnot(length(list.files(fastq_dir, pattern = paste0(run$Run, ".fastq"))))
#' output_dir <- file.path(config["bam"], basename(fastq_dir)) # Main output dir, change if wanted
#' trim_dir <- file.path(output_dir, "trim") # fastp / barcode output dir
#' dir.create(trim_dir, recursive = TRUE)
#' # Run in debug mode, press n + enter to step through, s + enter to step into sub function
#' debug(massiveNGSpipe:::barcode_detector_single)
#' #massiveNGSpipe:::barcode_detector_single(run,
#' #                                         fastq_dir, output_dir, trim_dir,
#' #                                         redownload_raw_if_needed = FALSE)
barcode_detector_single <- function(study_sample, fastq_dir, process_dir, trimmed_dir,
                                    redownload_raw_if_needed = TRUE, check_at_mean_size = 33,
                                    minimum_size = 26, max_barcode_left_size = 18) {
  sample <- study_sample$Run
  message("-- ", sample)
  stopifnot(is(study_sample, "data.table") && nrow(study_sample) == 1)
  # Pre-define variables
  barcode_sizes <- barcode5p_size <- barcode3p_size <- cut_left_rel_pos <- cut_right_rel_pos <- 0
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
  file <- try(massiveNGSpipe:::run_files_organizer(study_sample, fastq_dir)[[1]][1], silent = TRUE)
  raw_file_exists <- !is(file, "try-error")
  # Download
  if (!raw_file_exists & redownload_raw_if_needed) {
    massiveNGSpipe:::download_sra(study_sample, fastq_dir, compress = FALSE)
    file <- try(run_files_organizer(study_sample, fastq_dir)[[1]][1], silent = TRUE)
    raw_file_exists <- !is(file, "try-error")
  }
  # Trim
  if (!trimmed_file_exists) {
    detect_adapter_and_trim(file, process_dir)
  }
  file_trim <- run_files_organizer(study_sample, trimmed_dir)[[1]][1]

  json_file <- sub("/trimmed_", "/report_", sub("\\.fastq$", ".json", file_trim))
  if (!file.exists(json_file)) stop("Json file from fastp does not exist: ", json_file)
  a <- jsonlite::fromJSON(json_file)
  adapter <- a$adapter_cutting$read1_adapter_sequence
  if (is.null(adapter)) adapter <- "passed"
  adapter_trimmed_reads <- a$adapter_cutting$adapter_trimmed_reads
  if (is.null(adapter_trimmed_reads)) adapter_trimmed_reads <- NA
  reads_no_adapter_removed <- round(100 - (100* (adapter_trimmed_reads / a$summary$before_filtering$total_reads)), 1)
  max_size_before <- a$read1_before_filtering$total_cycles
  mean_size_before <- a$summary$before_filtering$read1_mean_length
  max_size_after <- a$summary$after_filtering$read1_mean_length
  max_size_after_all <- max_size_after

  reads_no_adapter_removed_ORFik <- NA
  fastq_cut <- NULL
  if (max_size_after >= check_at_mean_size) {
    manual_specified_barcodes_file <- file.path(trimmed_dir, "barcodes_manual.csv")
    manual_specified_barcodes_exists <- file.exists(manual_specified_barcodes_file)
    if (manual_specified_barcodes_exists) {
      dt_barcode <- fread(manual_specified_barcodes_file)
      stopifnot(ncol(dt_barcode) == 3 &&
                  all(colnames(dt_barcode) == c("Run", "barcode5p_size", "barcode3p_size")))
      stopifnot(sample %in% dt_barcode$Run)
      barcode_sizes <- unlist(dt_barcode[Run == sample, 2:3])
    } else { # Detect them
      curves <- a$read1_after_filtering[c("quality_curves", "content_curves")]
      curves$content_curves$max <- rowMaxs(as.matrix(setDT(curves$content_curves)), useNames = FALSE)

      barcode_sizes <- barcode_change_point(curves, max_barcode_left_size, max_size_before,
                                            max_size_after, minimum_size, z_score_normalize = FALSE)
      if (max(barcode_sizes) == 0) {
        barcode_sizes <- barcode_change_point(curves, max_barcode_left_size, max_size_before,
                                              max_size_after, minimum_size, z_score_normalize = TRUE)
      }
    }

    barcode5p_size <- barcode_sizes["barcode5p_size"]
    barcode3p_size <- barcode_sizes["barcode3p_size"]


    if (raw_file_exists) { # Internal fastq validations
      # Adapter checker
      if (adapter != "passed") {
        fastq_raw <- readDNAStringSet(file, format = "fastq", nrec = 1e6, use.names = FALSE)
        fastq_raw_trimmed <- remove_adapter_ORFik(fastq_raw, adapter, max.mismatch = 2)
        reads_no_adapter_removed_ORFik <- attr(fastq_raw_trimmed, "statistics")$reads_no_adapter_removed
      }
      # Barcode checker
      fastq <- readDNAStringSet(file_trim, format = "fastq", nrec = 1e5, use.names = FALSE)
      fastq_cut <- trim_flanks_ORFik(fastq, barcode5p_size, barcode3p_size)
    }
  }
  consensus_string_5p <- attr(fastq_cut, "consensus_string_5p")
  consensus_string_3p <- attr(fastq_cut, "consensus_string_3p")
  # subseq(dna, cut_left_rel_pos, width(dna) - cut_right_rel_pos)
  dt_stats_this <- data.table(id = sample, adapter = adapter,
                              barcode_detected = max(barcode_sizes) > 0,
                              max_length_raw = max_size_before,
                              mean_length_raw = mean_size_before,
                              mean_length_adapter_filtered = max_size_after,
                              mean_length_adapter_barcode_filtered = attr(fastq_cut, "mean_size_after_trim"),
                              barcode5p_size, barcode3p_size,
                              `reads_no_adapter_removed fastp(%)` = reads_no_adapter_removed,
                              `reads_no_adapter_removed ORFik(%)` = reads_no_adapter_removed_ORFik,
                              consensus_string_5p, consensus_string_3p)
  if (nrow(dt_stats_this) == 0) stop("Malformed values for barcode table, some are NULL")
  if (is.na(dt_stats_this$barcode_detected)) dt_stats_this[, barcode_detected := FALSE]
  print(dt_stats_this)

  return(dt_stats_this)
}

barcode_change_point <- function(curves, max_barcode_left_size,
                                 max_size_before, max_size_after,
                                 minimum_size, z_score_normalize = FALSE,
                                 Q = 3) {
  cpt_func <- if (z_score_normalize) {
    function(curve_vec, Q) {
      vec <- as.vector(scale(curve_vec))
      if (all(!is.finite(vec))) {
        vec <- rep(0, length(vec))
      } else vec[!is.finite(vec)] <- mean(vec, na.rm = TRUE)
      vec <- vec * abs(vec)
      data.table(changepoint::cpt.mean(vec, Q = Q)@cpts)
    }
  } else {
    function(curve_vec, Q) {
      data.table(changepoint::cpt.mean(curve_vec, Q = Q)@cpts)
    }
  }

  list <- lapply(curves, function(curve_type) {
    lapply(curve_type, function(curve_vec) cpt_func(curve_vec, Q))
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
  return(c(barcode5p_size = barcode5p_size, barcode3p_size = barcode3p_size))
}

run_barcode_detection_and_trim <- function(study_sample, source_dir, target_dir, trimmed_dir,
                                           mode, adapter) {
  # Step 1: Detect barcode
  barcode_dt <- barcode_detector_single(study_sample, source_dir, target_dir, trimmed_dir,
                                        redownload_raw_if_needed = mode == "online")

  if (barcode_dt$barcode_detected) {
    message("Barcode detected")
    # Step 2: Define directory to store original files before barcode removal
    barcode_dir <- file.path(trimmed_dir, "before_barcode_removal")
    # Step 3 & 4: Get and move relevant files
    file <- move_trimmed_files(study_sample, trimmed_dir, barcode_dir)

    # Step 4: Run STAR alignment on the trimmed file with barcode info
    ORFik::STAR.align.single(
      file,
      output.dir = target_dir,
      adapter.sequence = adapter,
      index.dir = "none",
      steps = "tr",
      trim.front = barcode_dt$barcode5p_size,
      trim.tail = barcode_dt$barcode3p_size
    )

    # Step 6: Delete intermediate file
    fs::file_delete(file)
  } else message("No barcodes detected")
  return(barcode_dt)
}

move_trimmed_files <- function(study_sample, trimmed_dir, barcode_dir) {
  # Ignores file 2 in pair for now!
  run <- study_sample$Run

  file_trim_all <- massiveNGSpipe:::run_files_organizer(study_sample, trimmed_dir)[[1]]
  file_trim <- file_trim_all[1]
  json_file <- grep(run, dir(trimmed_dir, "\\.json$", full.names = TRUE), value = TRUE)
  html_file <- grep(run, dir(trimmed_dir, "\\.html$", full.names = TRUE), value = TRUE)

  all_3_old_files <- c(file_trim, json_file, html_file)
  all_3_new_files <- sub("/trimmed_", "/", file.path(barcode_dir, basename(all_3_old_files))) # /trimmed2_
  stopifnot(length(file_trim) > 0)
  stopifnot(!any(duplicated(all_3_new_files)))
  # Step 4: Move old files
  dir.create(barcode_dir, showWarnings = FALSE, recursive = TRUE)
  fs::file_move(all_3_old_files, all_3_new_files)
  return(all_3_new_files[1])
}


detect_adapter_and_trim <- function(file, process_dir, file2 = NULL,
      adapter = try(fastqc_adapters_info(file, adapter_manual_dir = file.path(process_dir, "trim")))) {

  tempfile <- file
  is_fastq <- !grepl("\\.fasta$|\\.fasta\\.gz$", file)
  if (is_fastq) {

    polyN_adapter <- adapter == "NNNNNNNNNN"
    if (is(adapter, "try-error")) {
      message("This is a fasta file, fastqc adapter detection disabled")
      adapter <- "disable"
    } else if (polyN_adapter) {
      message("PolyN adapter found, removing polyN and setting adapter to 'disable'")
      seqs <- readDNAStringSet(file, "fastq", use.names = FALSE, with.qualities = TRUE)
      seqs_sub <- trimLRPatterns(Lpattern="", Rpattern=paste(rep("N", max(width(seqs))), collapse = ""), subject=seqs)
      seqs_sub <- seqs_sub[width(seqs_sub) > 20]
      seqs_sub@elementMetadata$qualities <- subseq(seqs_sub@elementMetadata$qualities, 1, width(seqs_sub))
      names(seqs_sub) <- as.character(seq_along(seqs_sub))
      tempfile <- file.path(tempdir(), basename(file))
      writeXStringSet(seqs_sub, tempfile, format = "fastq")
      adapter <- "AGATCGGAAGAG"
    }
  } else adapter <- "disable" # IF fasta

  ORFik::STAR.align.single(
    tempfile, file2,
    output.dir = process_dir,
    adapter.sequence = adapter,
    index.dir = "none", steps = "tr"
  )
  if (polyN_adapter) {file.remove(tempfile)}
  return(adapter)
}

#' @export
barcodes_manual_assign <- function(pipeline, barcode5p_size, barcode3p_size) {
  stopifnot(length(pipeline) > 0 & is.list(pipeline))
  study <- pipeline$study
  for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    index <- pipeline$organisms[[organism]]$index
    process_dir <- target_dir <- conf["bam"]
    trimmed_dir <- fs::path(process_dir, "trim")

    runs <- study[ScientificName == organism]
    barcodes_manual_assign_table(trimmed_dir, runs$Run, barcode5p_size,
                                 barcode3p_size)
  }
}

barcodes_manual_assign_table <- function(trimmed_dir, run_ids, barcode5p_size,
                                         barcode3p_size) {
  path <- file.path(trimmed_dir, "barcodes_manual.csv")
  dt_barcode <- data.table(Run = run_ids, barcode5p_size, barcode3p_size)
  fwrite(dt_barcode, path)
  message("Barcodes saved to: ", path)
  return(dt_barcode[])
}

#' @export
adapters_manual_assign <- function(pipeline, adapters) {
  stopifnot(length(pipeline) > 0 & is.list(pipeline))
  study <- pipeline$study
  for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    index <- pipeline$organisms[[organism]]$index
    process_dir <- target_dir <- conf["bam"]
    trimmed_dir <- fs::path(process_dir, "trim")

    runs <- study[ScientificName == organism]
    adapters_manual_assign_table(trimmed_dir, runs$Run, adapters)
  }
}

adapters_manual_assign_table <- function(trimmed_dir, run_ids, adapters) {
  path <- file.path(trimmed_dir, "adapters_manual.csv")
  dt_barcode <- data.table(Run = run_ids, adapter = adapters)
  fwrite(dt_barcode, path)
  message("Adapters saved to: ", path)
  return(dt_barcode[])
}

remove_adapter_ORFik <- function(fastq_raw, adapter, max.mismatch = 2,
                                 fixed = TRUE, add_statistics = TRUE) {

  hits <- find_adapter_start_pos(fastq_raw, adapter, max.mismatch, fixed)

  fastq_raw_trimmed <- subseqSafeNonEmptyTrim3p(fastq_raw, hits)

  if (add_statistics) {
    original_widths <- width(fastq_raw)
    trimmed_widths <- width(fastq_raw_trimmed)
    trimmed_too_short <- trimmed_widths < 20
    untrimmed_reads_raw_ORFik <- sum(trimmed_widths == original_widths)
    reads_no_adapter_removed_ORFik <- round((100* ( untrimmed_reads_raw_ORFik / length(fastq_raw_trimmed))), 1)
    attr(fastq_raw_trimmed, "statistics") <-
    data.table(original_reads = length(fastq_raw), trimmed_reads = length(fastq_raw_trimmed),
               reads_no_adapter_removed = reads_no_adapter_removed_ORFik,
               reads_too_short = sum(trimmed_too_short),
               zero_length_reads = sum(width(fastq_raw_trimmed) == 0),
               max_size_reads = sum(trimmed_widths == max(trimmed_widths)),
               mean_readlength = mean(trimmed_widths),
               median_readlength = median(trimmed_widths))
  }
  return(fastq_raw_trimmed)
}

trim_flanks_ORFik <- function(fastq, left = 0, right = 0) {
  stopifnot(!is.character(fastq) | is(fastq, "DNAString") | is(fastq, "DNAStringSet"))
  stopifnot(is.numeric(left) & is.numeric(right))
  stopifnot(min(c(left, right)) >= 0)

  cut_left_rel_pos <- left + 1
  cut_right_rel_pos <- right + 1


  fastq_width <- width(fastq)
  fastq_cut <- subseqSafe(fastq, cut_left_rel_pos, fastq_width - right)

  mean_size_after_trim <- round(mean(width(fastq_cut)), 0)
  attr(fastq_cut, "mean_size_after_trim") <- mean_size_after_trim
  string_5p <- string_3p <- ""
  if (left > 0) {
    string_5p <- consensusString(subseqSafe(fastq, 1, left))
  }
  if (right > 0) {
    right_start <- fastq_width - right + 1
    string_3p <- consensusString(subseqSafe(fastq, right_start, fastq_width))
  }
  attr(fastq_cut, "consensus_string_5p") <- ifelse(length(string_5p) == 0, "", string_5p)
  attr(fastq_cut, "consensus_string_3p") <- ifelse(length(string_3p) == 0, "", string_3p)

  return(fastq_cut)
}

subseqSafeNonEmptyTrim3p <- function(fastq_raw, hits) {
  empty <- lengths(hits) == 0
  h <- hits[!empty] - 1
  h <- unlist(heads(h, 1))
  whole_read_trimmed <- h == 0
  fastq_raw_trimmed <- fastq_raw
  fastq_raw_trimmed[!empty] <- subseqSafe(fastq_raw_trimmed[!empty], 1, h)
  attr(fastq_raw_trimmed, "empty") <- empty
  return(fastq_raw_trimmed)
}

subseqSafe <- function(x, start = NA, end = NA, include_invalid = TRUE) {
  stopifnot(length(start) == 1 | length(start) == length(x))
  stopifnot(length(end) == 1 | length(end) == length(x))
  if (identical(start, NA) & identical(end, NA)) return(subseq(x, start, end))
  width_x <- width(x)
  if (identical(start, NA)) start <- rep(1, length(x))
  if (identical(end, NA)) end <- width_x
  valid <- width_x >= start & width_x >= end & end != 0 & end >= start
  if (length(start) != 1) start <- start[valid]
  if (length(end) != 1) end <- end[valid]

  if (include_invalid) {
    if (any(!valid)) x[!valid] <- ""
    x[valid] <- subseq(x[valid], start, end)
  } else x <- subseq(x[valid], start, end)

  return(x)
}

find_adapter_start_pos <- function(reads, adapter, max.mismatch = 2, fixed = TRUE,
                                   minimum_adapter_start = NA, maximum_adapter_start = NA,
                                   bpparam = MulticoreParam(8, exportglobals = FALSE, log = FALSE)) {
  # if (!all(is.na(c(minimum_adapter_start, maximum_adapter_start)))) {
  #   reads <- subseqSafe(reads, minimum_adapter_start, maximum_adapter_start)
  # }
  adapter_time <- Sys.time()
  message("-- Detecting start positions of 3' adapter")
  if (is.na(minimum_adapter_start)) minimum_adapter_start <- 1

  hits <- bpvec(
    reads,
    FUN = function(chunk_reads, adapter_dna, max.mismatch, fixed) {
      start(IRangesList(vmatchPattern(adapter_dna,
                    chunk_reads,
                    max.mismatch = max.mismatch,
                    fixed = fixed)))
    },
    adapter_dna = DNAString(adapter),
    max.mismatch = max.mismatch,
    fixed = fixed,
    BPPARAM = bpparam
  )
  # vmatchPattern(DNAString(adapter), reads, max.mismatch = max.mismatch, fixed = fixed)
  message("---- Adapters detected per read before minimum pos filter:")
  print(table(lengths(hits)))
  hits <- hits[hits >= minimum_adapter_start]
  message("---- Adapters detected per read after minimum pos filter:")
  print(table(lengths(hits)))
  cat("--- Done Detecting adapter start positions "); print(round(Sys.time() - adapter_time, 1))
  return(hits)
}
