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

run_files_organizer <- function(runs, source_dir, exclude = c("json", "html")) {
  all_files <-
    lapply(seq_len(nrow(runs)), function(i) {
      # browser()
      filenames <-
        if (runs[i]$LibraryLayout == "PAIRED") {
          paste0(runs[i]$Run, c("_1", "_2"))
        } else {
          paste0(runs[i]$Run)
        }
      format <- c(".fastq", ".fq", ".fa", ".fasta")
      compressions <- c("", ".gz")
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
               file)
        } else {
          read_1_search <- paste0(filenames[1], format,
                                  rep(compressions, each = length(format)))
          read_1_search <- paste0(rep(prefix1, each = length(read_1_search)), read_1_search)

          temp_file <- file[basename(file) %in% read_1_search]
          if (length(temp_file) != 1) {
            stop("File format could not be detected",
                 file[1])
          } else file <- temp_file
        }
      }
      if (!is.null(file2) && length(file2) != 1) {


        if (length(file2) == 0) {
          stop("File does not exist to trim (both .gz and unzipped): ",
               file2)
        } else {
          read_2_search <- paste0(filenames[1], format,
                                  rep(compressions, each = length(format)))
          read_2_search <- paste0(rep(prefix2, each = length(read_2_search)), read_2_search)
          temp_file <- file2[basename(file2) %in% read_2_search]
          if (length(temp_file) != 1) {
            stop("File format could not be detected",
                 file[2])
          } else file2 <- temp_file
        }
      }
      return(c(file, file2))
    })
  no_duplicates <- length(unique(unlist(all_files, use.names = F))) == length(unlist(all_files, use.names = F))
  stopifnot(no_duplicates)
  return(all_files)
}
