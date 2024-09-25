#' Helper function to find the main adapter in a .fastq.gz file
#' @param file full path to fastq file to check with qc. Supports both
#' ".gz" and non compressed.
#' @return adapter candidates
fastqc_adapters_info <- function(file) {
  message("- Auto detecting adapter with fastQC candidate list:")
  # Create and save the candidate adapter table.

  candidates_file <- fs::path(tempdir(), "adapter_candidates.txt")
  if (!file.exists(candidates_file)) {
    candidates <- data.table::rbindlist(list(
      list(name = "Illumina Universal Adapter", value = "AGATCGGAAGAGC"),
      list(name = "Illumina Small RNA 3' Adapter", value = "TGGAATTCTCGG"),
      list(name = "Illumina Small RNA 5' Adapter", value = "GATCGTCGGACT"),
      list(name = "Nextera Transposase Sequence", value = "CTGTCTCTTATA"),
      list(name = "SOLID Small RNA Adapter", value = "CGCCTTGGCCGT"),
      list(name = "Ingolia 2012 adapter", value = "CTGTAGGCACCATCAAT"),
      list(name = "Illumina Uni. Adapter var2", value = "ATCGTAGATCGGAAG"),
      list(name = "polyA", value = "AAAAAAAAAA"),
      list(name = "Hakon 2", value = "CACTCGGGCACCAAGGA"),
      list(name = "Hakon 3", value = "GTGTCAGTCACTTCCAGCGG"),
      list(name = "Hakon 4", value = "TGTAGGCACCATC"),
      list(name = "Hakon 5", value = "TCGTATGCCGTCTTCTGCTTG")
    ))
    data.table::fwrite(
      candidates,
      file = candidates_file, sep = "\t", col.names = FALSE
    )
  } else {
    candidates <- fread(candidates_file, header = FALSE)
    colnames(candidates) <- c("name", "value")
  }


  # Run FastQC on first 1.5mil reads using the above table for adapter
  # detection.
  file_extension <- tools::file_ext(file)
  cat <- ifelse(file_extension == "gz", "zcat ", "cat ")
  temp_dir_unique <- tempfile()
  dir.create(temp_dir_unique)
  system(paste0(
    cat, file,
    " | head -n 6000000",
    " | fastqc stdin --extract",
    " --adapters ", candidates_file,
    " -o ", temp_dir_unique
  ))

  # Parse the report and return the adapter sequence with most hits
  # (or "disable" if none found).
  report <- readLines(fs::path(temp_dir_unique, "stdin_fastqc", "fastqc_data.txt"))
  status <- strsplit(report[grep("^>>Adapter Content", report)], "\t")[[1]][2]
  if (status == "pass") {
    return("disable")
  }
  adapters <- data.table::fread(
    text = head(report[-(1:grep("^>>Adapter Content", report))], -1)
  )
  if (nrow(adapters) == 0) {
    return("disable")
  }
  adapter_name <- names(which.max(colSums(adapters[, -c("#Position")])))
  return(candidates[name == adapter_name]$value)
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
