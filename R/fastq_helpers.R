#' Helper function to find the main adapter in a .fastq.gz file
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
      list(name = "Hakon 2", value = "CACTCGGGCACCAAGGA"),
      list(name = "Hakon 3", value = "GTGTCAGTCACTTCCAGCGG"),
      list(name = "Hakon 4", value = "TGTAGGCACCATC"),
      list(name = "Hakon 5", value = "AAAAAAAAAA"),
      list(name = "Hakon 6", value = "TCGTATGCCGTCTTCTGCTTG")
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
  adapter_name <- names(which.max(colSums(adapters[, -c("#Position")])))
  return(candidates[name == adapter_name]$value)
}
