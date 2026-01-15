
#' Convert Wiktoria template excel sheet to csv controlfile
#' @param path path to xlsx file
#' @param adapter NULL, else character, if specified add adapters to avoid needing to check it.
#' @return a data.table, also saved to disc at file.path(dirname(path), "controlfile.csv")
#' @examples
#' path <- "~/livemount/Bio_data/raw_data/vaccine_mod_bnt_1225-homo_sapiens_x_sars_cov2/multiplexed/Ribosome-profiling_251010.xlsx"
#' adapter <- c(Illumina_TRUESEQ = "AGATCGGAAGAG")
#' excel_Wiktoria_to_csv(path, adapter)
excel_Wiktoria_to_csv <- function(path = "~/livemount/Bio_data/raw_data/vaccine_mod_bnt_1225-homo_sapiens_x_sars_cov2/multiplexed/Ribosome-profiling_251010.xlsx",
                                  adapter = NULL) {
  dt <- as.data.table(readxl::read_xlsx(path, 2))
  dt_out <- dt[, .(`Sample type`, Cells, Condition, `Library name`,`3' Linker index`)]
  colnames(dt_out) <- c("LIBRARYTYPE", "CELLS","CONDITION", "LANE", "BARCODE")
  stopifnot(nrow(dt_out) > 0)
  stopifnot(ncol(dt_out) == 5)

  barcodes_per_lane <- dt_out[, .N, by = .(LANE, BARCODE)]
  if (any(barcodes_per_lane$N != 1)) {
    print(barcodes_per_lane[N != 1,])
    stop("You have duplicated barcodes within the same lane!")
  }
  dt_out[, sample := gsub(" |\\,", "_", gsub("|\\+|\\-", "", paste(LIBRARYTYPE, CELLS, CONDITION)))]
  dt_out[, sample := paste0(sample, "_rep", seq(.N)), by = sample]
  dt_out[, sample := gsub("__", "_", sample)]
  if (!is.null(adapter)) dt_out[, adapter := adapter]
  fwrite(dt_out, file.path(dirname(path), "controlfile.csv"))
  return(dt_out)
}


#' Demultiplex illumina lane fastq files to sample fastq
#'
#' Using a csv controlfile
#' @param path path to directory of multiplexed files, subfolder of raw_data dir
#' as . Files must be lane names like: /raw_data/study/multiplexed/L1_R1_001.fastq.gz.
#' Demultiplexed files will go to directory dirname(path)
#' @param path_out path to output directory for collapsed files
#' @param controlfile path to csv file of demultiplex info of samples for all lanes.
#' Default: file.path(path, "controlfile.csv"). Must contain at least
#' 3 columns: c("LANE", "BARCODE", "sample").
#' @param files character, relative path to all files. Default:
#'  \code{list.files(path, "R1_001\\.(fastq$|fastq.gz$)")}
#' @param umi_5p_size numeric, default 4. Size of 5' end UMI.
#' @param umi_3p_size numeric, default 5. Size of 3' end UMI.
#' @param minimum_bio_fragment_size numeric, default 19. Minimum size of biological fragment,
#' i.e. when can the barcode start for a read that has a biological signal,
#' it will then only use reads that has first adapter hit starting on position >=
#' minimum_bio_fragment_size + umi_5p_size + umi_3p_size + barcode_size.
#' @param compress logical, default TRUE, compress fastq output files
#' @param chunk_size numeric, default 2e7. Number of reads to process per stream
#' repeat.
#' @param bpparam MulticoreParam(8, exportglobals = FALSE, log = FALSE),
#' parallel processing setup for writing split fastq by barcodes.
#' @return invisible(NULL), files saved to disc
#' @importFrom ShortRead FastqStreamer yield sread writeFastq readFastq
#' @examples
#' # multiplexed is subfolder of raw_data study dir
#' path <- "~/livemount/Bio_data/raw_data/vaccine_mod_bnt_1225-homo_sapiens_x_sars_cov2/multiplexed/"
#' #devtools::load_all(); demultiplex(path)
demultiplex <- function(path, path_out = sub("multiplexed/", "trim/SINGLE", sub("/raw_data/", "/processed_data/", path)),
                        controlfile = file.path(path, "controlfile.csv"),
                        files = list.files(path, "R1_001\\.(fastq$|fastq.gz$)"),
                        umi_5p_size = 4, umi_3p_size = 5, minimum_bio_fragment_size = 19,
                        compress = TRUE, compress_collapsed = TRUE, chunk_size = 2e7,
                        bpparam = MulticoreParam(8, exportglobals = FALSE, log = FALSE)) {

  controlfile <- validate_demultiplex_input(path, path_out, files, controlfile)
  suffix <- paste0(".fasta", ifelse(compress_collapsed, ".gz", ""))
  suffix_fastq <- paste0(".fastq", ifelse(compress, ".gz", ""))
  raw_fastq_dir <- dirname(path)
  statistics <- data.table()
  for (file in files) {
    message("- ", file)
    controlfile[, READS_BY_BARCODE := 0]
    pool_name <- sub("_R1\\.001|_R2\\.001|\\.fastq\\.gz$", "", file)
    ctrl <- controlfile[filename == pool_name,]
    if (nrow(ctrl) == 0) stop("Found no matching file in controlfile to Lane: ", pool_name)
    barcodes <- unique(ctrl$BARCODE)
    barcode_sizes <- nchar(barcodes)
    stopifnot(length(unique(barcode_sizes)) == 1)
    barcodes_size <- barcode_sizes[1]

    #ShortRead::readFastq If we need quality
    file <- file.path(path, file)
    adapter <- ctrl$adapter[1]
    if (!isTruthy(adapter)) adapter <- massiveNGSpipe:::fastqc_adapters_info(file)
    if (adapter == "disable") stop("For demultiplexing adapter must be found!")

    # Create a streamer that reads 1e6 reads per chunk
    file_start_time <- Sys.time()
    streamer <- FastqStreamer(file, n = chunk_size)
    on.exit(try(close(streamer), silent = TRUE), add = TRUE)

    bg <- NULL
    chunk <- 1
    repeat{
      message("Chunk - ", chunk)
      chunk_time <- Sys.time()
      message("-- Loading chunk..")
      fq_chunk <- yield(streamer)
      cat("--- Chunk loaded, time used: "); print(round(Sys.time() - chunk_time, 1))
      if (length(fq_chunk) == 0) break

      append <- chunk > 1
      mode <- ifelse(append, "a", "w")
      # # Read the chunk (1e6 reads)
      reads <- sread(fq_chunk)
      # Find adapter start position
      minimum_adapter_start <- minimum_bio_fragment_size + umi_5p_size + umi_3p_size + barcodes_size
      hits <- find_adapter_start_pos(reads, adapter, max.mismatch = 2,
                                     fixed = TRUE,
                                     minimum_adapter_start = minimum_adapter_start,
                                     bpparam = bpparam)
      # Get barcode per read
      barcode_per_read_list <- get_multiplex_barcode_per_read(reads, hits, barcodes)
      # Save matched raw reads split by barcode
      while(chunk != 1 && bg$is_alive()) {
        message("- Waiting for write of raw fastq from last chunk")
        Sys.sleep(1)
      }
      bg <- save_raw_full_fastq_demultiplexed(fq_chunk, barcode_per_read_list, ctrl,
                                        compress, mode, chunk, raw_fastq_dir, bpparam)
      message("- BARCODE removal -> UMI removal -> TRIM -> COLLAPSE")
      # Get adapter trimmed reads and do: barcode removal -> UMI -> TRIM -> COLLAPSE
      barcode_per_sample <- barcode_per_read_list$barcode_per_read
      reads <- lapply(barcode_per_read_list$valid_barcode_unique,
                      function(barcode) {
                        barcode_per_read_list$adapter_reads_trimmed_adapter[barcode_per_sample == barcode]
                      })
      message("-- BARCODE removal")
      reads <- lapply(reads, function(x) subseq(x, end = -nchar(barcodes[1]) - 1))
      message("-- UMI removal")
      reads <- umi_processing(reads, umi_5p_size, umi_3p_size)
      message("-- collapsing")
      reads_collapsed <- lapply(reads, function(x) ORFik:::collapse.fastq.internal(x))
      # reads <- fastqc_mNGSp(reads_collapsed)
      names(reads_collapsed) <- sapply(barcode_per_read_list$valid_barcode_unique,
                                       function(barcode) ctrl[BARCODE == barcode]$sample)

      writeout.fastas(reads_collapsed, path_out, append = append)

      controlfile[filename == pool_name, READS_BY_BARCODE := READS_BY_BARCODE + table(barcode_per_sample)[BARCODE]]

      chunk <- chunk + 1
      cat("Chunk time: "); print(round(Sys.time() - chunk_time, 1))
      cat("Time of file: "); print(round(Sys.time() - file_start_time, 1))
    }

    message("-- Closing the file streamer..")
    close(streamer)

    # re-collapse collapsed chunks
    message("- All chunks done, recollapsing")

    files_to_recollapse <- file.path(out_dir, paste0(names(reads_collapsed), suffix))
    bplapply(files_to_recollapse, recollapse_fastq, BPPARAM = bpparam)


    message("- File done")
    # Full file mode:
    # fq_chunk <- readFastq(file)
    gc()
  }
  # Post processing

  stopifnot(length(list.files(raw_fastq_dir, "fastq")) >= nrow(controlfile))
  stopifnot(length(list.files(raw_fastq_dir, "fastq")) >= nrow(controlfile))
  final_csv_path <- file.path(path, "demultiplex_stats.csv")
  message("Writing final output csv to location: \n", final_csv_path)
  fwrite(controlfile, final_csv_path)
  # TODO: Move files to respective directories
  pool_across_lanes <- final_sample %in% colnames(controlfile)
  if (pool_across_lanes) {
    message("Merge raw fastq seperated by LANE by pool (final_sample column)")
    merge_lane_samples_per_pool_sample(raw_fastq_dir, controlfile, suffix_fastq)



    message("Merge collapsed fasta seperated by LANE by pool (final_sample column)")
    collapsed_pool_files_list <-
      merge_lane_samples_per_pool_sample(path_out, controlfile, suffix)

    bplapply(names(collapsed_pool_files_list), recollapse_fastq, BPPARAM = bpparam)

  }


  # Wait for fastqwrite of background job to finish
  i <- 0
  while(chunk != 1 && bg$is_alive()) {
    message("- Waiting for write of raw fastq from last chunk")
    Sys.sleep(i + 1)
    i <- i + 1
  }
  return(controlfile)
}

get_multiplex_barcode_per_read <- function(reads, adapter_start_pos, barcodes) {
  message("-- Extracting barcodes per read")
  barcode_extract_time <- Sys.time()
  reads_before_adapter <- subseqSafeNonEmptyTrim3p(reads, adapter_start_pos)
  adapter_found <- !attr(reads_before_adapter, "empty")
  barcodes <- as.character(barcodes)
  barcode_size <- nchar(barcodes[1])
  barcode_pos_start <- width(reads_before_adapter)[adapter_found] - barcode_size + 1

  barcode_per_read <- subseqSafe(reads_before_adapter[adapter_found], barcode_pos_start)
  message("- Top 10 barcode regions:")
  barcode_per_read_chr <- as.character(barcode_per_read)
  print(head(sort(table(as.character(barcode_per_read)), decreasing = TRUE), 10))

  valid_barcode <- barcode_per_read_chr %in% as.character(barcodes)
  valid_barcode_unique <- unique(barcodes[barcodes %in% barcode_per_read_chr])

  cat("--- Extracting barcodes per read"); print(round(Sys.time() - barcode_extract_time, 1))
  return(list(adapter_reads_trimmed_adapter = reads_before_adapter[adapter_found],
              barcode_per_read = barcode_per_read_chr,
              adapter_found = adapter_found,
              valid_barcode = valid_barcode,
              valid_barcode_unique = valid_barcode_unique))
}

validate_demultiplex_input <- function(path, path_out, files, controlfile) {
  stopifnot(is.character(path))
  stopifnot(is.character(path_out))
  stopifnot(length(path) == 1)
  stopifnot(length(path_out) == 1)
  stopifnot(dir.exists(path))
  stopifnot(identical(basename(path), "multiplexed"))
  stopifnot(length(files) > 0)
  stopifnot(is.character(controlfile))
  stopifnot(file.exists(controlfile))
  controlfile <- fread(controlfile, header = TRUE)
  stopifnot(all(c("LANE", "BARCODE", "sample") %in%
                  colnames(controlfile)))
  stopifnot(!any(duplicated(controlfile$sample)))
  stopifnot(length(files) != nrow(controlfile))
  stopifnot(length(unique(nchar(controlfile$BARCODE))) == 1)
  controlfile[, filename := paste0("L", LANE)]
  if (!dir.exists(path_out)) dir.create(path_out, FALSE, TRUE)
  return(controlfile)
}

save_raw_full_fastq_demultiplexed <- function(fq_chunk, barcode_per_read_list, ctrl,
                                              compress, mode, chunk, output_dir,
                                              bpparam) {
  message("-- Saving matched chunks by barcode")
  barcode_per_sample <- barcode_per_read_list$barcode_per_read
  adapter_found <- barcode_per_read_list$adapter_found
  valid_barcodes <- barcode_per_read_list$valid_barcode_unique

  # precompute output filenames
  out_files <- file.path(
    output_dir,
    paste0(ctrl[match(valid_barcodes, ctrl$BARCODE)]$sample,
           ".fastq", ifelse(compress, ".gz", ""))
  )
  # build a split lookup
  barcode_groups <- split(seq_along(barcode_per_sample), barcode_per_sample)
  # keep only barcodes we need
  barcode_groups <- barcode_groups[valid_barcodes]
  # now subset FASTQ *here*, once
  which_adapter_found <- which(adapter_found)
  fastq_subsets <- lapply(barcode_groups, function(idx) fq_chunk[which_adapter_found[idx]])

  #Save unmatched reads
  invalid_path <- file.path(output_dir, "unmatched")
  dir.create(invalid_path, showWarnings = FALSE, recursive = TRUE)
  valid_barcode <- barcode_per_read_list$valid_barcode
  invalid_set <- !adapter_found
  invalid_set[adapter_found][!valid_barcode] <- TRUE
  out_file_invalid <- file.path(invalid_path,
                                paste0("unmatched_", ctrl$filename[1], ".fastq", ifelse(compress, ".gz", "")))
  fastq_subsets <- c(fastq_subsets, unmatched = fq_chunk[invalid_set])
  out_files <- c(out_files, out_file_invalid)

  bg <- callr::r_bg(
    func = save_raw_full_fastq_demultiplexed_loop,
    args = list(
      fastq_subsets = fastq_subsets,
      out_files = out_files,
      chunk = chunk,
      compress = compress,
      mode = mode,
      bpparam = bpparam
    )
  )
  return(bg)
}

save_raw_full_fastq_demultiplexed_loop <- function(fastq_subsets,
                                                   out_files,
                                                   chunk,
                                                   compress,
                                                   mode, bpparam) {
  chunk_write_time <- Sys.time()
  BiocParallel::bpmapply(
    FUN = function(fq_subset, out_file, chunk, compress, mode) {
      if (chunk == 1 && file.exists(out_file))
        file.remove(out_file)

      ShortRead::writeFastq(fq_subset,
                            out_file,
                            mode = mode,
                            compress = compress)
      message(out_file)

      return(out_file)
    },
    fastq_subsets,
    as.list(out_files),
    MoreArgs = list(chunk = chunk, compress = compress, mode = mode),
    BPPARAM = bpparam,
    SIMPLIFY = FALSE
  )
  cat("--- Barcode splits saved, time used: "); print(round(Sys.time() - chunk_write_time, 1))
}


#' Process UMIs
#' @param reads a list of DNAStringSet
#' @param firstN integer, number of nt for 5' UMI part
#' @param lastN integer, number of nt for 3' UMI part
#' @param saveUMIs logical, default FALSE.
#' @return a list of DNAStringSet processed.
umi_processing <- function(reads, firstN, lastN, saveUMIs = FALSE) {
  if (saveUMIs) {
    reads <- lapply(reads, function(x) {
      umi1 <- character()
      umi2 <- character()
      if (firstN != 0) umi1 <- as.character(subseq(x, start = 1, end = firstN))
      if (lastN != 0) umi2 <- as.character(subseq(x, start = -lastN, end = -1))
      umi <- paste0(umi1,umi2)
      names(x) <- umi
      return(x)
    })
  }
  reads <- lapply(reads, function(x) subseq(x, start= firstN + 1, end = -lastN - 1))
  return(reads)
}

#' Write a list of fasta files
#'
#' @param reads a named list of DNAStringSets, names are file prefixes,
#' will append .fasta or .fasta.gz depending on if compress is TRUE.
#' @param out_dir output directory
#' @param append logical, default FALSE. Append to file (FALSE) or Overwrite (TRUE)
#' @param compress logical, default FALSE. Compress to .gz or not.
#' @return invisible(NULL), files saved to disc
writeout.fastas <- function(reads, out_dir, append = FALSE,
                            compress = FALSE) {
  suffix <- paste0(".fasta", ifelse(compress, ".gz", ""))
  mapply(function(x,y) writeXStringSet(x, file.path(out_dir, paste0(y, suffix)),
                                       append = append, compress = compress),
         reads, names(reads))
}


#' Given ..., demultiplex to fastq
#' @param path path to directory of
#' @param csv path to csv file of
#' @return invisible(NULL), files saved to disc
demultiplex_from_bin <- function(path, csv, demultiplex_tool_bin = "~/bin/demultiplexer") {
  # TODO: copy over:
  #
  call <- paste(demultiplex_tool_bin, "-i", path, "-c", csv)
  ret <- system(call)

  if (ret != 0) stop("Demultiplexing failed")
  return(invisible(NULL))
}

merge_lane_samples_per_pool_sample <- function(raw_fastq_dir, controlfile, suffix) {
  raw_split_by_lane_dir <- file.path(raw_fastq_dir, "split_by_lane")
  dir.create(raw_split_by_lane_dir, showWarnings = FALSE, recursive = TRUE)
  lane_files_raw_base <- paste0(controlfile$sample, suffix)
  lane_files_raw <- file.path(
    raw_fastq_dir,
    lane_files_raw_base
  )
  fs::file_move(lane_files_raw, raw_split_by_lane_dir)
  lane_files_raw <- file.path(raw_split_by_lane_dir, lane_files_raw_base)
  pool_files_raw_base <- paste0(controlfile$final_sample, suffix)
  pool_files_raw <- file.path(raw_fastq_dir, pool_files_raw_base)

  pool_list_from_lanes <- split(lane_files_raw, pool_files_raw)
  message("Lane to pool split list:")
  print(pool_list_from_lanes)
  ORFik::mergeFastq(pool_list_from_lanes)
  return(pool_list_from_lanes)
}

recollapse_fastq <- function(file) {
  seqs <- readDNAStringSet(file, use.names = TRUE)
  replicates <- data.table(seqs = as.character(seqs), N = as.numeric(sub(".*_x", "", names(seqs))))

  ORFik:::local_DTthreads(1)
  replicates <- replicates[, .(N = sum(N)), by = seqs][order(N, decreasing = TRUE),
  ]
  new_seqs <- replicates$seqs
  names(new_seqs) <- ORFik:::collapse_header_set(replicates$N, seq.int(nrow(replicates)))
  new_seqs <- new_seqs[nchar(new_seqs) >=20]
  writeXStringSet(DNAStringSet(new_seqs, use.names = TRUE), file)
}
