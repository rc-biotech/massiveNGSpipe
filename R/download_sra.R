download_sra <- function(info, outdir, compress = TRUE,
                         sratoolkit_path =
                           fs::path_dir(ORFik::install.sratoolkit()),
                         delete_srr_preformat = TRUE) {
  fs::dir_create(outdir)
  outdir <- fs::path_real(outdir)

  accessions <- info$Run
  if (is.null(accessions) | (length(accessions) == 0)) {
    stop("Could not find SRR numbers in 'info' i.e. the study table")
  }

  message("Downloading SRA runs:")
  for (accession in accessions) {
    PAIRED_END <- info$LibraryLayout[info$Run == accession] == "PAIRED"
    download_raw_srr(accession, outdir, compress, sratoolkit_path, PAIRED_END)
  }
  delete_existing_preformat_files(outdir, accessions, delete_srr_preformat)
}

download_raw_srr <- function(accession, outdir, compress = TRUE,
                             sratoolkit_path =
                                 fs::path_dir(ORFik::install.sratoolkit()),
                             PAIRED_END) {
  # Check if the run is available on AWS (not always the case, e.g.
  # SRR4000302)
  if (download_sra_aws(accession, outdir) || download_sra_ascp(accession, outdir)) {
    sra_to_fastq(accession, outdir, compress, sratoolkit_path, PAIRED_END)
  } else {
    download_sra_ebi(accession, outdir, compress)
  }
  return(cleanup_and_validate_fastq_download(accession, outdir, PAIRED_END,
                                             compress))
}

download_sra_aws <- function(run_accession, outdir, aws_bin = "aws",
                             aws_mirror = "s3://sra-pub-run-odp/sra") {

  file_aws_url <- file.path(aws_mirror, run_accession)
  check_file_on_aws <- system2(aws_bin, c(
    "--no-sign-request", "s3", "ls",
    file_aws_url
  ), stdout = FALSE)
  run_exists_on_aws <- check_file_on_aws == 0

  if (run_exists_on_aws) {
    download_status <- system2(aws_bin, c(
      "--no-sign-request", "s3", "sync",
      file_aws_url,
      outdir
    ))
    succesful_download_aws <- download_status == 0
    stopifnot("Run found on aws, but awscli exited with non-zero exit code" = succesful_download_aws)
    run_exists_on_aws <- succesful_download_aws
    # TODO: Add file is downloaded check too (this is not enough)
  }
  return(run_exists_on_aws)
}

download_sra_ascp <- function(run_accession, outdir,
                              ascp_bin = fs::path_home(".aspera/connect/bin/ascp"),
                              ascp_private_key = "~/.aspera/connect/etc/asperaweb_id_dsa.openssh",
                              resume_level = 1,
                              ssh_port = 33001,
                              max_transfer_rate = "300m",
                              ebi_server_run_path = find_ascp_srr_url(run_accession)) {
  # EBI servers via Aspera.
  message("Falling back to ascp")
  out_path <- fs::path_join(c(outdir, run_accession))
  ret <- system2(ascp_bin, c(
    "-T",
    "-k", resume_level,
    "-l", max_transfer_rate,
    "-P", ssh_port,
    "-i", ascp_private_key,
    ebi_server_run_path, out_path
  ))
  stopifnot("ascp exited with non-zero exit code" = ret == 0)
}

download_sra_ebi <- function(accession, outdir, compress) {
  message("Using ebi fallback")
  res <- try(ORFik::download.SRA(accession, outdir, rename = FALSE))
  if (is(res, "try-error")) {
    message("Using fastq-dump fallback")
    res <- ORFik::download.SRA(accession, outdir, use.ebi.ftp = FALSE,
                               rename = FALSE, compress = compress)
  } else {
    if (!compress) {
      lapply(res, function(fastq) R.utils::gunzip(fastq, overwrite = TRUE))
    }
  }
}

sra_to_fastq <- function(accession, outdir, compress = TRUE,
                                  sratoolkit_path = fs::path_dir(ORFik::install.sratoolkit()),
                                  PAIRED_END, progress_bar = TRUE,
                                  fasterq_temp_dir = check_tempdir_has_space_sra_to_fastq(accession, outdir, tempdir(), PAIRED_END)) {
  message("Extracting SRA run:", accession)
  sra_path <- fs::path_join(c(outdir, accession))
  if (!file.exists(sra_path)) {
    stop(".SRA file does not exist to convert to fastq")
  }
  binary <- "fasterq-dump"
  sratoolkit_binary <- fs::path_join(c(sratoolkit_path, binary))
  stopifnot(file.exists(sratoolkit_binary))

  progress_bar <- ifelse(progress_bar, "-p", "")
  using_temp_dir <- fasterq_temp_dir == tempdir()
  threads <- ifelse(using_temp_dir, "--threads 24", "--threads 6")
  temp_path <- sra_path
  if (using_temp_dir) {
    message("- Copying file to temp drive..")
    temp_path <- file.path(tempdir(), accession)
    file.copy(sra_path, temp_path)
  }
  ret <- system2(sratoolkit_binary,
                 c(temp_path, "-f", progress_bar, "-t", fasterq_temp_dir,
                   threads, "--split-files", "--skip-technical",
                   "--outdir", fasterq_temp_dir)
  )
  if (using_temp_dir) {
    message("- Copying file back to main drive..")
    out_file_temp <- paste0(temp_path, if(PAIRED_END) c("_1", "_2"), ".fastq")
    stopifnot(all(file.exists(out_file_temp)))
    out_file <- paste0(sra_path, if(PAIRED_END) c("_1", "_2"), ".fastq")
    file.copy(out_file_temp, out_file, overwrite = TRUE)
    file.remove(c(out_file_temp, temp_path))
    lite_file <- paste0(temp_path, "_1", ".fastq")
    if (file.exists(lite_file)) file.remove(lite_file)
  }

  # Compress the output if needed
  if (compress) {
    filenames <- paste0(sra_path, if(PAIRED_END) c("_1", "_2"), ".fastq")
    for (filename in filenames) {
      ret <- system2("pigz", c("-f", "--best",
                               fs::path_join(c(outdir, filename))
      ))
      stopifnot("pigz exited with non-zero exit code" = ret == 0)
    }
  }
  return(ret == 0)
}

#'
#' Check 2 files for paired end etc.
#' Validate no other files exist that should not be there, try to remove
#' @noRd
cleanup_and_validate_fastq_download <- function(accession, outdir, PAIRED_END,
                                                compress, suffix = ".fastq",
                                                compression_format = ".gz") {
  if(!compress) {
    compression_format <- compression_format_regex <- ""
  } else compression_format_regex <- paste0("\\", compression_format)
  filenames_all <- filenames <- list.files(outdir, pattern = accession)
  if (length(filenames) == 0) stop("There are no files with pattern ", accession, "in directory: ",
                                   outdir)
  suffix_regex <- paste0("\\", suffix, compression_format_regex)
  filenames <- grep(suffix_regex, filenames, value = TRUE)
  if (length(filenames) == 0) stop("Accession found, but no file with valid suffix: ", suffix_regex,
                                   "Found accession files: ", paste(filenames_all, collapse = ", "))
  # TODO: Check this is correct, If special format .lite fastq, remove the larger file
  lite_file_format <- paste0(accession, c("", "_1"), suffix, compression_format)
  if (setequal(filenames, lite_file_format)) {
    lite_file_regex <- paste0(accession, "_1", suffix_regex)
    lite_file <- list.files(outdir, lite_file_regex, full.names = TRUE)
    message("Removing detected 'lite file' from SRA: ", lite_file)
    file.remove(lite_file)
    filenames <- filenames[grep(lite_file_regex, filenames, invert = TRUE)]
    stopifnot(length(filenames) == 1)
  }
  return(validate_fastq_download(filenames, accession, PAIRED_END, compress,
                                 outdir))
}

#' Check if fastq-dump returned either a single accession.fastq file or
#' two accession_(1,2).fastq files and whether it matches LibraryLayout
#' TODO: figure out what to do otherwise
#' @noRd
validate_fastq_download <- function(filenames, accession, PAIRED_END, is_compressed = FALSE,
                                    outdir, suffix = ".fastq", compression_format = ".gz") {

  if (is_compressed) {
    suffix_old <- suffix
    non_compressed_version <- paste0(accession, suffix_old)
    if (non_compressed_version %in% filenames) {
      warning("Non compressed version of ", suffix, " found in directory,
              deleting noncompressed version!")
      file.remove(file.path(outdir, non_compressed_version))
      filenames <- filenames[!(filenames %in% non_compressed_version)]
    }
    suffix <- paste0(suffix, compression_format)
  } else {
    suffix_old <- paste0(suffix, compression_format)
    compressed_version <- paste0(accession, suffix_old)
    if (compressed_version %in% filenames) {
      warning("Compressed version of ", suffix, " found in directory,
              deleting compressed version!")
      file.remove(file.path(outdir, compressed_version))
      filenames <- filenames[!(filenames %in% compressed_version)]
    }
  }
  if (setequal(filenames, paste0(accession, suffix))) {
    stopifnot(
      "extracted one fastq, but LibraryLayout is not SINGLE" =
        !PAIRED_END
    )
  } else if (setequal(filenames, sapply(
    paste0(c("_1", "_2"), suffix), function(ext) paste0(accession, ext)
  ))) {
    stopifnot(
      "extracted two fastqs, but LibraryLayout is not PAIRED" =
        PAIRED_END
    )
  } else {
    fs::file_delete(fs::dir_ls(outdir,
                               glob = paste0("**/", accession, paste0("*", suffix))
    ))
    stop(
      "Deleting .fastq suffix files, invalid filenames returned from fasterq-dump: ",
      paste(filenames, collapse = " ")
    )
  }
  return(return(filenames))
}

find_ascp_srr_url <- function(run_accession) {
  ORFik::find_url_ebi(run_accession, ebi_file_format = "sra_ftp",
                      convert_to_ascp = TRUE)
}

check_tempdir_has_space_sra_to_fastq <- function(accession, outdir, tempdir, PAIRED_END) {
  file <- fs::path_join(c(outdir, accession))
  if (!file.exists(file)) {
    file_extracted_gb <- 5 # Set a guess, TODO: make safe
  } else {
    file_extracted_gb <- 3.7*(file.size(file) / (1024^3))
  }
  file_extracted_gb <- ifelse(PAIRED_END, file_extracted_gb*2, file_extracted_gb)

  tempdir_free_gb <- get_system_usage(drive = ORFik:::detect_drive(tempdir))$Drive_Free
  tempdir_free_gb <- as.numeric(sub("[a-z]+|[A-Z]+", "", tempdir_free_gb))

  dir_to_use <- ifelse(file_extracted_gb > tempdir_free_gb, outdir, tempdir)
  if (is.na(dir_to_use)) dir_to_use <- outdir
  return(dir_to_use)
}

delete_existing_preformat_files <- function(outdir, accessions, delete_srr_preformat = TRUE) {
  preformat_files <- fs::path(outdir, accessions)
  delete_srr_preformat <- delete_srr_preformat && any(file.exists(preformat_files))
  if(delete_srr_preformat) {
    fs::file_delete(preformat_files[file.exists(preformat_files)])
  }
}

install_ascp <- function(path = ".aspera/connect/bin/ascp") {
  install_path <- fs::path_home(path)
  if (file.exists(install_path)) return(install_path)
  message("Installing ascp to default location:", install_path)
  if (.Platform$OS.type != "unix")
    stop("On windows OS, run through WSL!")
  is_linux <- Sys.info()[1] == "Linux"
  if (is_linux) {
    base_url <- "ibm-aspera-connect_4.1.0.46-linux_x86_64.tar.gz"
  } else base_url <- stop("Implement")
  tempfile <- file.path(tempdir(), "aspera.tar.gz")
  download.file(paste0("https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/0a07f/0/", base_url), tempfile)
  file <- untar(tempfile, exdir = dirname(tempfile))
  sh <- gsub(".tar.gz$", ".sh", tempfile)
  system(paste("chmod +x", sh))
  system(sh)
  message("done")
  return(install_path)
}

install_aws <- function(path = "~/bin/aws") {

  install_path <- fs::path_home(path)
  if (file.exists(install_path)) return(install_path)
  message("Installing aws to default location:", install_path)
  if (.Platform$OS.type != "unix")
    stop("On windows OS, run through WSL!")
  is_linux <- Sys.info()[1] == "Linux"
  if (is_linux) {
    base_url <- "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip"
  } else base_url <- stop("Implement")
  tempfile <- file.path(tempdir(), "aws-cli.zip")
  download.file(base_url, tempfile)
  file <- unzip(tempfile, exdir = dir)
  sh <- paste0(file.path(dir, "install"), "-i ~/bin -b ~/bin")
  system(sh)
  message("done")
  return(install_path)
}
