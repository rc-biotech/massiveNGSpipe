download_raw_srr <- function(run_accession, outdir, compress = TRUE,
                             sratoolkit_path =
                                 fs::path_dir(ORFik::install.sratoolkit()),
                             PAIRED_END) {
    # Check if the run is available on AWS (not always the case, e.g.
    # SRR4000302)
     if (!download_sra_aws(run_accession, outdir)) {
       download_sra_ascp(run_accession, outdir)
     }
  extract_srr_preformat(run_accession, outdir, compress,
                        sratoolkit_path, PAIRED_END)
}

download_sra <- function(info, outdir, compress = TRUE,
                         sratoolkit_path =
                             fs::path_dir(ORFik::install.sratoolkit()),
                         delete_srr_preformat = TRUE) {
    fs::dir_create(outdir)
    outdir <- fs::path_real(outdir)

    run_accessions <- info$Run
    if (is.null(run_accessions) | (length(run_accessions) == 0)) {
        stop("Could not find SRR numbers in 'info'")
    }

    message("Downloading SRA runs:")
    for (accession in run_accessions) {
      PAIRED_END <- info$LibraryLayout[info$Run == accession] == "PAIRED"
      download_raw_srr(accession, outdir, compress, sratoolkit_path, PAIRED_END)
    }

    preformat_files <- fs::path(outdir, info$Run)
    delete_srr_preformat <- delete_srr_preformat & any(file.exists(preformat_files))
    if(delete_srr_preformat) {
      fs::file_delete(preformat_files[file.exists(preformat_files)])
    }

    message("Extracting SRA runs:")
}

download_sra_aws <- function(run_accession, outdir, aws_bin = "aws",
                             aws_mirror = "s3://sra-pub-run-odp/sra") {

  file_aws_url <- file.path(aws_mirror, run_accession)
  ret <- system2(aws_bin, c(
    "--no-sign-request", "s3", "ls",
    file_aws_url
  ), stdout = FALSE)

  if (ret == 0) {
    ret <- system2(aws_bin, c(
      "--no-sign-request", "s3", "sync",
      file_aws_url,
      outdir
    ))
    stopifnot("Run found on aws, but awscli exited with non-zero exit code" = ret == 0)
    # TODO: Add file is downloaded check too (this is not enough)
  }
  return(ret == 0)
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

find_ascp_srr_url <- function(run_accession) {
  ORFik::find_url_ebi(run_accession, ebi_file_format = "sra_ftp",
                      convert_to_ascp = TRUE)
}

extract_srr_preformat <- function(accession, outdir, compress = TRUE,
                                  sratoolkit_path = fs::path_dir(ORFik::install.sratoolkit()),
                                  PAIRED_END, progress_bar = TRUE,
                                  fasterq_temp_dir = tempdir()) {

  message(accession)
  binary <- "fasterq-dump"
  already_compressed <- FALSE
  sratoolkit_binary <- fs::path_join(c(sratoolkit_path, binary))
  stopifnot(file.exists(sratoolkit_binary))

  srr_path <- fs::path_join(c(outdir, accession))
  progress_bar <- ifelse(progress_bar, "-p", "")
  ret <- system2(sratoolkit_binary, c(
    "-f", progress_bar, "-t", fasterq_temp_dir,
    "--split-files", "--skip-technical",
    "--outdir", outdir,
    srr_path
  ))
  fasterq_dump_failed <- ret != 0
  if (fasterq_dump_failed) {
    message("Using ebi fallback")
    res <- try(ORFik::download.SRA(accession, outdir, rename = FALSE))
    already_compressed <- TRUE
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

  filenames <- fs::path_file(fs::dir_ls(outdir,
                                        glob = paste0("**/", accession, "*.fastq")
  ))
  # TODO: Check this is correct, If special format .lite fastq, remove the larger file
  lite_file_format <- paste0(accession, c("", "_1"), ".fastq")
  if (setequal(filenames, lite_file_format)) {
    lite_file <- list.files(outdir, paste0(accession, "_1.fastq$"), full.names = TRUE)
    file.remove(lite_file)
    filenames <- filenames[grep("_1\\.fastq", filenames, invert = TRUE)]
    stopifnot(length(filenames) == 1)
  }

  validate_fastq_download(filenames, accession, PAIRED_END, already_compressed,
                          outdir)
  # Compress the output if needed
  if (compress & !already_compressed) {
    for (filename in filenames) {
      ret <- system2("pigz", c("-f", "--best",
                               fs::path_join(c(outdir, filename))
      ))
      stopifnot("pigz exited with non-zero exit code" = ret == 0)
    }
  }

  return(filenames)
}

#' Check if fastq-dump returned either a single accession.fastq file or
#' two accession_(1,2).fastq files and whether it matches LibraryLayout
#' TODO: figure out what to do otherwise
#' @noRd
validate_fastq_download <- function(filenames, accession, PAIRED_END, is_compressed = FALSE,
                                    outdir) {
  suffix <- ".fastq"
  if (is_compressed) suffix <- paste0(suffix, ".gz")
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
      "invalid filenames returned from fasterq-dump: ",
      paste(filenames, collapse = " ")
    )
  }
  return(invisible(NULL))
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
