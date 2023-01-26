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
                             fs::path_dir(ORFik::install.sratoolkit())) {
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

    fs::file_delete(
      fs::path(outdir, info$Run)
    )

    message("Extracting SRA runs:")
}

download_sra_aws <- function(run_accession, outdir) {
  ret <- system2("aws", c(
    "--no-sign-request", "s3", "ls",
    paste("s3://sra-pub-run-odp/sra", run_accession, sep = "/")
  ), stdout = FALSE)

  if (ret == 0) {
    ret <- system2("aws", c(
      "--no-sign-request", "s3", "sync",
      paste("s3://sra-pub-run-odp/sra", run_accession, sep = "/"),
      outdir
    ))
    stopifnot("Run found on aws, but awscli exited with non-zero exit code" = ret == 0)
  }
  return(ret == 0)
}

download_sra_ascp <- function(run_accession, outdir) {
  # If not, fall back to EBI servers via Aspera.
  # Note: accessions longer than 6 digits are grouped into another layer
  # of directories.
  vol1_path <-
    if (nchar(run_accession) == 9) {
      paste("vol1/srr",
            stringr::str_sub(run_accession, 1, 6),
            run_accession,
            sep = "/"
      )
    } else {
      paste("vol1/srr",
            stringr::str_sub(run_accession, 1, 6),
            paste0("00", stringr::str_sub(run_accession, -1)),
            run_accession,
            sep = "/"
      )
    }
  srr_path <- fs::path_join(c(outdir, run_accession))
  ret <- system2(fs::path_home(".aspera/connect/bin/ascp"), c(
    "-k1", "-QT", "-l", "300m", "-P33001",
    "-i", "~/.aspera/connect/etc/asperaweb_id_dsa.openssh",
    paste0("era-fasp@fasp.sra.ebi.ac.uk:", vol1_path),
    srr_path
  ))
  stopifnot("ascp exited with non-zero exit code" = ret == 0)
}

extract_srr_preformat <- function(accession, outdir, compress = TRUE,
                                  sratoolkit_path =
                                    fs::path_dir(ORFik::install.sratoolkit()),
                                  PAIRED_END) {

  message(accession)
  srr_path <- fs::path_join(c(outdir, accession))
  ret <- system2(fs::path_join(c(sratoolkit_path, "fasterq-dump")), c(
    "-f", "-p", "--split-files", "--skip-technical",
    "--outdir", outdir, srr_path
  ))
  stopifnot("fasterq-dump exited with non-zero exit code" = ret == 0)

  # Check if fastq-dump returned either a single <accession>.fastq file or
  # two <accession>_{1,2}.fastq files and whether it matches LibraryLayout
  # TODO: figure out what to do otherwise
  filenames <- fs::path_file(fs::dir_ls(outdir,
                                        glob = paste0("**/", accession, "*.fastq")
  ))
  if (setequal(filenames, paste0(accession, ".fastq"))) {
    stopifnot(
      "extracted one fastq, but LibraryLayout is not SINGLE" =
        !PAIRED_END
    )
  } else if (setequal(filenames, sapply(
    c("_1.fastq", "_2.fastq"), function(ext) paste0(accession, ext)
  ))) {
    stopifnot(
      "extracted two fastqs, but LibraryLayout is not PAIRED" =
        PAIRED_END
    )
  } else {
    fs::file_delete(fs::dir_ls(outdir,
                               glob = paste0("**/", accession, "*.fastq")
    ))
    stop(
      "invalid filenames returned from fasterq-dump: ",
      paste(filenames, collapse = " ")
    )
  }

  # Compress the output if needed
  if (compress) {
    for (filename in filenames) {
      ret <- system2("pigz", c("-f", "--best",
                               fs::path_join(c(outdir, filename))
      ))
      stopifnot("pigz exited with non-zero exit code" = ret == 0)
    }
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
