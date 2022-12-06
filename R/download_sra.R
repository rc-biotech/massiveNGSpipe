download_raw_srr <- function(run_accession, outdir,
                             sratoolkit_path =
                                 fs::path_dir(ORFik::install.sratoolkit())) {
    srr_path <- fs::path_join(c(outdir, run_accession))

    # Check if the run is available on AWS (not always the case, e.g.
    # SRR4000302)
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
        stopifnot("awscli exited with non-zero exit code" = ret == 0)
    } else {
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
        ret <- system2(fs::path_home(".aspera/connect/bin/ascp"), c(
            "-k1", "-QT", "-l", "300m", "-P33001",
            "-i", "~/.aspera/connect/etc/asperaweb_id_dsa.openssh",
            paste0("era-fasp@fasp.sra.ebi.ac.uk:", vol1_path),
            srr_path
        ))
        stopifnot("ascp exited with non-zero exit code" = ret == 0)
    }
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
        download_raw_srr(accession, outdir, sratoolkit_path)
    }

    message("Extracting SRA runs:")
    for (accession in run_accessions) {
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
                info$LibraryLayout[info$Run == accession] == "SINGLE"
            )
        } else if (setequal(filenames, sapply(
            c("_1.fastq", "_2.fastq"), function(ext) paste0(accession, ext)
        ))) {
            stopifnot(
                "extracted two fastqs, but LibraryLayout is not PAIRED" =
                info$LibraryLayout[info$Run == accession] == "PAIRED"
            )
        } else {
            fs::file_delete(fs::dir_ls(outdir,
                glob = paste0("**/", run_accession, "*.fastq")
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
}