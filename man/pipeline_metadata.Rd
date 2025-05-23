% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_metadata.R
\name{pipeline_metadata}
\alias{pipeline_metadata}
\title{Download all metadata for all accessions}
\usage{
pipeline_metadata(
  accessions,
  config,
  force = FALSE,
  max_attempts = 7,
  max_attempts_loop = 1,
  BPPARAM = if (length(accessions) > 5) {
     bpparam()
 } else SerialParam()
)
}
\arguments{
\item{accessions}{character vector, all candidate accession numbers, allowed types:\cr
\itemize{
\item Biorpoject ID (Only numbers)\cr
\item Bioproject accession (PRJ)\cr
\item SRA study (SRA)\cr
\item ENA study (ERA)\cr
\item GEO study (GSE)\cr
\item NULL (not "NULL", but NULL), don't add new accession, only look if you
current google sheet or local csv is ready.
}}

\item{config}{path, default \code{ORFik::config()}, where will
fastq, bam, references and ORFik experiments go}

\item{force}{logical, default FALSE. Force redownload of all, if something failed}

\item{max_attempts}{numeric, default 7. How many attempts to download a specific study beforing
failing? This call flood SRA with calls, so set it higher if you use a lot
of samples. Will use a random increasing wait time between study attempts.}

\item{max_attempts_loop}{numeric, default 7. How many attempts for full download loop beforing
failing? This call flood SRA with calls, so set it higher if you use a lot
of studies. Will use a random increasing wait time between study attempts. See log files for more error info.}
}
\value{
a data.table of all metadata for all accessions
}
\description{
Download all metadata for all accessions
}
