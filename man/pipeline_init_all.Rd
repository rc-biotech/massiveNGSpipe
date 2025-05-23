% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_init_helpers.R
\name{pipeline_init_all}
\alias{pipeline_init_all}
\title{Init all studies to pipeline objects}
\usage{
pipeline_init_all(
  config,
  complete_metadata = config$complete_metadata,
  progress_report = TRUE,
  gene_symbols = TRUE,
  only_complete_genomes = FALSE,
  dbs = c("ensembl", "refseq"),
  show_status_per_exp = TRUE
)
}
\arguments{
\item{config}{path, default \code{ORFik::config()}, where will
fastq, bam, references and ORFik experiments go}

\item{complete_metadata}{path, default file.path(project_dir, "RFP_FINAL_LIST.csv").
The list of final candidates that are checked and have unique rows per bioproject}

\item{only_complete_genomes}{logical, default FALSE. If TRUE, will only init the subset
with complete genome/annotation directories. Will fail if 0 species are done.}

\item{dbs}{character, prioritized order of genome databases to search, default:
c("ensembl", "refseq"). Options are c("ensembl", "refseq", "genbank") in any wanted
priority order.}

\item{simple_progress_report}{logical, default TRUE. Display current progress.}
}
\value{
a list of pipelines
}
\description{
Init all studies to pipeline objects
}
