% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_flags.R
\name{libtype_flags}
\alias{libtype_flags}
\title{Define the library specific steps}
\usage{
libtype_flags(preset, mode = "online", contam = FALSE)
}
\arguments{
\item{preset}{character, default "Ribo-seq".
Alternatives: c("Ribo-seq", "RNA-seq", "disome", "empty")}

\item{mode}{either of c("online", "local"). Local will disable some steps
like the fetch step (as data is already presumed to exist locally etc)}

\item{contam}{logical, FALSE, do contamint removal, using a seperat made
STAR index during genome preparation step.}
}
\value{
a character vector with names being the grouping in functions
}
\description{
This is called flags
}
