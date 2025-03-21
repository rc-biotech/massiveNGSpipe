% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_run.R
\name{pipelines_names}
\alias{pipelines_names}
\title{Pipeline experiment names}
\usage{
pipelines_names(pipelines, recursive = TRUE)
}
\arguments{
\item{pipelines}{a list, the pipelines object}

\item{recursive}{logical, default TRUE, If false return as list}
}
\value{
character vector of experiment names, recursive TRUE gives list.
}
\description{
Get all pipeline experiment names from pipeline
Remember all species per pipeline object is combined, so this
function unlists the whole when recursive is TRUE
}
