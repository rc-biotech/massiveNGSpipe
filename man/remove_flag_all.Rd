% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_flags.R
\name{remove_flag_all}
\alias{remove_flag_all}
\title{Remove flags for all pipelines specified}
\usage{
remove_flag_all(config, steps = names(config$flag), pipelines)
}
\arguments{
\item{config}{list, output from pipeline_config(), the global config for your
NGS pipeline}

\item{steps, }{which flags to set to FALSE (not done)
for given pipeline objects}

\item{pipelines}{list, output of pipeline_init_all}
}
\value{
invisible(NULL)
}
\description{
Set all flags specified as not done for all experiments in
for all pipeline objects specified.
}
