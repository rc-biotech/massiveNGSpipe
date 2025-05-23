% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_progress.R
\name{progress_report}
\alias{progress_report}
\title{Report summary of current progress}
\usage{
progress_report(
  pipelines,
  config,
  show_status_per_exp = TRUE,
  show_done = TRUE,
  status_plot = FALSE,
  show_stats = FALSE,
  return_progress_vector = FALSE,
  check_merged_org = FALSE,
  named_progress_vector = FALSE,
  system_usage_stats = TRUE
)
}
\arguments{
\item{pipelines}{list, output of pipeline_init_all}

\item{config}{list, output from pipeline_config(), the global config for your
NGS pipeline}

\item{show_status_per_exp}{logical, default TRUE.
Show a message per experiment with progression status.}

\item{show_done}{logical, default TRUE. If FALSE, display only status
of projects that are not done. Stats will still show for all.}

\item{status_plot}{plot an  interactive plot of total status}

\item{show_stats}{logical, default TRUE, output trim/alignment stats
\itemize{
\item plots, set to FALSE if you only want progress report.
}}

\item{return_progress_vector}{logical, default FALSE. If true,
return integer of current flag step per experiment, set
named_progress_vector to TRUE, to get experiment names as names.}

\item{check_merged_org}{logical, default FALSE. Check status for merged organisms file
relative to pipeline.}

\item{named_progress_vector}{logical, default FALSE. Add experiment names as names
to progress_vector if returned, ignored otherwise.}

\item{system_usage_stats}{logical, default TRUE. Show system usage stats.}
}
\value{
invisible(NULL) / or progress vector
}
\description{
How many studies have finished etc
}
