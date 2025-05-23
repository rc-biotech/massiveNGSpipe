% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/google_API.R
\name{plot_all_versions}
\alias{plot_all_versions}
\title{Export ggplots to multiple platforms}
\usage{
plot_all_versions(
  plot,
  file_prefix,
  formats = c("jpg", "svg"),
  send_to_google_drive = FALSE,
  send_to_discord = FALSE,
  width = 8,
  height = 6,
  dpi = 600,
  discord_connection = discord_connection_default_cached(),
  google_drive_dir = google_drive_dir_links(1),
  formats_discord = formats[!(formats \%in\% "svg")],
  formats_google = formats[!(formats \%in\% c("jpg", "png"))],
  preview_image = FALSE,
  discord_message = NULL
)
}
\arguments{
\item{plot}{the ggplot}

\item{file_prefix}{the name of plot without the file extension}

\item{formats}{c("jpg", "svg")}

\item{send_to_google_drive}{= FALSE,}

\item{send_to_discord}{= FALSE}

\item{width}{= 8}

\item{height}{= 6}

\item{dpi}{= 600}

\item{discord_connection}{= discord_connection_default_cached()}

\item{google_drive_dir}{= google_drive_dir_links(1)}

\item{formats_discord}{= formats\link{!(formats \%in\% "svg")}}

\item{formats_google}{= formats\link{!(formats \%in\% c("jpg", "png"))}}

\item{preview_image}{= FALSE, if TRUE, will open browser and display image, then wait for you to press enter,
then you can cancel if ratios are wrong etc.}

\item{discord_message}{= NULL, if character, sends the message to discord after images are sent.}
}
\value{
invisible(NULL)
}
\description{
Supports: Disc, google drive and discord
}
