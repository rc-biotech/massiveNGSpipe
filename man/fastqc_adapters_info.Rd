% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastq_helpers.R
\name{fastqc_adapters_info}
\alias{fastqc_adapters_info}
\title{Helper function to find the main adapter in a .fastq.gz file}
\usage{
fastqc_adapters_info(
  file,
  nreads = 2e+06,
  adapters_file = fs::path(tempdir(), "adapter_candidates.txt")
)
}
\arguments{
\item{file}{full path to fastq file to check with qc. Supports both
".gz" and non compressed.}

\item{nreads}{integer, default 6000000 (6e6),
number of reads to read in for check}

\item{adapters_file}{= fs::path(tempdir(), "adapter_candidates.txt"),
will use a internal set of known adapters, saved to a tempfile.}
}
\value{
adapter candidates
}
\description{
Helper function to find the main adapter in a .fastq.gz file
}
\examples{
# Define parameters
n_strings <- 500000  # Number of strings
fixed_total_size <- 51
string_length <- seq(20, 33)  # Length of each string
string_lengths <- sample(string_length, n_strings, replace = TRUE)  # Length of each string
letters <- c("A", "T", "C", "G")  # DNA letters
# Generate random DNA strings
all_chars <- DNAString(paste(sample(letters, sum(string_lengths), replace = TRUE), collapse = ""))
dna_strings <- DNAStringSet(successiveViews(all_chars, width = string_lengths))
# View the first few strings
head(dna_strings)

# Make barcode, dna, adapter
illumina_true_seq <- DNAStringSet(DNAString("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"))
barcode5p <- DNAStringSet(DNAString("AATTGG"))
barcode3p <- DNAStringSet(DNAString("AATTGG"))
dna <- DNAStringSet(paste0(barcode5p, dna_strings, barcode3p, illumina_true_seq))
dna <- subseq(dna, 1, fixed_total_size)
names(dna) <- paste0("seq_", seq(length(dna)))
qualities <- substr("#<GGGGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFGGGGGGGGG", 1, fixed_total_size)
mcols(dna)$qualities <- rep(BStringSet(qualities), length(dna))

tempfile <- tempfile(fileext = ".fastq")
writeXStringSet(dna, tempfile, format = "fastq")

fastqc_adapters_info(tempfile, 1e6)
}
