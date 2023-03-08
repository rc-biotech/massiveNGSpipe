# massive_NGS_pipe

#### About

massive_NGS_pipe is a R package for full process integration of
NGS studies, it can handle thousands of samples from multiple organisms
in a single run.

It currently supports to set up easily a fail safe and fallback supported 
pipeline that does the specific steps:
- Find candidate accessions for studies you might be interested in
- A clean way of helping with annotating samples using
either google sheet integration or local csv files.
- Automatic download of genome and gff for each organism (supports
fixing malformed gffs, adding pseudo 5' UTRs etc)
- Automatic download of samples with fallback options:
 1. AWS (amazone, fastest), 2. ENA (fast), 3. fastq-dump (slowest)
- Automatic detection of 3' adapters and trimming (fastp)
- Collapsing of duplicated reads for faster processing
- barcode / UMI detection / removal (To be implemented)
- Autmoatic removal of contaminants: phix, rRNA, tRNA,..
- Genome alignment with STAR (with QC)
- Create ORFik.experiment object for easy analysis in R
- Convert output files to any format: bigwig, wig, bed, ofst etc.
- Merging of samples by type

Ribo-seq specific:
- automatic pshifting and quality validation

Visulization specific:
Our data output is directly supported to be browsed by the
RiboCrypt visuzliation tool which supports among others:
- Both genomic and transciptomic view browser
- Differential expression
- Read length heatmaps
- Codon usage
- And much more

#### Primary pipeline (Data from online repository)

Package is also available here on github
```r
library(massiveNGSpipe)
## Setup
config <- pipeline_config() # <- set up paths
## Curate metadata
accessions <- accessions_to_use("GSE152850", "Saccharomyces cerevisiae", FALSE)
curate_metadata(accessions, config, "Saccharomyces cerevisiae")
## Init pipeline and show current progress report
pipelines <- pipeline_init_all(config) # Initialize pipeline configuration for all experiments
progress_report(pipelines, config)
# Start the pipeline (In RStudio start this step with the 'background jobs' tab)
run_pipeline(pipelines, config, wait = 100, BPPARAM = bpparam())
```  

#### Installation

Package is currently only available here on github
```r
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("rc-biotech/massive_NGS_pipe")
```  

#### Local pipeline (Data from personal hard drive / local server)

To run your local data, you run the pipeline, but skips the steps of 
downloading. You also need to annotate the basic sample metadata yourself,
since it does not exist online yet. 

```r
library(massiveNGSpipe)
config_dirs <- ORFik::config() # Where all data is stored (all your pipelines)
# -> project_dir specifies log, status and metadata folder for specific pipeline
project_dir <- file.path(dirname(config_dirs)[1], "test_pipeline")
mode <- "local" # Data already exists on your computer / server
google_url <- NULL # Do not use google sheet for metadata storage
preset <- "Ribo-seq" # What kind of pipeline (use Ribo-seq preset)

# Init config (set up paths and pipeline settings)
config <- pipeline_config(project_dir, config_dirs, mode = mode,
                          preset = preset, google_url = google_url)

# Now, for local data you need to define a data.table with local file metadata
# Do this for each experiment (split experiment for each organism)
name <- "local_mNGSp_run_test" # This is the folder name in the fastq folder
organism <- "Homo sapiens"
organism_caps <- gsub(" ", "_", trimws(tolower(organism)))
path <-  paste(name, organism_caps, sep = "-")
path_full <- file.path(config$config["fastq"], path)
stopifnot(basename(path) == path) # relative path only!
stopifnot(dir.exists(path_full)) # Must exist!

save_file <- file.path(config[["metadata"]],
                       paste0("SraRunInfo_", name, ".csv"))
# Add all information needed here
dt <- local_study_csv(path_full, name, organism = organism,
                      paired = "SINGLE", AUTHOR = "Håkon",
                      sample_title = c("WT_ribo_cycloheximide_rep1",
                                       "WT_ribo_cycloheximide_rep2"))
dir.create(dirname(save_file), FALSE, TRUE)
fwrite(dt, save_file)

# Now for all experimental directories, validate metadata
all_names <- c(name) # Add all here (here we did only 1)
# Now run this, when you get to validation either open your csv and set
# the KEEP column to TRUE, for sample you want to run.
# Or do as bellow to run and update KEEP column in R
curate_metadata(all_names, config)

# R way of setting samples to run.
temp <- fread(config$temp_metadata)
temp$KEEP <- TRUE
fwrite(temp, config$complete_metadata)
curate_metadata(all_names, config)

# Now init pipeline objects and run pipeline
pipelines <- pipeline_init_all(config)

run_pipeline(pipelines, config, wait = 20, BPPARAM = SerialParam())
```  
