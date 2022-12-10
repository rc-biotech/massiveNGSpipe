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

#### Primary pipeline

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
Package is also available here on github
```r
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("rc-biotech/massive_NGS_pipe")
```  
