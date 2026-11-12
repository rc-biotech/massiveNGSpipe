#' Add a reference genome from files or an NCBI assembly accession
#'
#' Builds a complete pipeline reference (as \code{get_annotation_and_index()} does for Ensembl/RefSeq downloads) from a
#' genome fasta and a GTF/GFF annotation that you supply: local paths or URLs (gzip is handled), or an NCBI assembly
#' accession (GCA_/GCF_), whose \code{_genomic.fna.gz} and \code{_genomic.gtf.gz} (else \code{_genomic.gff.gz}) are
#' fetched from the NCBI FTP site. Use it for assemblies biomartr cannot find (new GenBank-only assemblies such as
#' T2T-goat2.0), lab assemblies, or annotation you curated yourself.\cr
#' Steps: files are copied/downloaded into \code{<ref_dir>/<reference_folder_name(organism)>/}; ORFik builds the
#' optimised annotation (txdb, regions), contaminants (phiX, ncRNA, tRNA, rRNA) and \code{outputs.rds};
#' \code{assembly_fragmentation_flag()} is checked; the STAR index is built. \code{reference_source.rds} records the
#' sources, md5 sums and date. Afterwards \code{pipeline_init_all()} finds a complete reference and uses it like any other,
#' and \code{make_mixed_genome()} can combine it with other species.
#' @param organism scientific name, e.g. "Capra hircus"
#' @param genome path or URL of the genome fasta (.fa/.fna/.fasta, optionally .gz); NULL when \code{accession} is given
#' @param annotation path or URL of the GTF/GFF3 (optionally .gz); NULL when \code{accession} is given
#' @param accession NCBI assembly accession (e.g. "GCA_041735815.1"), used when genome/annotation are NULL. Default: the
#' assembly pinned for the organism in \code{reference_assembly()} (study strain), when there is one
#' @param ref_dir reference directory, default \code{ORFik::config()["ref"]}
#' @param star_index logical, default TRUE: build the STAR index
#' @param overwrite logical, default FALSE: stop if the folder already holds a complete reference (outputs.rds)
#' @return the annotation paths (named character, as \code{ORFik::getGenomeAndAnnotation()}), invisibly
#' @export
#' @examples
#' # add_custom_genome("Capra hircus", accession = "GCA_041735815.1")
#' # add_custom_genome("Mus musculus", genome = "~/my_mm39.fa", annotation = "~/my_annotation.gtf")
add_custom_genome <- function(organism, genome = NULL, annotation = NULL, accession = pinned_accession(organism),
                              ref_dir = ORFik::config()["ref"], star_index = TRUE, overwrite = FALSE) {
  if (is.null(genome) != is.null(annotation)) stop("Give both genome and annotation, or an accession")
  if (!is.null(genome) && missing(accession)) accession <- NULL   # own files: the pinned accession does not describe them
  if (is.null(genome)) {
    if (is.null(accession)) stop("Give genome + annotation, or an NCBI accession")
    urls <- ncbi_assembly_urls(accession)
    genome <- urls[["genome"]]; annotation <- urls[["annotation"]]
  }
  dir <- file.path(ref_dir, reference_folder_name(organism))
  if (file.exists(file.path(dir, "outputs.rds")) && !overwrite)
    stop("Complete reference exists already: ", dir, " (move it away or use overwrite = TRUE)")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  stem <- gsub(" ", "_", organism)
  label <- if (!is.null(accession)) accession else "custom"
  fa <- custom_reference_fetch(genome, file.path(dir, paste0(stem, "_", label, "_genomic.fa")))
  ext <- if (grepl("\\.gff3?(\\.gz)?$", annotation, ignore.case = TRUE)) "gff3" else "gtf"
  gtf <- custom_reference_fetch(annotation, file.path(dir, paste0(stem, "_", label, "_genomic.", ext)))
  message("- Building reference for ", organism, " in ", dir)
  db <- if (!is.null(accession) && startsWith(accession, "GCF_")) "refseq" else "genbank"   # steers ORFik lookups only
  ann <- get_annotation(organism, output.dir = dir, genome = fa, GTF = gtf, db = db, gene_symbols = FALSE,
                        assembly = if (is.null(accession)) organism else accession)
  if (is.null(ann)) stop("ORFik could not build the reference from ", fa, " and ", gtf)
  assembly_fragmentation_flag(ann["genome"])
  saveRDS(list(organism = organism, accession = accession, genome_source = genome, annotation_source = annotation,
               genome_md5 = unname(tools::md5sum(fa)), annotation_md5 = unname(tools::md5sum(gtf)),
               created = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file.path(dir, "reference_source.rds"))
  if (star_index)
    ORFik::STAR.index(ann, notify_load_existing = FALSE, max.ram = file.size(ann["genome"]) / 1024^3 * 20)
  invisible(ann)
}

#' Assembly accession pinned for an organism in reference_assembly(), else NULL
#' @noRd
pinned_accession <- function(organism) {
  a <- reference_assembly(organism)
  if (grepl("^GC[AF]_[0-9]{9}\\.[0-9]+$", a)) a else NULL
}

#' NCBI FTP URLs of an assembly's genome fasta and annotation
#'
#' Resolves the assembly directory (\code{genomes/all/GCA/041/735/815/GCA_041735815.1_<name>/}) and returns the
#' \code{_genomic.fna.gz} and \code{_genomic.gtf.gz} URLs (\code{_genomic.gff.gz} when there is no GTF).
#' @param accession "GCA_" or "GCF_" accession with version
#' @return named character: genome, annotation
#' @export
ncbi_assembly_urls <- function(accession) {
  if (!grepl("^GC[AF]_[0-9]{9}\\.[0-9]+$", accession)) stop("Not an NCBI assembly accession with version: ", accession)
  digits <- sub("^GC[AF]_([0-9]{9})\\..*$", "\\1", accession)
  base <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/", substr(accession, 1, 3), "/",
                 substr(digits, 1, 3), "/", substr(digits, 4, 6), "/", substr(digits, 7, 9), "/")
  listing <- tryCatch(readLines(base, warn = FALSE), error = function(e) stop("Cannot list ", base, ": ", conditionMessage(e)))
  dirs <- unique(regmatches(listing, regexpr(paste0(accession, "_[^/\"]+"), listing)))
  if (length(dirs) != 1) stop("Assembly directory for ", accession, " not found under ", base)
  files <- tryCatch(readLines(paste0(base, dirs, "/"), warn = FALSE), error = function(e) character())
  has <- function(suffix) any(grepl(paste0(dirs, suffix), files, fixed = TRUE))
  ann <- if (has("_genomic.gtf.gz")) "_genomic.gtf.gz" else if (has("_genomic.gff.gz")) "_genomic.gff.gz" else
    stop(accession, " has no GTF/GFF annotation on the NCBI FTP site")
  c(genome = paste0(base, dirs, "/", dirs, "_genomic.fna.gz"), annotation = paste0(base, dirs, "/", dirs, ann))
}

#' Copy or download one reference file to its destination, decompressing gzip
#' @noRd
custom_reference_fetch <- function(source, dest) {
  if (file.exists(dest)) return(dest)
  gz <- grepl("\\.gz$", source)
  tmp <- paste0(dest, if (gz) ".gz" else "", ".part")
  if (grepl("^(https?|ftp)://", source)) {
    # R's default timeout (60 s) cuts large genomes (gorilla T2T: 974 MB, stopped at 738 MB, 2026-11-12)
    old_to <- options(timeout = max(7200, getOption("timeout"))); on.exit(options(old_to), add = TRUE)
    if (utils::download.file(source, tmp, mode = "wb", quiet = TRUE) != 0) stop("Download failed: ", source)
  } else {
    if (!file.exists(source)) stop("File not found: ", source)
    stopifnot(file.copy(source, tmp, overwrite = TRUE))
  }
  if (gz) {
    R.utils::gunzip(tmp, paste0(dest, ".part2"), remove = TRUE)
    stopifnot(file.rename(paste0(dest, ".part2"), dest))
  } else stopifnot(file.rename(tmp, dest))
  dest
}
