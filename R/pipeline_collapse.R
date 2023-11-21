#' Collapse trimmed single-end reads and move them into "trim/SINGLE"
#' subdirectory. PE reads are moved into "trim/PAIRED" without modification.
#' @param pipeline a pipeline object
pipe_collapse <- function(pipeline, config) {
  study <- pipeline$study
  for (organism in names(pipeline$organisms)) {
    conf <- pipeline$organisms[[organism]]$conf
    if (!step_is_next_not_done(config, "collapsed", conf["exp"])) next
    trimmed_dir <- fs::path(conf["bam"], "trim")
    runs_paired <- study[ScientificName == organism &
                           LibraryLayout == "PAIRED"]
    runs_single <- study[ScientificName == organism &
                           LibraryLayout != "PAIRED"]
    any_paired_libs <- nrow(runs_paired) > 0
    if (any_paired_libs) {
      outdir <- fs::path(trimmed_dir, runs_paired[1]$LibraryLayout)
      fs::dir_create(outdir)
      filenames <- paste0(
        "trimmed" , c("1_", "2_"), runs_paired[i]$Run, "_1", ".fastq"
      )
      message("Paired end data is now collapsed into 1 file,
                    and 2nd file is reverse complimented before merging!")
      full_filenames <- file.path(trimmed_dir, filenames)
      BiocParallel::bplapply(full_filenames, function(filename) {
        ORFik::collapse.fastq(
          filename, outdir,
          compress = TRUE
        )
        fs::file_delete(filename)
      }, BPPARAM = BiocParallel::MulticoreParam(16))
      # Read in and reverse second file, then merge back into 1.
      full_filenames <- file.path(outdir, paste0("collapsed_", filenames))
      first_files <- filenames[seq_along(filenames) %% 2 == 1]
      second_files <- filenames[seq_along(filenames) %% 2 == 0]

      for (i in seq_along(first_files)) {
        first_file <- first_files[(i*2)-1]
        second_file <- second_files[(i*2)]
        a <- readDNAStringSet(first_file, format = "fasta", use.names = TRUE)
        b <- readDNAStringSet(second_file, format = "fasta", use.names = TRUE)
        b <- reverseComplement(b)
        writeXStringSet(c(a,b), first_file, format = "fasta")
        fs::file_delete(second_file)
      }
    }

    any_single_libs <- nrow(runs_single) > 0
    if (any_single_libs) {
      outdir <- fs::path(trimmed_dir, "SINGLE")
      fs::dir_create(outdir)
      BiocParallel::bplapply(runs_single$Run, function(srr) {
        filename <- list.files(trimmed_dir, paste0(srr, "\\."),
                               full.names = TRUE)
        if (length(filename) != 1) {
          filename <- file.path(trimmed_dir,
                                paste0("trimmed_", srr, ".fastq"))
        }
        ORFik::collapse.fastq(
          filename, outdir,
          compress = TRUE
        )
        fs::file_delete(filename)
      }, BPPARAM = BiocParallel::MulticoreParam(16))
    }

    set_flag(config, "collapsed", conf["exp"])
  }
}

pipe_cigar_collapse_single <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "cigar_collapse", name(df))) next
    files <- filepath(df, "ofst")
    for (f in files) {
      ofst <- fimport(f)
      ORFik:::export.ofst(collapseDuplicatedReads(ofst), f)
    }
    set_flag(config, "cigar_collapse", name(df))
  }
}

pipe_cigar_collapse <- function(pipelines, config) {
  exp <- get_experiment_names(pipelines)
  done_exp <- unlist(lapply(exp, function(e) step_is_next_not_done(config, "cigar_collapse", e)))

  for (experiments in exp[done_exp]) {
    try <- try({
      df_list <- lapply(experiments, function(e)
        read.experiment(e, validate = FALSE, output.env = new.env()))
      pipe_cigar_collapse_single(df_list, config)
    })
    if (is(try, "try-error"))
      warning("Failed at step, cigar_collapse, study: ", experiments[1])
  }
}
