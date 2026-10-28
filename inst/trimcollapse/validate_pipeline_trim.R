# Validate the fast trim + collapse path inside the real pipeline stage (pipe_trim_collapse), against fastp + ORFik.
# Builds a throw-away local project with the given single-end fastq files (default: two small samples made from any
# fastq you pass), runs pipe_trim_collapse() with fast_trim = FALSE and = TRUE, and compares everything later steps use:
# collapsed reads (sequence/count multisets), report summaries/filtering/adapter fields, adapter_barcode_table.csv,
# sample markers and stage flags.
# Usage: Rscript validate_pipeline_trim.R <massiveNGSpipe source dir> <fastq1> [<fastq2> ...]
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else "~/livemount/forks/massiveNGSpipe"
inputs <- if (length(args) > 1) args[-1] else stop("give one or more single-end fastq(.gz) files")
suppressMessages(devtools::load_all(src, quiet = TRUE))
library(data.table)

base <- file.path(tempdir(), "validate_fast_trim")
unlink(base, recursive = TRUE)
runs <- sprintf("VAL%02d", seq_along(inputs))
run_once <- function(fast) {
  root <- file.path(base, if (fast) "fast" else "fastp")
  raw <- file.path(root, "raw_data"); bam <- file.path(root, "processed_data"); proj <- file.path(root, "local_pipeline_val")
  study <- "validation_study"; exp <- paste0(study, "-homo_sapiens")
  dir.create(file.path(raw, exp), recursive = TRUE); dir.create(bam, recursive = TRUE); dir.create(proj, recursive = TRUE)
  for (i in seq_along(inputs)) file.copy(inputs[i], file.path(raw, exp, paste0(runs[i], ".fastq.gz")))
  conf_orfik <- c(fastq = raw, bam = bam, ref = file.path(root, "ref"), exp = file.path(root, "exp"))
  config <- pipeline_config(project_dir = proj, config = conf_orfik, mode = "local", google_url = NULL, fast_trim = fast)
  config$error_dir <- file.path(root, "error_logs"); dir.create(config$error_dir, showWarnings = FALSE)
  config$session_dir <- file.path(root, "session"); dir.create(config$session_dir, showWarnings = FALSE)
  conf <- c(exp = exp, bam = file.path(bam, exp), fastq = file.path(raw, exp), ref = file.path(root, "ref"))
  sdt <- data.table(Run = runs, ScientificName = "Homo sapiens", LIBRARYTYPE = "RFP", LibraryLayout = "SINGLE", BioProject = study)
  pipelines <- list(list(accession = study, study = sdt, organisms = list(`Homo sapiens` = list(conf = conf, index = NULL))))
  names(pipelines) <- study
  t0 <- Sys.time(); pipe_trim_collapse(pipelines, config)
  list(root = root, trim = file.path(bam, exp, "trim"), config = config, exp = exp,
       seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
read_collapsed <- function(f) {
  x <- Biostrings::readDNAStringSet(f)
  dt <- data.table(seq = as.character(x), n = as.integer(sub(".*_x", "", names(x))))
  dt[, .(n = sum(n)), by = seq][order(seq)]
}
a <- run_once(FALSE); b <- run_once(TRUE)
ok <- TRUE
check <- function(label, cond) { cat(if (cond) "PASS " else "FAIL ", label, "\n", sep = ""); ok <<- ok && cond }
for (r in runs) {
  fa <- file.path(a$trim, "SINGLE", paste0("collapsed_trimmed_", r, ".fasta.gz"))
  fb <- file.path(b$trim, "SINGLE", paste0("collapsed_trimmed_", r, ".fasta.gz"))
  check(paste(r, "collapsed file exists in trim/SINGLE (both paths)"), file.exists(fa) && file.exists(fb))
  check(paste(r, "collapsed reads identical"), isTRUE(all.equal(read_collapsed(fa), read_collapsed(fb))))
  ja <- jsonlite::fromJSON(file.path(a$trim, paste0("report_", r, ".json")))
  jb <- jsonlite::fromJSON(file.path(b$trim, paste0("report_", r, ".json")))
  check(paste(r, "report summary identical"), isTRUE(all.equal(ja$summary[c("before_filtering", "after_filtering")],
                                                               jb$summary[c("before_filtering", "after_filtering")], tolerance = 1e-4)))
  check(paste(r, "report filtering identical"), isTRUE(all.equal(ja$filtering_result[order(names(ja$filtering_result))],
                                                                 jb$filtering_result[order(names(jb$filtering_result))])))
  check(paste(r, "adapter fields identical"), identical(ja$adapter_cutting$read1_adapter_sequence, jb$adapter_cutting$read1_adapter_sequence) &&
          identical(ja$adapter_cutting$adapter_trimmed_reads, jb$adapter_cutting$adapter_trimmed_reads))
  check(paste(r, "no trimmed fastq left by the fast path"), !file.exists(file.path(b$trim, paste0("trimmed_", r, ".fastq"))))
}
ta <- fread(file.path(a$trim, "adapter_barcode_table.csv")); tb <- fread(file.path(b$trim, "adapter_barcode_table.csv"))
cols <- c("id", "adapter", "barcode_detected", "barcode5p_size", "barcode3p_size", "mean_length_adapter_filtered")
check("adapter_barcode_table identical (adapter, barcode sizes, lengths)", isTRUE(all.equal(ta[, ..cols], tb[, ..cols])))
check("trim and collapsed markers + flags set (fast path)",
      setequal(samples_done(b$config, "trim", b$exp), runs) && setequal(samples_done(b$config, "collapsed", b$exp), runs) &&
        step_is_done(b$config, "trim", b$exp) && step_is_done(b$config, "collapsed", b$exp))
print(ta[, ..cols])
cat(sprintf("pipe_trim_collapse wall time: fastp + ORFik %.0f s, fast path %.0f s\n", a$seconds, b$seconds))
cat(if (ok) "ALL PASS\n" else "SOME FAILED\n")
