# Optional, source-reviewed anchored sequence cleanup AFTER fixed front UMIs and
# adapter/post-adapter-tail removal. Unconfigured studies retain their old path.
anchored_front_settings <- function(trimmed_dir, run, paired = FALSE) {
  path <- file.path(trimmed_dir, 'anchored_front_manual.csv')
  if (!file.exists(path)) return(NULL)
  x <- data.table::fread(path)
  needed <- c('Run', 'base', 'min', 'max')
  if (!all(needed %in% names(x)) || anyDuplicated(x$Run))
    stop('Anchored front table requires unique Run,base,min,max columns')
  z <- x[Run == run]
  if (!nrow(z)) return(NULL)
  if (paired) stop('Anchored front cleanup currently requires single-end footprint input; paired handling must be explicit')
  if (nrow(z) != 1L || is.na(z$base) || !z$base %in% c('A','C','G','T') ||
      is.na(z$min) || is.na(z$max) || !is.numeric(z$min) || !is.numeric(z$max) ||
      z$min != as.integer(z$min) || z$max != as.integer(z$max) || z$min < 1 ||
      z$max < z$min || z$max > 50)
    stop('Invalid anchored front base/min/max contract')
  list(base = z$base, min = as.integer(z$min), max = as.integer(z$max))
}

# Separate helper makes the exact longest-anchored-match rule directly testable.
anchored_front_chunk <- function(fq, base, min, max) {
  sequence <- ShortRead::sread(fq)
  first <- as.character(Biostrings::subseq(sequence, start = 1L,
                        width = pmin(Biostrings::width(sequence), max)))
  cut <- attr(regexpr(paste0('^',base,'{',min,',',max,'}'), first), 'match.length')
  cut <- pmax(0L, as.integer(cut))
  # BStringSet retains the original ASCII qualities without a Phred round trip.
  quality <- Biostrings::quality(Biostrings::quality(fq))
  output <- ShortRead::ShortReadQ(
    sread = Biostrings::subseq(sequence, start = cut + 1L),
    quality = Biostrings::subseq(quality, start = cut + 1L), id = ShortRead::id(fq))
  list(reads = output, cut = cut)
}

anchored_front_trim_fastq <- function(file, settings, report_json,
                                      min.length = 20L, chunk_size = 100000L,
                                      fastp_path = install.fastp(), threads = 4L) {
  stopifnot(file.exists(file), file.exists(report_json), min.length >= 1L,
            chunk_size >= 1L, threads >= 1L)
  # The canonical FASTQ/report are published only after all streaming/stats work.
  folder <- dirname(file)
  intermediate <- tempfile('anchored_input_', tmpdir = folder, fileext = '.fastq')
  final <- tempfile('anchored_output_', tmpdir = folder, fileext = '.fastq')
  combined <- tempfile('anchored_report_', tmpdir = folder, fileext = '.json')
  on.exit(unlink(c(intermediate, final, combined)), add = TRUE)
  file.create(intermediate)
  streamer <- ShortRead::FastqStreamer(file, n = as.integer(chunk_size))
  on.exit(close(streamer), add = TRUE)
  histogram <- numeric(settings$max + 1L); records <- 0
  repeat {
    fq <- ShortRead::yield(streamer)
    if (!length(fq)) break
    z <- anchored_front_chunk(fq, settings$base, settings$min, settings$max)
    histogram <- histogram + tabulate(z$cut + 1L, nbins = length(histogram))
    records <- records + length(fq)
    ShortRead::writeFastq(z$reads, intermediate, mode = 'a', compress = FALSE)
  }
  if (!records) stop('No input reads available for anchored front cleanup')
  metrics <- file.path(dirname(report_json), 'trim_order_metrics')
  dir.create(metrics, showWarnings = FALSE, recursive = TRUE)
  id <- sub('^report_|\\.json$', '', basename(report_json))
  id <- sub('\\.json$', '', id)
  stats <- file.path(metrics, paste0(id, '_anchored_front.json'))
  html <- file.path(dirname(report_json), sub('\\.json$', '.html', basename(report_json)))
  args <- c('--in1', intermediate, '--out1', final, '--disable_adapter_trimming',
            '--disable_quality_filtering', '--disable_trim_poly_g', '--length_required', as.character(min.length),
            '--thread', as.character(threads), '--json', stats, '--html', html)
  status <- system2(path.expand(fastp_path), shQuote(args))
  if (status != 0L || !file.exists(final) || !file.exists(stats))
    stop('Anchored front final FASTQ/statistics pass failed')
  before <- jsonlite::fromJSON(report_json, simplifyVector = FALSE)
  after <- jsonlite::fromJSON(stats, simplifyVector = FALSE)
  if (as.numeric(after$summary$before_filtering$total_reads) != records)
    stop('Anchored front FASTQ record count does not match final statistics')
  # Preserve original-input -> final-insert flow, not just the last-stage flow.
  after$summary$before_filtering <- before$summary$before_filtering
  after$summary$sequencing <- before$summary$sequencing
  after$read1_before_filtering <- before$read1_before_filtering
  after$adapter_cutting <- before$adapter_cutting
  for (k in setdiff(union(names(before$filtering_result), names(after$filtering_result)), 'passed_filter_reads')) {
    a <- before$filtering_result[[k]]; b <- after$filtering_result[[k]]
    if (is.numeric(a) || is.numeric(b)) after$filtering_result[[k]] <- sum(c(a,b))
  }
  after$command <- paste(before$command, '[anchored_front_cleanup]',
                        paste(c(path.expand(fastp_path), args), collapse = ' '))
  after$massiveNGSpipe_trim <- before$massiveNGSpipe_trim
  after$massiveNGSpipe_anchored_front <- list(
    schema = 1L, base = settings$base, min = settings$min, max = settings$max,
    position = 'after_fixed_front_and_adapter_tail', input_records = records,
    matched_records = sum(histogram[-1L]), removed_bases = sum((0:settings$max)*histogram),
    cut_length_counts = as.list(stats::setNames(histogram, 0:settings$max)),
    minimum_final_length = min.length, stage_json = stats)
  jsonlite::write_json(after, combined, auto_unbox = TRUE, pretty = TRUE)
  if (!file.copy(report_json, file.path(metrics, paste0(id, '_before_anchored_front.json')), overwrite = TRUE))
    stop('Cannot preserve pre-anchored trimming report')
  if (!file.rename(final, file)) stop('Cannot publish anchored front FASTQ')
  if (!file.rename(combined, report_json)) stop('Cannot publish anchored front report; trim stage must rerun')
  invisible(after$massiveNGSpipe_anchored_front)
}
