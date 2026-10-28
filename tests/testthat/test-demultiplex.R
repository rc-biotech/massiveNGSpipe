library(data.table)

# Synthetic pool: [UMI5 2][insert][UMI3 5][barcode 14][adapter][tail], read length 150
make_pool <- function(dir, file, barcodes, n = 4000, seed = 1) {
  set.seed(seed)
  rnd <- function(k) paste(sample(c("A", "C", "G", "T"), k, replace = TRUE), collapse = "")
  adapter <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  truth <- list(); lines <- character(4 * n)
  for (i in seq_len(n)) {
    s <- sample(names(barcodes), 1); ins <- rnd(sample(10:40, 1)); bc <- barcodes[[s]]
    if (i %% 10 == 0) { p <- sample(14, 1); substr(bc, p, p) <- setdiff(c("A", "C", "G", "T"), substr(bc, p, p))[1] }
    seq <- substr(paste0(rnd(2), ins, rnd(5), bc, adapter, rnd(150)), 1, 150)
    if (i %% 25 == 0) seq <- rnd(150)                      # no adapter
    else if (nchar(ins) >= 15) truth[[s]] <- c(truth[[s]], ins)
    lines[(4 * i - 3):(4 * i)] <- c(paste0("@r", i), seq, "+", strrep("I", 150))
  }
  con <- gzfile(file.path(dir, file), "w"); writeLines(lines, con); close(con)
  truth
}

barcodes <- c(s01 = "CATATTCCTGGTGG", s02 = "ACAACAAGGTTCGG", s03 = "TCGGTGAGGATTGG", s04 = "GACAACATCAGCGG")

test_that("demultiplex assigns, trims UMIs, merges pool files and collapses exactly", {
  skip_if(Sys.which("cc") == "", "no C compiler")
  root <- file.path(tempdir(), "raw_data", "poolA"); pool <- file.path(root, "multiplexed")
  dir.create(pool, recursive = TRUE, showWarnings = FALSE)
  t1 <- make_pool(pool, "run1_1.fastq.gz", barcodes, seed = 1)
  t2 <- make_pool(pool, "run2_1.fastq.gz", barcodes, seed = 2)
  ctrl <- rbind(data.table(sample = paste0(names(barcodes), "_r1"), BARCODE = barcodes, file = "run1_1.fastq.gz",
                           final_sample = names(barcodes)),
                data.table(sample = paste0(names(barcodes), "_r2"), BARCODE = barcodes, file = "run2_1.fastq.gz",
                           final_sample = names(barcodes)))
  ctrl[, `:=`(umi_5p = 2, umi_3p = 5)]
  fwrite(ctrl, file.path(pool, "controlfile.csv"))
  out <- file.path(tempdir(), "processed_data", "poolA", "trim")
  st <- demultiplex(pool, out_dir = out, tmp_dir = file.path(tempdir(), "demux_tmp"), parallel_files = 2, collapse_workers = 2)
  for (s in names(barcodes)) {
    f <- file.path(out, "SINGLE", paste0("collapsed_trimmed_", s, ".fasta.gz"))
    expect_true(file.exists(f))
    fa <- Biostrings::readDNAStringSet(f)
    counts <- as.integer(sub(".*_x", "", names(fa)))
    got <- unname(rep(as.character(fa), counts))
    expect_equal(sort(got), sort(c(t1[[s]], t2[[s]])))
    expect_true(all(diff(counts) <= 0))                  # sorted by count
    js <- jsonlite::read_json(file.path(out, paste0("report_", s, ".json")))
    expect_equal(js$summary$after_filtering$total_reads, length(c(t1[[s]], t2[[s]])))
    expect_gte(js$summary$before_filtering$total_reads, js$summary$after_filtering$total_reads)
  }
  expect_true(file.exists(file.path(pool, "demultiplex_stats.csv")))
  ps <- fread(file.path(pool, "demultiplex_pool_stats.csv"))
  expect_equal(sum(ps[key == "total"]$reads), 8000)
  expect_true(all(st$reads_rescued > 0))                 # 1-mismatch barcodes rescued
})

test_that("demultiplex resumes completed pool files and validates the control file", {
  skip_if(Sys.which("cc") == "", "no C compiler")
  root <- file.path(tempdir(), "raw_data", "poolB"); pool <- file.path(root, "multiplexed")
  dir.create(pool, recursive = TRUE, showWarnings = FALSE)
  make_pool(pool, "L1_R1_001.fastq.gz", barcodes, n = 500)
  fwrite(data.table(sample = names(barcodes), BARCODE = barcodes, LANE = 1), file.path(pool, "controlfile.csv"))
  tmp <- file.path(tempdir(), "demux_tmp_B"); out <- file.path(tempdir(), "processed_data", "poolB", "trim")
  demultiplex(pool, out_dir = out, tmp_dir = tmp, keep_tmp = TRUE)            # legacy LANE -> L1_R1_001.fastq.gz
  expect_true(file.exists(file.path(tmp, "stats.L1_R1_001.tsv")))
  t0 <- file.mtime(file.path(tmp, "stats.L1_R1_001.tsv")); Sys.sleep(1.1)
  expect_message(demultiplex(pool, out_dir = out, tmp_dir = tmp, keep_tmp = TRUE), "reused L1_R1_001")   # resumed
  expect_equal(file.mtime(file.path(tmp, "stats.L1_R1_001.tsv")), t0)
  fwrite(data.table(sample = c("a", "b"), BARCODE = c("ACGT", "ACGT"), file = "x.fastq.gz"), file.path(pool, "bad.csv"))
  expect_error(demultiplex(pool, controlfile = file.path(pool, "bad.csv")), "not found|Duplicated")
})

test_that("demultiplex stage is off by default and first when enabled", {
  expect_false("demultiplexed" %in% libtype_flags("Ribo-seq", "local"))
  f <- libtype_flags("Ribo-seq", "local", demultiplex = TRUE)
  expect_equal(unname(f[1]), "demultiplexed")
  expect_equal(names(f)[1], "pipe_demultiplex")
  expect_equal(unname(libtype_flags("Ribo-seq", "online", demultiplex = TRUE)[1:3]), c("start", "fetch", "demultiplexed"))
})

test_that("pipe_demultiplex delivers pooled runs to their experiments and sets markers and flags", {
  skip_if(Sys.which("cc") == "", "no C compiler")
  base <- file.path(tempdir(), "mngsp_demux_stage"); unlink(base, recursive = TRUE)
  raw <- file.path(base, "raw_data"); bam <- file.path(base, "processed_data"); proj <- file.path(base, "local_pipeline_test")
  pool <- file.path(raw, "poolC", "multiplexed"); dir.create(pool, recursive = TRUE)
  dir.create(bam, recursive = TRUE); dir.create(proj, recursive = TRUE)
  truth <- make_pool(pool, "lane1_1.fastq.gz", barcodes, n = 1000)
  fwrite(data.table(sample = names(barcodes), BARCODE = barcodes, file = "lane1_1.fastq.gz", umi_5p = 2, umi_3p = 5),
         file.path(pool, "controlfile.csv"))
  conf_orfik <- c(fastq = raw, bam = bam, ref = file.path(base, "ref"), exp = file.path(base, "exp"))
  config <- pipeline_config(project_dir = proj, config = conf_orfik, mode = "local", demultiplex = TRUE, google_url = NULL)
  config$error_dir <- file.path(base, "error_logs")
  expect_true("demultiplexed" %in% names(config$flag))
  # two experiments from one pool: s01/s02 yeast study, s03 a second organism; s04 not curated
  mk <- function(acc, runs, org, exp) list(accession = acc, study = data.table(Run = runs, ScientificName = org),
                                          organisms = setNames(list(list(conf = c(exp = exp, bam = file.path(bam, exp), fastq = file.path(raw, acc)))), org))
  pipelines <- list(studyA = mk("studyA", c("s01", "s02"), "Saccharomyces cerevisiae", "studyA-saccharomyces_cerevisiae"),
                    studyB = mk("studyB", "s03", "Homo sapiens", "studyB-homo_sapiens"))
  pipe_demultiplex(pipelines, config)
  fa <- file.path(bam, "studyA-saccharomyces_cerevisiae", "trim", "SINGLE", "collapsed_trimmed_s02.fasta.gz")
  expect_true(file.exists(fa))
  expect_true(file.exists(file.path(bam, "studyB-homo_sapiens", "trim", "SINGLE", "collapsed_trimmed_s03.fasta.gz")))
  expect_false(file.exists(file.path(bam, "studyA-saccharomyces_cerevisiae", "trim", "SINGLE", "collapsed_trimmed_s03.fasta.gz")))
  expect_setequal(samples_done(config, "collapsed", "studyA-saccharomyces_cerevisiae"), c("s01", "s02"))
  expect_setequal(samples_done(config, "trim", "studyB-homo_sapiens"), "s03")
  expect_true(step_is_done(config, "demultiplexed", "studyA-saccharomyces_cerevisiae"))
  expect_true(step_is_done(config, "demultiplexed", "studyB-homo_sapiens"))
  n_expected <- length(truth[["s02"]]); got <- Biostrings::readDNAStringSet(fa)
  expect_equal(sum(as.integer(sub(".*_x", "", names(got)))), n_expected)
  t0 <- file.mtime(fa); Sys.sleep(1.1)
  pipe_demultiplex(pipelines, config)                                       # second call: nothing to do
  expect_equal(file.mtime(fa), t0)
})
