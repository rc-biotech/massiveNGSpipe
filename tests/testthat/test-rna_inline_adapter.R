test_that("RNA adapters are cached by Run and actual input, including negative decisions", {
  d <- tempfile(); dir.create(d)
  a <- file.path(d, "long.fastq"); b <- file.path(d, "short.fastq")
  writeLines("long", a); writeLines("short", b)
  called <- character()
  testthat::local_mocked_bindings(short_insert_adapter = function(file) {
    called <<- c(called, basename(file))
    if (basename(file) == "long.fastq") NA_character_ else "AGATCGGAAGAGC"
  })
  expect_true(is.na(rna_inline_adapter(a, "SRR1", d)))
  expect_identical(rna_inline_adapter(b, "SRR2", d), "AGATCGGAAGAGC")
  expect_true(is.na(rna_inline_adapter(a, "SRR1", d)))
  expect_identical(rna_inline_adapter(b, "SRR2", d), "AGATCGGAAGAGC")
  expect_identical(called, c("long.fastq", "short.fastq"))
  # Changing the actual selected mate/input invalidates the same Run's cache.
  expect_identical(rna_inline_adapter(b, "SRR1", d), "AGATCGGAAGAGC")
  writeLines("short changed", b)
  rna_inline_adapter(b, "SRR2", d)
  expect_length(called, 4L)
})

test_that("raw RNA retention is independent of delete_collapsed_files", {
  d <- tempfile(); dir.create(file.path(d, "aligned"), recursive = TRUE)
  raw <- file.path(d, "SRR1.fastq"); writeLines("raw", raw)
  writeLines("bam", file.path(d, "aligned", "SRR1.bam"))
  config <- fake_config(preset = "RNA-seq", extra = list(
    delete_raw_files = FALSE, delete_collapsed_files = TRUE))
  testthat::local_mocked_bindings(system2 = function(...) 0L, .package = "base")
  testthat::local_mocked_bindings(STAR.allsteps.multiQC = function(...) NULL)
  testthat::local_mocked_bindings(check_alignment_rate = function(...) NULL)
  alignment_final_checks(d, d, data.table::data.table(Run = "SRR1"), list(raw), config, "ge")
  expect_true(file.exists(raw))
  config$delete_raw_files <- TRUE
  alignment_final_checks(d, d, data.table::data.table(Run = "SRR1"), list(raw), config, "ge")
  expect_false(file.exists(raw))
})

test_that("one mixed RNA study aligns long and short inputs with their own adapter and steps", {
  d <- tempfile(); dir.create(d)
  raw <- file.path(d, "raw"); dir.create(raw)
  files <- file.path(raw, c("SRR1.fastq", "SRR2.fastq"))
  for (f in files) writeLines(c("@r", "ACGTACGT", "+", "IIIIIIII"), f)
  config <- fake_config(preset = "RNA-seq", mode = "online")
  p <- fake_pipelines(bam_dir = file.path(d, "processed"), runs = data.table::data.table(
    Run = c("SRR1", "SRR2"), LibraryLayout = "SINGLE", LIBRARYTYPE = "RNA",
    ScientificName = "Homo sapiens"))[[1]]
  p$organisms[[1]]$conf <- c(p$organisms[[1]]$conf, fastq = raw)
  fake_mark_all_done(config, c("start", "fetch"), p$organisms[[1]]$conf['exp'])
  aligned <- list()
  testthat::local_mocked_bindings(
    STAR.install = function(...) "STAR", install.fastp = function(...) "fastp",
    STAR.align.single = function(file1, ..., steps, adapter.sequence, resume) {
      aligned[[basename(file1)]] <<- list(steps=steps, adapter=adapter.sequence, resume=resume)
    }, .package = "ORFik")
  testthat::local_mocked_bindings(
    short_insert_adapter = function(...) NA_character_,
    rna_inline_adapter = function(file, run, ...) if (run == "SRR1") NA_character_ else "AGATCGGA",
    run_experiment_subprocess = function(func, args, ...) do.call(func, args),
    estimate_bam_sort_ram = function(...) 1,
    bam_sort_ram_decision = function(...) list(decision="proceed"),
    resolve_align_scratch = function(file1, output_dir, ...) list(file1=file1, output.dir=output_dir, used_ssd=FALSE),
    alignment_final_checks = function(...) NULL)
  pipeline_align_one_organism(p, "Homo sapiens", config)
  expect_identical(aligned[['SRR1.fastq']], list(steps="ge", adapter="auto", resume="ge"))
  expect_identical(aligned[['SRR2.fastq']], list(steps="tr-ge", adapter="AGATCGGA", resume=NULL))
})

test_that("real source-tree subprocess resolves the FASTQ alphabet guard", {
  d <- tempfile(); dir.create(d)
  base <- file.path(d, "base.fastq"); color <- file.path(d, "color.fastq")
  writeLines(c("@r", "ACGTACGT", "+", "IIIIIIII"), base)
  writeLines(c("@r", "T0123012", "+", "IIIIIIII"), color)
  f <- function(path) { assert_base_space_fastq(path, "fixture"); TRUE }
  expect_true(run_experiment_subprocess(f, list(path = base),
    file.path(d, "base.out"), file.path(d, "base.err"), poll_interval = .1))
  expect_error(run_experiment_subprocess(f, list(path = color),
    file.path(d, "color.out"), file.path(d, "color.err"), poll_interval = .1), "colour-space")
})
