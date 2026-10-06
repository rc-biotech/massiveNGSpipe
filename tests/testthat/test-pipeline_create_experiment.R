# Unit tests for pipeline_create_experiment() (R/pipeline_preset_steps_sub.R)
# -- specifically the author-value sanitization fix. ORFik::create.experiment()'s
# own `if (author != "")` check errors ("missing value where TRUE/FALSE
# needed") on a literal NA author, not just an empty string. Confirmed
# live, PRJEB50305-saccharomyces_cerevisiae, 2026-10-06: AUTHOR was
# genuinely NA for every sample in that study's metadata, and
# unique(study$AUTHOR) passed that NA straight through to
# ORFik::create.experiment().

#' Minimal study data.table with every column cleanup_metadata_for_exp()
#' and pipeline_create_experiment() itself read
#' @noRd
fake_study_for_experiment <- function(author = NA_character_, run_ids = c("SRR001", "SRR002")) {
  n <- length(run_ids)
  data.table::data.table(
    Run = run_ids, ScientificName = "Homo sapiens", LibraryLayout = "SINGLE",
    LIBRARYTYPE = "RFP", REPLICATE = seq_len(n), CONDITION = "WT", GENE = NA_character_,
    FRACTION = "", TIMEPOINT = "", BATCH = "", INHIBITOR = "chx",
    CELL_LINE = "NONE", TISSUE = "NONE", AUTHOR = author
  )
}

test_that("pipeline_create_experiment() never passes a literal NA author to ORFik::create.experiment()", {
  config <- fake_config(preset = "Ribo-seq")
  prior_steps <- names(config$flag); prior_steps <- prior_steps[seq_len(which(prior_steps == "exp") - 1)]
  fake_mark_all_done(config, prior_steps, "PRJNA000001-homo_sapiens")
  bam_dir <- tempfile("bam_")
  pipeline <- list(
    accession = "PRJNA000001",
    study = fake_study_for_experiment(author = NA_character_),
    organisms = list(`Homo sapiens` = list(
      conf = c(exp = "PRJNA000001-homo_sapiens", bam = bam_dir),
      annotation = c(gtf = "fake.gtf", genome = "fake.fa")
    ))
  )
  captured_author <- NULL
  testthat::local_mocked_bindings(
    match_bam_to_metadata = function(...) c("SRR001.bam", "SRR002.bam")
  )
  testthat::local_mocked_bindings(
    create.experiment = function(..., author) { captured_author <<- author; invisible(NULL) },
    read.experiment = function(...) "fake_df",
    .package = "ORFik"
  )

  pipeline_create_experiment(pipeline, config)

  expect_false(is.na(captured_author))
  expect_identical(captured_author, "")
})

test_that("pipeline_create_experiment() passes through a real, non-NA author unchanged", {
  config <- fake_config(preset = "Ribo-seq")
  prior_steps <- names(config$flag); prior_steps <- prior_steps[seq_len(which(prior_steps == "exp") - 1)]
  fake_mark_all_done(config, prior_steps, "PRJNA000001-homo_sapiens")
  bam_dir <- tempfile("bam_")
  pipeline <- list(
    accession = "PRJNA000001",
    study = fake_study_for_experiment(author = "Smith"),
    organisms = list(`Homo sapiens` = list(
      conf = c(exp = "PRJNA000001-homo_sapiens", bam = bam_dir),
      annotation = c(gtf = "fake.gtf", genome = "fake.fa")
    ))
  )
  captured_author <- NULL
  testthat::local_mocked_bindings(
    match_bam_to_metadata = function(...) c("SRR001.bam", "SRR002.bam")
  )
  testthat::local_mocked_bindings(
    create.experiment = function(..., author) { captured_author <<- author; invisible(NULL) },
    read.experiment = function(...) "fake_df",
    .package = "ORFik"
  )

  pipeline_create_experiment(pipeline, config)

  expect_identical(captured_author, "Smith")
})

test_that("pipeline_create_experiment() collapses more than one distinct author into a single scalar", {
  config <- fake_config(preset = "Ribo-seq")
  prior_steps <- names(config$flag); prior_steps <- prior_steps[seq_len(which(prior_steps == "exp") - 1)]
  fake_mark_all_done(config, prior_steps, "PRJNA000001-homo_sapiens")
  bam_dir <- tempfile("bam_")
  study <- fake_study_for_experiment(author = c("Smith", "Jones"))
  pipeline <- list(
    accession = "PRJNA000001",
    study = study,
    organisms = list(`Homo sapiens` = list(
      conf = c(exp = "PRJNA000001-homo_sapiens", bam = bam_dir),
      annotation = c(gtf = "fake.gtf", genome = "fake.fa")
    ))
  )
  captured_author <- NULL
  testthat::local_mocked_bindings(
    match_bam_to_metadata = function(...) c("SRR001.bam", "SRR002.bam")
  )
  testthat::local_mocked_bindings(
    create.experiment = function(..., author) { captured_author <<- author; invisible(NULL) },
    read.experiment = function(...) "fake_df",
    .package = "ORFik"
  )

  pipeline_create_experiment(pipeline, config)

  expect_length(captured_author, 1)
  expect_identical(captured_author, "Smith; Jones")
})
