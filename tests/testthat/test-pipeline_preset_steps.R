# Unit tests for the pure helpers in R/pipeline_preset_steps.R
# (libtype_long_to_short/short_to_long, library_strategies_allowed,
# safe_se_cbind). The pipe_*()/pipeline_collection_*()/pipeline_merge_*()
# orchestrators are out of scope -- they need real ORFik experiments,
# BAMs, and STAR/BiocParallel infra.

test_that("libtype_long_to_short / libtype_short_to_long recode the known pairs and pass through unknowns", {
  expect_identical(libtype_long_to_short("Ribo-seq"), "RFP")
  expect_identical(libtype_long_to_short("RNA-seq"), "RNA")
  expect_identical(libtype_long_to_short("disome"), "disome") # passthrough
  expect_identical(libtype_long_to_short(c("Ribo-seq", "RNA-seq")), c("RFP", "RNA"))

  expect_identical(libtype_short_to_long("RFP"), "Ribo-seq")
  expect_identical(libtype_short_to_long("RNA"), "RNA-seq")
  expect_identical(libtype_short_to_long("disome"), "disome")
})

test_that("libtype_long_to_short / libtype_short_to_long round-trip for the two known pairs", {
  expect_identical(libtype_short_to_long(libtype_long_to_short("Ribo-seq")), "Ribo-seq")
  expect_identical(libtype_short_to_long(libtype_long_to_short("RNA-seq")), "RNA-seq")
})

test_that("library_strategies_allowed returns the expected fixed set", {
  expect_setequal(library_strategies_allowed(),
                  c("RNA-Seq", "Ribo-seq", "ssRNA-seq", "ncRNA-Seq", "miRNA-Seq", "OTHER"))
})

test_that("safe_se_cbind cbinds matching-nrow SummarizedExperiments on the happy path", {
  se1 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2)))
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(5:8, nrow = 2)))
  result <- safe_se_cbind(list(a = se1, b = se2))
  expect_identical(ncol(result), 4L)
})

test_that("safe_se_cbind reports a clear diagnostic error for mismatched row counts", {
  se_2row <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2)))
  se_3row <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:6, nrow = 3)))
  err <- tryCatch(
    safe_se_cbind(list(good1 = se_2row, good2 = se_2row, bad = se_3row)),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "bad", fixed = TRUE)
  expect_match(conditionMessage(err), "nrow=3", fixed = TRUE)
})

test_that("pipeline_cleanup() calls save_expanded_alignment_metrics() as its final step when this preset collapses reads", {
  config <- fake_config(preset = "Ribo-seq", mode = "online")
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(organism = "Homo sapiens", bam_dir = bam_dir,
                              runs = data.table::data.table(
                                Run = "SRR001", LibraryLayout = "SINGLE",
                                LIBRARYTYPE = "RFP", ScientificName = "Homo sapiens"))
  pipeline <- pipelines[[1]]
  exp <- pipeline$organisms[["Homo sapiens"]]$conf["exp"]

  fake_mark_all_done(config, c("start", "fetch", "trim", "collapsed", "aligned"), exp)

  aligned_dir <- file.path(bam_dir, "aligned")
  dir.create(aligned_dir, recursive = TRUE)
  file.create(file.path(aligned_dir, "SRR001_Aligned.sortedByCoord.out.bam"))

  seen_args <- NULL
  testthat::local_mocked_bindings(
    save_expanded_alignment_metrics = function(bam_dir, collapsed_dir, study_org, ...) {
      seen_args <<- list(bam_dir = bam_dir, collapsed_dir = collapsed_dir, study_org = study_org)
    }
  )

  pipeline_cleanup(pipeline, config)

  expect_true(step_is_done(config, "cleanbam", exp))
  expect_true(file.exists(file.path(aligned_dir, "SRR001.bam")))
  expect_false(is.null(seen_args))
  expect_identical(as.character(seen_args$bam_dir), as.character(aligned_dir))
  expect_identical(as.character(seen_args$collapsed_dir), file.path(bam_dir, "trim", "SINGLE"))
  expect_identical(seen_args$study_org$Run, "SRR001")
})

test_that("pipeline_cleanup() does NOT call save_expanded_alignment_metrics() for a preset without collapsing", {
  config <- fake_config(preset = "RNA-seq", mode = "online") # no "collapsed" flag for RNA-seq
  bam_dir <- tempfile("bam_")
  pipelines <- fake_pipelines(organism = "Homo sapiens", bam_dir = bam_dir,
                              runs = data.table::data.table(
                                Run = "SRR001", LibraryLayout = "SINGLE",
                                LIBRARYTYPE = "RNA", ScientificName = "Homo sapiens"))
  pipeline <- pipelines[[1]]
  exp <- pipeline$organisms[["Homo sapiens"]]$conf["exp"]

  fake_mark_all_done(config, c("start", "fetch", "aligned"), exp)

  aligned_dir <- file.path(bam_dir, "aligned")
  dir.create(aligned_dir, recursive = TRUE)
  file.create(file.path(aligned_dir, "SRR001_Aligned.sortedByCoord.out.bam"))

  called <- FALSE
  testthat::local_mocked_bindings(
    save_expanded_alignment_metrics = function(...) called <<- TRUE
  )

  pipeline_cleanup(pipeline, config)

  expect_true(step_is_done(config, "cleanbam", exp))
  expect_false(called)
})
