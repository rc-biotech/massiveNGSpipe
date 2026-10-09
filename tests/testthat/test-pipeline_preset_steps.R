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

# pipe_align_clean()'s flattened, organism-sorted dispatch -- see
# R/star_index_lock.R and this function's own header comment for why:
# grouping every study of the same organism contiguously means an
# organism's STAR shared-memory index is loaded from disk once per
# organism per run, not once per study, and lets one organism's failure
# in a mixed-organism study stop blocking a sibling organism that
# already succeeded.

#' Build a two-study, mixed-organism `pipelines` list matching the
#' exact scenario this feature was designed for: study1 spans
#' organisms A and B, study2 spans organisms A and C.
fake_mixed_organism_pipelines <- function() {
  mk_pipeline <- function(accession, organisms) {
    list(accession = accession,
        organisms = stats::setNames(
          lapply(organisms, function(o) list(conf = c(exp = paste0(accession, "-", o), bam = tempfile()),
                                             index = tempfile())),
          organisms),
        study = data.table::data.table(Run = "SRR1", ScientificName = organisms[1]))
  }
  stats::setNames(
    list(mk_pipeline("study1", c("A", "B")), mk_pipeline("study2", c("A", "C"))),
    c("study1", "study2"))
}

test_that("pipe_align_clean() groups every study of the same organism contiguously", {
  config <- fake_config(extra = list(error_dir = tempfile("err_")))
  pipelines <- fake_mixed_organism_pipelines()

  seen <- character()
  testthat::local_mocked_bindings(
    pipeline_align_contaminants_one_organism = function(...) invisible(NULL),
    pipeline_align_one_organism = function(pipeline, organism, config, pipelines, keep_loaded_after) {
      seen <<- c(seen, paste0(pipeline$accession, "-", organism))
      invisible(NULL)
    },
    pipeline_cleanup_one_organism = function(...) invisible(NULL)
  )

  pipe_align_clean(pipelines, config)

  expect_length(seen, 4)
  # Every organism's units (however many) must be contiguous -- no
  # other organism's unit interleaved between them.
  organism_of <- sub("^.*-", "", seen)
  for (org in unique(organism_of)) {
    positions <- which(organism_of == org)
    expect_identical(positions, seq(positions[1], positions[1] + length(positions) - 1),
                     info = paste("organism", org, "not contiguous in", paste(seen, collapse = ", ")))
  }
  # Both studies' organism A units must still both be present.
  expect_true("study1-A" %in% seen)
  expect_true("study2-A" %in% seen)
})

test_that("pipe_align_clean() passes keep_loaded_after = TRUE only when the NEXT unit shares this organism", {
  config <- fake_config(extra = list(error_dir = tempfile("err_")))
  pipelines <- fake_mixed_organism_pipelines()

  keep_loaded_seen <- list()
  testthat::local_mocked_bindings(
    pipeline_align_contaminants_one_organism = function(...) invisible(NULL),
    pipeline_align_one_organism = function(pipeline, organism, config, pipelines, keep_loaded_after) {
      keep_loaded_seen[[paste0(pipeline$accession, "-", organism)]] <<- keep_loaded_after
      invisible(NULL)
    },
    pipeline_cleanup_one_organism = function(...) invisible(NULL)
  )

  pipe_align_clean(pipelines, config)

  # Organism A has two units (study1, study2); whichever runs FIRST
  # among them must keep the index loaded for the one right after it.
  # Organisms B and C have exactly one unit each -- never kept loaded.
  expect_false(keep_loaded_seen[["study1-B"]])
  expect_false(keep_loaded_seen[["study2-C"]])
  a_keep_flags <- c(keep_loaded_seen[["study1-A"]], keep_loaded_seen[["study2-A"]])
  expect_equal(sum(a_keep_flags), 1) # exactly one of the two A-units keeps it loaded for the other
})

test_that("pipe_align_clean() isolates one organism's align failure from a sibling organism of the SAME study", {
  # error_dir defaults to NULL in fake_config() (report_failed_pipe_path()
  # then returns a fresh random tempfile() every call, never actually
  # written) -- a real path is needed here since this test checks
  # file.exists() against it.
  config <- fake_config(extra = list(error_dir = tempfile("err_")))
  pipelines <- fake_mixed_organism_pipelines()

  cleaned <- character()
  testthat::local_mocked_bindings(
    pipeline_align_contaminants_one_organism = function(...) invisible(NULL),
    pipeline_align_one_organism = function(pipeline, organism, config, pipelines, keep_loaded_after) {
      if (pipeline$accession == "study1" && organism == "B") stop("simulated STAR crash for study1-B")
      invisible(NULL)
    },
    pipeline_cleanup_one_organism = function(pipeline, organism, config) {
      cleaned <<- c(cleaned, paste0(pipeline$accession, "-", organism))
      invisible(NULL)
    }
  )

  suppressWarnings(pipe_align_clean(pipelines, config))

  # study1-A, study2-A, study2-C all succeed and reach cleanup;
  # study1-B fails and must NOT reach cleanup, but must also NOT have
  # blocked any of the other three organism units.
  expect_setequal(cleaned, c("study1-A", "study2-A", "study2-C"))
  expect_true(file.exists(report_failed_pipe_path(config, c(exp = "study1-B"))))
  expect_false(file.exists(report_failed_pipe_path(config, c(exp = "study1-A"))))
})
