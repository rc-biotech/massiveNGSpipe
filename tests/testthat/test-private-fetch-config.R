test_that("failed extraction preserves source archive on both scratch paths", {
  for (scratch in c(FALSE, TRUE)) {
    env <- new.env(parent = globalenv())
    sys.source(testthat::test_path('../../R/download_sra.R'), env)
    d <- tempfile(); dir.create(d)
    archive <- file.path(d, 'SRR001'); writeLines('cached archive', archive)
    original <- readBin(archive, 'raw', file.info(archive)$size)
    binary <- tempfile(); dir.create(binary)
    file.create(file.path(binary, 'fasterq-dump'))
    work <- tempfile()
    env$tmp_claim_root <- function() tempdir()
    env$fasterq_tmp_need_gb <- function(...) 1
    env$tmp_space_claim <- function(...) if (scratch) 'claim' else NULL
    env$tmp_space_release <- function(...) unlink(work, recursive = TRUE)
    env$system2 <- function(command, args, ...) {
      input <- args[1]
      writeLines('partial', paste0(input, '_1.fastq'))
      3L
    }
    # Supply a deterministic work directory through tempdir() in this scope.
    env$tempdir <- function() dirname(work)
    work <- file.path(dirname(work), 'fasterq_SRR001')
    expect_error(env$sra_to_fastq('SRR001', d, FALSE, binary, TRUE),
                 'Fasterq dump failed')
    expect_true(file.exists(archive))
    expect_identical(readBin(archive, 'raw', file.info(archive)$size), original)
    expect_false(file.exists(paste0(archive, '_1.fastq')))
    unlink(c(d, binary, work), recursive = TRUE)
  }
})

test_that("pipeline passes only explicitly configured private toolkit", {
  env <- new.env(parent = globalenv())
  sys.source(testthat::test_path('../../R/pipeline_preset_steps_sub.R'), env)
  env$step_is_done <- function(...) FALSE
  env$set_flag <- env$set_sample_flag <- function(...) invisible(NULL)
  env$samples_done <- function(...) character()
  env$pipeline_log_base <- function(...) tempdir()
  env$run_experiment_subprocess <- function(func, args, ...) do.call(func, args)
  calls <- list()
  env$download_sra <- function(...) { calls[[length(calls) + 1L]] <<- list(...) }
  info <- data.table::data.table(Run = 'SRR001', ScientificName = 'human', LibraryLayout = 'PAIRED')
  pipeline <- list(study = info, organisms = list(human = list(conf = c(exp = 'study', fastq = tempdir()))))
  env$pipeline_download(pipeline, list(compress_raw_data = FALSE))
  env$pipeline_download(pipeline, list(compress_raw_data = FALSE, sratoolkit_path_extra = 'wrong'))
  env$pipeline_download(pipeline, list(compress_raw_data = FALSE, sratoolkit_path = '/private/bin'))
  expect_false('sratoolkit_path' %in% names(calls[[1]]))
  expect_false('sratoolkit_path' %in% names(calls[[2]]))
  expect_identical(calls[[3]]$sratoolkit_path, '/private/bin')
})

test_that("both fallback attempts run serially inside the fetch subprocess", {
  seen <- list()
  testthat::local_mocked_bindings(download.SRA = function(...) {
    args <- list(...); seen[[length(seen) + 1L]] <<- args
    if (length(seen) == 1L) stop('simulated EBI failure')
    character()
  }, .package = 'ORFik')
  expect_output(download_sra_ebi('SRR001', tempdir(), FALSE), NA)
  expect_length(seen, 2L)
  expect_true(all(vapply(seen, function(x) inherits(x$BPPARAM, 'SerialParam'), logical(1))))
  expect_false(seen[[2]]$use.ebi.ftp)
})
