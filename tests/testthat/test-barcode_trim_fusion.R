# Adapter/front cuts may share a pass, but fastp fixed tail cuts occur before
# adapter matching. Tail barcodes require the ordered streaming backend.
# Biological sequence preservation is covered by the real fastp regression.

test_that("tail barcodes use post-adapter trimming rather than the fused raw cut", {
  seen <- NULL
  testthat::local_mocked_bindings(
    barcode_detector_single = function(...) data.table::data.table(adapter="AGATCGGAAGAGC",barcode_detected=TRUE,barcode5p_size=4,barcode3p_size=6),
    fastp_post_adapter_tail = function(file, file2, output.dir, adapter, trim.front, trim.tail) {
      seen <<- list(file=file,file2=file2,adapter=adapter,front=trim.front,tail=trim.tail)
    }
  )
  testthat::local_mocked_bindings("STAR.align.single"=function(...)stop("Unsafe fused trim called"),.package="ORFik")
  result<-run_barcode_detection_and_trim(data.table::data.table(Run="SRR001",LibraryLayout="SINGLE"),"fastq_dir","target_dir","trimmed_dir",mode="local",file="raw.fastq")
  expect_equal(seen,list(file="raw.fastq",file2=NULL,adapter="AGATCGGAAGAGC",front=4,tail=6))
  expect_true(result$barcode_detected)
})

test_that("run_barcode_detection_and_trim() still does exactly one pass when NO barcode is detected", {
  # This is the key behavior change: previously, "no barcode" meant the
  # earlier (separate) adapter-only pass's output was just left as the
  # final result, with no second pass. Now there is no earlier pass at
  # all -- the single real pass here (with 0/0 sizes, a harmless no-op)
  # IS the adapter trim, so it must still happen exactly once.
  calls <- list()
  testthat::local_mocked_bindings(
    barcode_detector_single = function(...) {
      data.table::data.table(adapter = "AGATCGGAAGAGC", barcode_detected = FALSE,
                             barcode5p_size = 0, barcode3p_size = 0)
    }
  )
  testthat::local_mocked_bindings(
    "STAR.align.single" = function(file1, file2 = NULL, output.dir, adapter.sequence,
                                   index.dir, steps, trim.front = 0, trim.tail = 0, ...) {
      calls[[length(calls) + 1]] <<- list(trim.front = trim.front, trim.tail = trim.tail)
      invisible(NULL)
    },
    .package = "ORFik"
  )

  study_sample <- data.table::data.table(Run = "SRR002", LibraryLayout = "SINGLE")
  result <- run_barcode_detection_and_trim(study_sample, "fastq_dir", "target_dir", "trimmed_dir",
                                           mode = "local", file = "raw.fastq")

  expect_length(calls, 1)
  expect_equal(calls[[1]]$trim.front, 0)
  expect_equal(calls[[1]]$trim.tail, 0)
  expect_false(result$barcode_detected)
})

test_that("run_barcode_detection_and_trim() passes detect_reads_to_process through to detection", {
  seen_reads_to_process <- NULL
  testthat::local_mocked_bindings(
    barcode_detector_single = function(study_sample, fastq_dir, process_dir, trimmed_dir,
                                       redownload_raw_if_needed, detect_reads_to_process) {
      seen_reads_to_process <<- detect_reads_to_process
      data.table::data.table(adapter = "disable", barcode_detected = FALSE,
                             barcode5p_size = 0, barcode3p_size = 0)
    }
  )
  testthat::local_mocked_bindings("STAR.align.single" = function(...) invisible(NULL), .package = "ORFik")

  study_sample <- data.table::data.table(Run = "SRR003", LibraryLayout = "SINGLE")
  run_barcode_detection_and_trim(study_sample, "fastq_dir", "target_dir", "trimmed_dir",
                                 mode = "online", file = "raw.fastq",
                                 detect_reads_to_process = 250000)

  expect_equal(seen_reads_to_process, 250000)
})

test_that("run_barcode_detection_and_trim() never passes the 'passed' sentinel as adapter.sequence", {
  # Regression test for a real bug found via live end-to-end testing:
  # barcode_dt$adapter can be the literal string "passed" (fastp's own
  # "nothing was cut" JSON sentinel -- see barcode_detector_single()),
  # which STAR.align.single()'s shell script rejects outright ("adapter
  # can only have bases in {A, T, C, G}"). Must be sanitized to "auto".
  seen_adapter <- NULL
  testthat::local_mocked_bindings(
    barcode_detector_single = function(...) {
      data.table::data.table(adapter = "passed", barcode_detected = FALSE,
                             barcode5p_size = 0, barcode3p_size = 0)
    }
  )
  testthat::local_mocked_bindings(
    "STAR.align.single" = function(file1, file2 = NULL, output.dir, adapter.sequence, ...) {
      seen_adapter <<- adapter.sequence
      invisible(NULL)
    },
    .package = "ORFik"
  )

  study_sample <- data.table::data.table(Run = "SRR005", LibraryLayout = "SINGLE")
  run_barcode_detection_and_trim(study_sample, "fastq_dir", "target_dir", "trimmed_dir",
                                 mode = "local", file = "raw.fastq")

  expect_identical(seen_adapter, "auto")
})

test_that("run_barcode_detection_and_trim() passes a real detected adapter through unchanged", {
  seen_adapter <- NULL
  testthat::local_mocked_bindings(
    barcode_detector_single = function(...) {
      data.table::data.table(adapter = "AGATCGGAAGAGC", barcode_detected = TRUE,
                             barcode5p_size = 4, barcode3p_size = 0)
    }
  )
  testthat::local_mocked_bindings(
    "STAR.align.single" = function(file1, file2 = NULL, output.dir, adapter.sequence, ...) {
      seen_adapter <<- adapter.sequence
      invisible(NULL)
    },
    .package = "ORFik"
  )

  study_sample <- data.table::data.table(Run = "SRR006", LibraryLayout = "SINGLE")
  run_barcode_detection_and_trim(study_sample, "fastq_dir", "target_dir", "trimmed_dir",
                                 mode = "local", file = "raw.fastq")

  expect_identical(seen_adapter, "AGATCGGAAGAGC")
})

test_that("run_barcode_detection_and_trim() passes file2 through for paired-end input", {
  seen_file2 <- NULL
  testthat::local_mocked_bindings(
    barcode_detector_single = function(...) {
      data.table::data.table(adapter = "disable", barcode_detected = FALSE,
                             barcode5p_size = 0, barcode3p_size = 0)
    }
  )
  testthat::local_mocked_bindings(
    "STAR.align.single" = function(file1, file2 = NULL, ...) {
      seen_file2 <<- file2
      invisible(NULL)
    },
    .package = "ORFik"
  )

  study_sample <- data.table::data.table(Run = "SRR004", LibraryLayout = "PAIRED")
  run_barcode_detection_and_trim(study_sample, "fastq_dir", "target_dir", "trimmed_dir",
                                 mode = "local", file = "raw_1.fastq", file2 = "raw_2.fastq")

  expect_identical(seen_file2, "raw_2.fastq")
})


test_that("an explicit adapter override survives an empty fastp adapter self-report", {
  folder<-tempfile("adapter_override_");dir.create(folder)
  data.table::fwrite(data.table::data.table(Run="SRR010",adapter="AGATCGGAAGAG"),file.path(folder,"adapters_manual.csv"))
  seen<-NULL
  testthat::local_mocked_bindings(
    barcode_detector_single=function(...)data.table::data.table(adapter="passed",barcode_detected=TRUE,barcode5p_size=2,barcode3p_size=10),
    fastp_post_adapter_tail=function(file,file2,output.dir,adapter,...)seen<<-adapter
  )
  run_barcode_detection_and_trim(data.table::data.table(Run="SRR010",LibraryLayout="SINGLE"),"raw","output",folder,mode="local",file="raw.fastq")
  expect_identical(seen,"AGATCGGAAGAG")
})
