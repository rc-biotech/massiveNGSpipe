# Unit tests for R/pipeline_align_memory.R -- the dynamic BAM-sort RAM
# estimate/decision/defer-marker logic built from a real production
# failure (PRJNA599943-homo_sapiens_RNA-seq, SRR11294057, 2026-10-09).
# get_system_usage() (from ORFik, real system calls) is mocked throughout,
# matching the pattern already used in test-pipeline_collapse.R.

test_that("estimate_bam_sort_ram scales with file size and has a 30GB floor", {
  big <- tempfile(); writeLines("x", big)
  # file.size is tiny for a real temp file -- mock it to a realistic size.
  testthat::local_mocked_bindings(
    file.info = function(...) data.frame(size = 20e9, row.names = c(...)[1]), .package = "base"
  )
  expect_equal(estimate_bam_sort_ram(big), 80e9) # 20GB * 4

  small <- tempfile(); writeLines("x", small)
  testthat::local_mocked_bindings(
    file.info = function(...) data.frame(size = 1e6, row.names = c(...)[1]), .package = "base"
  )
  expect_equal(estimate_bam_sort_ram(small), 30e9) # floor, not 1e6*4
})

test_that("estimate_bam_sort_ram_from_error parses STAR's real error text, NA otherwise", {
  real_error <- paste0(
    "EXITING because of fatal ERROR: not enough memory for BAM sorting:\n",
    "SOLUTION: re-run STAR with at least --limitBAMsortRAM 62884790985"
  )
  expect_equal(estimate_bam_sort_ram_from_error(real_error), 62884790985)
  expect_true(is.na(estimate_bam_sort_ram_from_error("some unrelated alignment error")))
})

test_that("bam_sort_ram_decision: plenty free -> proceed", {
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Memory_Total_GB = 500, Memory_Usage_GB = 50)
  )
  d <- bam_sort_ram_decision(80e9) # 80GB needed, 450GB free, 500GB total
  expect_identical(d$decision, "proceed")
})

test_that("bam_sort_ram_decision: free short but total sufficient -> defer", {
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Memory_Total_GB = 500, Memory_Usage_GB = 470)
  )
  d <- bam_sort_ram_decision(80e9) # 80GB needed, only 30GB free, but 450GB usable total
  expect_identical(d$decision, "defer")
})

test_that("bam_sort_ram_decision: exceeds even the safety-capped total -> crash", {
  testthat::local_mocked_bindings(
    get_system_usage = function(...) list(Memory_Total_GB = 500, Memory_Usage_GB = 10)
  )
  d <- bam_sort_ram_decision(480e9) # 480GB needed > 500*0.9=450GB usable total
  expect_identical(d$decision, "crash")
})

test_that("waiting-for-RAM marker: round-trips, and a second mark doesn't reset the clock", {
  config <- fake_config()
  mark_waiting_for_ram(config, "exp1", "SRR001")
  expect_false(is.na(waiting_for_ram_minutes(config, "exp1", "SRR001")))

  first <- readRDS(waiting_for_ram_path(config, "exp1", "SRR001"))
  Sys.sleep(1.1)
  mark_waiting_for_ram(config, "exp1", "SRR001") # should be a no-op, timestamp unchanged
  expect_identical(readRDS(waiting_for_ram_path(config, "exp1", "SRR001")), first)

  clear_waiting_for_ram(config, "exp1", "SRR001")
  expect_true(is.na(waiting_for_ram_minutes(config, "exp1", "SRR001")))
})

test_that("waiting_for_ram_minutes is NA for a sample never marked", {
  config <- fake_config()
  expect_true(is.na(waiting_for_ram_minutes(config, "exp1", "SRR999")))
})
