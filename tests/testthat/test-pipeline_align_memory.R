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

  first <- readRDS(waiting_for_ram_path(config, "exp1", "SRR001"))$first
  Sys.sleep(1.1)
  mark_waiting_for_ram(config, "exp1", "SRR001") # continuous re-check: first-defer time unchanged, last updated
  m <- readRDS(waiting_for_ram_path(config, "exp1", "SRR001"))
  expect_identical(m$first, first)
  expect_true(m$last > first)

  clear_waiting_for_ram(config, "exp1", "SRR001")
  expect_true(is.na(waiting_for_ram_minutes(config, "exp1", "SRR001")))
})

test_that("waiting_for_ram_minutes is NA for a sample never marked", {
  config <- fake_config()
  expect_true(is.na(waiting_for_ram_minutes(config, "exp1", "SRR999")))
})

test_that("waiting-for-RAM: only continuous waiting counts; a long gap restarts the clock", {
  config <- fake_config(); config$ram_wait_gap_minutes <- 30
  t0 <- as.POSIXct("2026-10-22 02:00:00", tz = "UTC")
  mark_waiting_for_ram(config, "e", "R1", now = t0)
  mark_waiting_for_ram(config, "e", "R1", now = t0 + 20 * 60)          # re-checked within the gap
  expect_equal(waiting_for_ram_minutes(config, "e", "R1", now = t0 + 25 * 60), 25)
  # not re-checked for hours (stage loop busy elsewhere): the wait counts as restarted, not as 10 h
  expect_equal(waiting_for_ram_minutes(config, "e", "R1", now = t0 + 10 * 3600), 0)
  mark_waiting_for_ram(config, "e", "R1", now = t0 + 10 * 3600)
  expect_equal(readRDS(waiting_for_ram_path(config, "e", "R1"))$first, t0 + 10 * 3600)
  # legacy marker (single POSIXct) is read as first = last
  saveRDS(t0, waiting_for_ram_path(config, "e", "R2"))
  expect_equal(waiting_for_ram_minutes(config, "e", "R2", now = t0 + 5 * 60), 5)
})

test_that("bam_sort_ram_decision keeps the reserve free and prefers MemAvailable", {
  u <- list(Memory_Total_GB = 500, Memory_Usage_GB = 400, Memory_Available_GB = 150)
  expect_identical(bam_sort_ram_decision(52e9, usage = u)$decision, "proceed")                  # 150 free
  expect_identical(bam_sort_ram_decision(80e9, reserve_gb = 80, usage = u)$decision, "defer")   # 150 - 80 = 70 < 80
  expect_identical(bam_sort_ram_decision(52e9, reserve_gb = 80, usage = u)$decision, "proceed")
  expect_identical(bam_sort_ram_decision(400e9, reserve_gb = 80, usage = u)$decision, "crash")  # > 500*0.9 - 80
})

test_that("ram_reserve_gb: config, then shared policy file, then 0", {
  expect_equal(ram_reserve_gb(list(ram_reserve_gb = 12)), 12)
  pol <- tempfile(fileext = ".json"); writeLines('{"reserve_gb": 80}', pol)
  withr::local_envvar(RB_POLICY = pol)
  expect_equal(ram_reserve_gb(list()), 80)
  withr::local_envvar(RB_POLICY = tempfile())
  expect_equal(ram_reserve_gb(list()), 0)
})

test_that("bam_sort_plan: STAR sort when it fits, Rsamtools sort when STAR's would defer or crash", {
  u <- list(Memory_Total_GB = 500, Memory_Available_GB = 150)
  cfg <- list(ram_reserve_gb = 80)
  p <- bam_sort_plan(52e9, cfg, usage = u, external_available = TRUE)
  expect_identical(p$sorter, "STAR"); expect_identical(p$limit_bytes, 52e9); expect_identical(p$decision$decision, "proceed")
  p <- bam_sort_plan(323e9, cfg, usage = u, external_available = TRUE)      # STAR: defer (323 > 150 - 80 free)
  expect_identical(p$sorter, "Rsamtools"); expect_identical(p$limit_bytes, 16e9); expect_identical(p$decision$decision, "proceed")
  p <- bam_sort_plan(400e9, cfg, usage = u, external_available = TRUE)      # STAR: crash
  expect_identical(p$sorter, "Rsamtools"); expect_identical(p$decision$decision, "proceed")
  p <- bam_sort_plan(400e9, list(ram_reserve_gb = 80, external_bam_sort_gb = 100), usage = u, external_available = TRUE)
  expect_identical(p$decision$decision, "defer")                            # 100 > 150 - 80
  expect_identical(bam_sort_plan(400e9, cfg, usage = u, external_available = FALSE)$decision$decision, "crash")
  expect_identical(bam_sort_plan(400e9, list(ram_reserve_gb = 80, external_bam_sort = FALSE), usage = u,
                                 external_available = TRUE)$sorter, "STAR")
})
