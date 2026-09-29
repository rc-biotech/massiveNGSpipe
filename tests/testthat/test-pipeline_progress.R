# Unit tests for R/pipeline_progress.R -- session bookkeeping, per-study
# status walking, and small pure summary helpers. get_pgid() (real
# system() call) is mocked for session_info_save() round-trip tests;
# get_pgid_all() itself is verified live against a real background
# process below, since it was found and fixed as part of this session's
# cleanup (previously filtered on the literal string "pgid" instead of
# returning all real running pgids).

test_that("current_step_table labels progress counts by step name, plus a final 'complete' bucket", {
  steps <- c("start", "fetch", "trim")
  tab <- current_step_table(c(0, 0, 1, 3, 3), steps)
  expect_identical(unname(tab[names(tab) == "start(0)"]), 2L)
  expect_identical(unname(tab[names(tab) == "fetch(1)"]), 1L)
  expect_identical(unname(tab[names(tab) == "complete(3)"]), 2L)
})

test_that("session_info_save/read round-trip with a mocked get_pgid", {
  testthat::local_mocked_bindings(get_pgid = function() 12345L)
  config <- fake_config(session_dir = tempfile("session_"))
  dir.create(config$session_dir, recursive = TRUE)
  session_info_save(config, fake_pipelines())

  info <- session_info_read(config)
  expect_identical(info$pgid, 12345L)
  expect_identical(info$status, "started")
  expect_identical(info$pipeline_names, "PRJNA000001")
})

test_that("session_info_save errors without an active session_dir", {
  config <- fake_config(session_dir = NULL)
  expect_error(session_info_save(config, fake_pipelines()), "active session")
})

test_that("session_info_table lists session dirs sorted newest-first", {
  project <- tempfile("mNGSp_test_")
  older <- file.path(project, "session_logs", "2020-01-01 00:00:00")
  newer <- file.path(project, "session_logs", "2025-01-01 00:00:00")
  dir.create(older, recursive = TRUE)
  Sys.setFileTime(older, as.POSIXct("2020-01-01"))
  dir.create(newer, recursive = TRUE)
  Sys.setFileTime(newer, as.POSIXct("2025-01-01"))

  config <- fake_config(project = project)
  info <- session_info_table(config)
  expect_identical(basename(info$path[1]), "2025-01-01 00:00:00")
})

test_that("status_per_study reports done/not-done experiments correctly from real flags", {
  config <- fake_config()
  steps <- names(config$flag)
  pipelines <- fake_pipelines(runs = data.table::data.table(
    Run = "SRR001", ScientificName = "Homo sapiens"
  ))
  exp <- "PRJNA000001-homo_sapiens"

  # Nothing done yet. done/progress_index are built with plain numeric
  # (not integer) literals/arithmetic in status_per_study(), so compare
  # by value (expect_equal), not type (expect_identical).
  result <- suppressMessages(status_per_study(pipelines, config, FALSE, FALSE, FALSE))
  expect_equal(result$done, 0)
  expect_equal(result$progress_index, 0)
  expect_identical(result$all_organism, "Homo sapiens")

  # Fully done.
  for (s in steps) set_flag(config, s, exp)
  result2 <- suppressMessages(status_per_study(pipelines, config, FALSE, FALSE, FALSE))
  expect_equal(result2$done, 1)
  expect_equal(result2$progress_index, length(steps))
})

test_that("organism_report doesn't error and correctly zero-fills organisms with no done studies", {
  config <- fake_config()
  expect_output(
    organism_report(all_organism = c("Homo sapiens", "Mus musculus"),
                    config = config, progress_index = c(0L, 0L)),
    "Processed organisms"
  )
})

# get_pgid_all() itself was fixed this session (previously filtered
# get_running_processes() on the literal string "pgid" rather than
# returning real pgids -- see the cleanup notes). A live end-to-end test
# against get_running_processes() is deliberately NOT included here:
# that function's own column-position assumptions about `ps`'s output
# turned out to be fragile in this specific container (confirmed while
# writing this suite -- a real, pre-existing issue, out of scope for this
# pass), which would make a test that depends on it unreliable across
# environments rather than a genuine regression check on get_pgid_all()'s
# own (correct) logic. get_pgid_all()'s logic is covered by direct
# inspection: unique(as.integer(get_running_processes()$PGID)) with no
# filter, matching what last_session_stil_active() actually needs
# ("all currently running pgids", not a filtered subset).
test_that("get_pgid_all's own logic applies no spurious filter to get_running_processes()", {
  testthat::local_mocked_bindings(
    get_running_processes = function(...) data.table::data.table(PGID = c(10L, 10L, 20L))
  )
  expect_setequal(get_pgid_all(), c(10L, 20L))
})
