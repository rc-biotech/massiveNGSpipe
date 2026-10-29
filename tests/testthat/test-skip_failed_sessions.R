test_that("studies with errors in selected previous sessions are skipped", {
  proj <- file.path(tempdir(), "skip_failed_test"); unlink(proj, recursive = TRUE)
  s1 <- file.path(proj, "error_logs", "2026-10-28 10:00:00"); s2 <- file.path(proj, "error_logs", "2026-10-29 10:00:00")
  dir.create(s1, recursive = TRUE); dir.create(s2, recursive = TRUE)
  saveRDS("A-homo_sapiens: (pshift) failed", file.path(s2, "A-homo_sapiens.rds"))     # last session
  saveRDS("B-homo_sapiens: (trim) failed", file.path(s1, "B-homo_sapiens.rds"))       # second last
  mk <- function(s) list(organisms = list(hs = list(conf = c(exp = paste0(s, "-homo_sapiens")))))
  pl <- list(A = mk("A"), B = mk("B"), C = mk("C"))
  cfg <- list(project = proj, skip_failed_sessions = 0)
  expect_equal(names(massiveNGSpipe:::drop_studies_failed_in_sessions(pl, cfg)), c("A", "B", "C"))   # default off
  cfg$skip_failed_sessions <- 1
  expect_equal(names(massiveNGSpipe:::drop_studies_failed_in_sessions(pl, cfg)), c("B", "C"))
  cfg$skip_failed_sessions <- 2
  expect_equal(names(massiveNGSpipe:::drop_studies_failed_in_sessions(pl, cfg)), c("A", "C"))
  cfg$skip_failed_sessions <- c(1, 2)
  expect_equal(names(massiveNGSpipe:::drop_studies_failed_in_sessions(pl, cfg)), "C")
})

test_that("skip_failed_sessions accepts a start time (all sessions since then)", {
  proj <- file.path(tempdir(), "skip_failed_since"); unlink(proj, recursive = TRUE)
  for (s in c("2026-10-20 10:00:00.1", "2026-10-25 00:00:34.5", "2026-10-29 11:00:00.2")) dir.create(file.path(proj, "error_logs", s), recursive = TRUE)
  saveRDS("x", file.path(proj, "error_logs", "2026-10-20 10:00:00.1", "OLD-homo_sapiens.rds"))
  saveRDS("x", file.path(proj, "error_logs", "2026-10-25 00:00:34.5", "V1-homo_sapiens.rds"))
  saveRDS("x", file.path(proj, "error_logs", "2026-10-29 11:00:00.2", "V2-homo_sapiens.rds"))
  mk <- function(s) list(organisms = list(hs = list(conf = c(exp = paste0(s, "-homo_sapiens")))))
  pl <- list(OLD = mk("OLD"), V1 = mk("V1"), V2 = mk("V2"), OK = mk("OK"))
  cfg <- list(project = proj, skip_failed_sessions = "2026-10-25")
  expect_equal(names(massiveNGSpipe:::drop_studies_failed_in_sessions(pl, cfg)), c("OLD", "OK"))
})
