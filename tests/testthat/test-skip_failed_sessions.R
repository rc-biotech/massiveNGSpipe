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
