# Temp-space claims (R/tmp_space.R): pre-emptive reservations on a shared temp
# drive, stale-claim purge with cleanup, and the stale R tempdir report.

fake_free <- function(gb) function(dir) gb

test_that("tmp_space_claim grants while headroom remains and counts other live claims", {
  d <- withr::local_tempdir()
  id1 <- tmp_space_claim(d, 40, file.path(d, "a"), min_free_gb = 10, free_fun = fake_free(100))
  expect_type(id1, "character")
  # 100 free - 40 outstanding - 40 requested = 20 >= 10
  id2 <- tmp_space_claim(d, 40, file.path(d, "b"), min_free_gb = 10, free_fun = fake_free(100))
  expect_type(id2, "character")
  # 100 - 80 - 20 = 0 < 10: refused
  expect_message(id3 <- tmp_space_claim(d, 20, file.path(d, "c"), min_free_gb = 10, free_fun = fake_free(100)),
                 "refused")
  expect_null(id3)
  tmp_space_release(d, id1); tmp_space_release(d, id2)
  expect_type(tmp_space_claim(d, 20, min_free_gb = 10, free_fun = fake_free(100)), "character")
})

test_that("outstanding GB shrinks as the claim's own files are written", {
  d <- withr::local_tempdir()
  w <- file.path(d, "work"); dir.create(w)
  writeBin(raw(2 * 1024^2), file.path(w, "x"))  # 2 MiB already on disk
  id <- tmp_space_claim(d, 0.005, w, min_free_gb = 0, free_fun = fake_free(1))
  cl <- readRDS(file.path(tmp_space_claim_dir(d), paste0(id, ".rds")))
  expect_lt(claims_outstanding_gb(list(cl)), 0.005)
  expect_gt(claims_outstanding_gb(list(cl)), 0)
})

test_that("tmp_space_release deletes the claimed paths and the claim file", {
  d <- withr::local_tempdir()
  w <- file.path(d, "work"); dir.create(w); writeLines("x", file.path(w, "f"))
  id <- tmp_space_claim(d, 1, w, min_free_gb = 0, free_fun = fake_free(10))
  tmp_space_release(d, id)
  expect_false(dir.exists(w))
  expect_length(list.files(tmp_space_claim_dir(d), "\\.rds$"), 0)
})

test_that("a claim of a dead process is purged with its files, never paths outside dir", {
  d <- withr::local_tempdir()
  outside <- withr::local_tempfile(); writeLines("keep", outside)
  w <- file.path(d, "dead_work"); dir.create(w); writeLines("x", file.path(w, "f"))
  dir.create(tmp_space_claim_dir(d))
  dead_pid <- 2^22 + 12345  # above pid_max: never alive
  saveRDS(list(id = "dead", pid = dead_pid, host = host_name(), gb = 50, paths = c(w, outside),
               time = Sys.time()), file.path(tmp_space_claim_dir(d), "dead.rds"))
  # the dead claim's 50 GB no longer blocks a new claim
  expect_message(id <- tmp_space_claim(d, 50, min_free_gb = 10, free_fun = fake_free(70)), "dead process")
  expect_type(id, "character")
  expect_false(dir.exists(w))
  expect_true(file.exists(outside))
  expect_false(file.exists(file.path(tmp_space_claim_dir(d), "dead.rds")))
})

test_that("claims of other hosts are kept (pids are per container)", {
  d <- withr::local_tempdir(); dir.create(tmp_space_claim_dir(d))
  saveRDS(list(id = "other", pid = 2^22 + 1, host = "another-host", gb = 50, paths = character(),
               time = Sys.time()), file.path(tmp_space_claim_dir(d), "other.rds"))
  expect_null(suppressMessages(tmp_space_claim(d, 30, min_free_gb = 10, free_fun = fake_free(70))))
  expect_true(file.exists(file.path(tmp_space_claim_dir(d), "other.rds")))
})

test_that("tmp_space_claim rejects paths outside dir", {
  d <- withr::local_tempdir()
  expect_error(tmp_space_claim(d, 1, "/etc/passwd", min_free_gb = 0, free_fun = fake_free(10)))
  expect_error(tmp_space_claim(d, 1, d, min_free_gb = 0, free_fun = fake_free(10)))
})

test_that("tmp_cleanup_stale only deletes old, closed dirs holding pipeline artefacts", {
  d <- withr::local_tempdir()
  old_time <- Sys.time() - 5 * 86400
  mk <- function(name, files) {
    p <- file.path(d, name); dir.create(p)
    for (f in files) { dir.create(dirname(file.path(p, f)), showWarnings = FALSE, recursive = TRUE)
      writeLines("x", file.path(p, f)) }
    all <- c(list.files(p, recursive = TRUE, full.names = TRUE, include.dirs = TRUE), p)
    Sys.setFileTime(all, old_time)
    p
  }
  sra <- mk("RtmpAAAA", c("SRR123", "SRR123_1.fastq", "fasterq.tmp.host.1/x"))
  fqc <- mk("RtmpBBBB", c("file1a2b/stdin_fastqc.html", "adapter_candidates.txt"))
  user <- mk("RtmpCCCC", c("my_analysis.qs"))
  fresh <- file.path(d, "RtmpDDDD"); dir.create(fresh); writeLines("x", file.path(fresh, "SRR1"))
  rep <- tmp_cleanup_stale(d, older_than_days = 2, delete = FALSE)
  expect_setequal(rep[candidate == TRUE]$path, c(sra, fqc))
  expect_true(all(dir.exists(c(sra, fqc, user, fresh))))
  rep <- tmp_cleanup_stale(d, older_than_days = 2, delete = TRUE)
  expect_false(any(dir.exists(c(sra, fqc))))
  expect_true(all(dir.exists(c(user, fresh))))
})

test_that("fasterq_tmp_need_gb covers sra copy, scratch and output", {
  d <- withr::local_tempdir()
  writeBin(raw(1024^2), file.path(d, "SRR1"))  # 1 MiB .sra
  se <- fasterq_tmp_need_gb("SRR1", d, PAIRED_END = FALSE)
  pe <- fasterq_tmp_need_gb("SRR1", d, PAIRED_END = TRUE)
  sra_gb <- 1 / 1024
  expect_equal(se, 1.1 * (sra_gb + 2 * 5 * sra_gb))
  expect_gt(pe, se)
  info <- data.table::data.table(Run = "SRR1", spots = 1e6, bases = 5e7, avgLength = 50, LibraryLayout = "SINGLE")
  expect_equal(fasterq_tmp_need_gb("SRR1", d, FALSE, info),
               1.1 * (sra_gb + 2 * estimate_fastq_tmp_gb(info)$total_GiB))
})

test_that("low_disk_hint names a full drive and is NULL otherwise", {
  cfg <- list(config = c(fastq = tempdir(), bam = tempdir()))
  local_mocked_bindings(tmp_free_gb = function(dir) 0.5)
  expect_match(low_disk_hint(cfg), "disk full")
  local_mocked_bindings(tmp_free_gb = function(dir) 500)
  expect_null(low_disk_hint(cfg))
})

test_that("fastqc_parse_adapters removes the FastQC report dir after parsing", {
  report_dir <- file.path(withr::local_tempdir(), "fileabc")
  dir.create(file.path(report_dir, "stdin_fastqc"), recursive = TRUE)
  report <- file.path(report_dir, "stdin_fastqc", "fastqc_data.txt")
  writeLines("dummy", report)
  local_mocked_bindings(run_fastqc = function(...) report,
                        read_fastq_adapter_output = function(p) { expect_true(file.exists(p)); "parsed" })
  expect_identical(fastqc_parse_adapters("x.fastq", 10, "adapters.txt"), "parsed")
  expect_false(dir.exists(report_dir))
})

test_that("error_hint prefers a specific remote error over the generic BiocParallel worker message", {
  msg <- "BiocParallel errors\n 1 remote errors\n 5 unevaluated and other errors\n first remote error: Library contained no periodic or accepted read-lengths"
  expect_match(error_hint(msg), "periodicity")
  expect_match(error_hint("5 unevaluated and other errors"), "worker died")
})
