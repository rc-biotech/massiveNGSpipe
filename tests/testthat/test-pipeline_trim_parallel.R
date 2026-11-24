library(data.table)
# Parallel trim stage (R/pipeline_trim_parallel.R): worker counts, the fork pool, and identical outputs to serial

test_that("trim_workers_dynamic: fixed option wins, CPU and memory limits, floor and cap", {
  withr::with_options(list(massiveNGSpipe.trim_workers = 3), expect_identical(trim_workers_dynamic(10, 48, 0, 500), 3L))
  withr::with_options(list(massiveNGSpipe.trim_workers = 0), expect_identical(trim_workers_dynamic(10, 48, 0, 500), 1L))
  withr::with_options(list(massiveNGSpipe.trim_workers = NULL), {
    expect_identical(trim_workers_dynamic(20, 48, 0, 500), 8L)       # idle host: cap 8
    expect_identical(trim_workers_dynamic(20, 48, 28, 500), 5L)      # 20 idle cores / 4 per worker
    expect_identical(trim_workers_dynamic(20, 48, 60, 500), 2L)      # overloaded host: floor 2
    expect_identical(trim_workers_dynamic(3, 48, 0, 500), 3L)        # never more than libraries
    expect_identical(trim_workers_dynamic(20, 48, 0, 104), 3L)       # (104 - 80 GB reserve) / 8 GB
    expect_identical(trim_workers_dynamic(20, 48, 60, 70), 1L)       # below the reserve: serial, floor does not apply
    expect_identical(trim_workers_dynamic(20, 48, NA, 500), 1L)      # no /proc: serial
  })
  withr::with_options(list(massiveNGSpipe.trim_workers = NULL, massiveNGSpipe.trim_max_workers = 12,
                           massiveNGSpipe.trim_gb_per_worker = 20),
                      expect_identical(trim_workers_dynamic(20, 48, 0, 300), 11L))
})

test_that("trim method, fixed workers and stage slots come from the config, else options", {
  expect_identical(trim_parallel_method(list()), "staged")
  expect_identical(trim_parallel_method(list(trim_parallel = FALSE)), "serial")
  expect_identical(trim_parallel_method(list(trim_parallel = "naive")), "naive")
  withr::with_options(list(massiveNGSpipe.trim_parallel = "serial"), expect_identical(trim_parallel_method(list()), "serial"))
  expect_error(trim_parallel_method(list(trim_parallel = "fast")))
  expect_null(trim_workers_fixed(list()))
  expect_identical(trim_workers_fixed(list(trim_workers = 4)), 4L)
  withr::with_options(list(massiveNGSpipe.trim_workers = 2), expect_identical(trim_workers_fixed(list()), 2L))
  expect_identical(trim_slots(4), c(detect = 4L, trim = 2L))
  expect_identical(trim_slots(1), c(detect = 1L, trim = 1L))
  expect_identical(trim_slots(6, list(trim_stage_slots = c(detect = 2, trim = 1))), c(detect = 2L, trim = 1L))
})

test_that("run_fork_pool chains stages in order, respects slots and returns results by item", {
  skip_on_os("windows")
  d <- tempfile("pool_"); dir.create(d)
  funs <- list(a = function(i) { Sys.sleep(0.2); i * 10 },
               b = function(i, prev) { file.create(file.path(d, paste0("run_", i, "_", Sys.getpid()))); Sys.sleep(0.3)
                 n_now <- length(list.files(d, "^run_")); unlink(file.path(d, paste0("run_", i, "_", Sys.getpid())))
                 c(prev + 1, n_now) })
  res <- run_fork_pool(6, funs, slots = c(a = 3, b = 2), run_ids = paste0("R", 1:6), poll = 0.05)
  expect_identical(vapply(res, `[[`, numeric(1), 1), c(11, 21, 31, 41, 51, 61))
  expect_true(max(vapply(res, `[[`, numeric(1), 2)) <= 2)              # never more than 2 jobs in stage b
})

test_that("run_fork_pool starts items in start_order (largest library first)", {
  skip_on_os("windows")
  d <- tempfile("poolorder_"); dir.create(d)
  f <- list(all = function(i) { cat(i, "\n", file = file.path(d, "order"), append = TRUE); i })
  res <- run_fork_pool(4, f, slots = c(all = 1), start_order = c(3L, 1L, 4L, 2L), poll = 0.05)
  expect_identical(unlist(res), 1:4)                                  # results stay by item
  expect_identical(as.integer(readLines(file.path(d, "order"))), c(3L, 1L, 4L, 2L))
})

test_that("run_fork_pool: a failure stops new items, running ones finish, all errors are reported", {
  skip_on_os("windows")
  d <- tempfile("poolerr_"); dir.create(d)
  funs <- list(all = function(i) {
    if (i == 2) stop("broken input")
    Sys.sleep(if (i == 1) 1 else 0.1)
    file.create(file.path(d, paste0("done_", i))); i })
  expect_error(run_fork_pool(6, funs, slots = c(all = 2), run_ids = paste0("R", 1:6), poll = 0.05),
               "R2 \\(all\\): .*broken input")
  done <- sort(list.files(d))
  expect_true("done_1" %in% done)                                     # started before the failure: finished
  expect_false(any(c("done_4", "done_5", "done_6") %in% done))        # not started after it
})

# ---- end to end on a small fixture: serial and parallel give the same files and markers ----
make_lib <- function(file, n = 4000, seed = 1, umi5 = 0) {
  set.seed(seed)
  rnd <- function(k) paste(sample(c("A", "C", "G", "T"), k, replace = TRUE), collapse = "")
  adapter <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
  inserts <- vapply(1:300, function(i) rnd(sample(26:34, 1)), "")       # duplicated footprints, as in Ribo-seq
  lines <- character(4 * n)
  for (i in seq_len(n)) {
    seq <- substr(paste0(rnd(umi5), sample(inserts, 1), adapter, strrep("A", 80)), 1, 75)
    lines[(4 * i - 3):(4 * i)] <- c(sprintf("@r%d", i), seq, "+", strrep("F", 75))
  }
  con <- gzfile(file, "w"); writeLines(lines, con); close(con)
  file
}

trim_fixture <- function(root, method, workers, broken = NULL) {
  raw <- file.path(root, "raw"); target <- file.path(root, "proc"); trimmed <- file.path(target, "trim")
  dir.create(raw, recursive = TRUE); dir.create(trimmed, recursive = TRUE)
  runs <- data.table(Run = sprintf("SRRT%02d", 1:4), LIBRARYTYPE = c("RFP", "RFP", "RNA", "RNA"), LibraryLayout = "SINGLE")
  files <- lapply(seq_len(nrow(runs)), function(i) {
    f <- make_lib(file.path(raw, paste0(runs$Run[i], ".fastq.gz")), seed = i, umi5 = if (i == 2) 4 else 0)
    if (identical(runs$Run[i], broken)) {   # colour-space reads: refused at once by assert_base_space_fastq()
      f <- file.path(raw, paste0(runs$Run[i], "_cs.fastq"))
      writeLines(rep(c("@r", "T0123012301230123012301230123", "+", strrep("F", 29)), 50), f)
    }
    c(f, NA)
  })
  config <- fake_config(project = file.path(root, "project"))
  ctx <- list(all_files = files, runs = runs, mode = "local", target_dir = target, trimmed_dir = trimmed, source_dir = raw,
              config = config, experiment = "fixture-homo_sapiens")
  withr::with_options(list(massiveNGSpipe.fast_trim = TRUE, massiveNGSpipe.adapter_detection = "scan",
                           massiveNGSpipe.barcode_detector = "guarded"),
                      trim_samples_run(ctx, method = method, workers = workers))
  ctx
}
collapsed_set <- function(f) {
  x <- Biostrings::readDNAStringSet(f)
  data.table(seq = as.character(x), n = as.integer(sub(".*_x", "", names(x))))[order(seq)]
}

test_that("staged and naive parallel trim write the same files and markers as serial", {
  skip_on_os("windows")
  skip_if(Sys.which("cc") == "", "no C compiler")
  base <- tempfile("ptrim_")
  s <- trim_fixture(file.path(base, "serial"), "serial", 1)
  for (m in c("staged", "naive")) {
    p <- trim_fixture(file.path(base, m), m, 3)
    for (r in s$runs$Run) {
      a <- file.path(s$trimmed_dir, "SINGLE", paste0("collapsed_trimmed_", r, ".fasta.gz"))
      b <- file.path(p$trimmed_dir, "SINGLE", paste0("collapsed_trimmed_", r, ".fasta.gz"))
      expect_true(file.exists(b))
      expect_identical(collapsed_set(b), collapsed_set(a))
      ja <- jsonlite::fromJSON(file.path(s$trimmed_dir, paste0("report_", r, ".json")))
      jb <- jsonlite::fromJSON(file.path(p$trimmed_dir, paste0("report_", r, ".json")))
      expect_identical(jb$summary, ja$summary)
      expect_identical(jb$filtering_result, ja$filtering_result)
    }
    for (step in c("trim", "collapsed")) {
      va <- sample_flag_values(s$config, step, s$experiment); vb <- sample_flag_values(p$config, step, p$experiment)
      expect_identical(vb[order(names(vb))], va[order(names(va))])
    }
    expect_length(list.files(file.path(p$trimmed_dir, "barcode_detect_sample")), 0)   # detection scratch removed
  }
})

test_that("parallel trim: a failing library fails the stage but finished libraries keep their markers", {
  skip_on_os("windows")
  skip_if(Sys.which("cc") == "", "no C compiler")
  root <- tempfile("ptrim_err_")
  expect_error(trim_fixture(root, "naive", 2, broken = "SRRT01"), "SRRT01")
  config <- fake_config(project = file.path(root, "project"))
  done <- samples_done(config, "trim", "fixture-homo_sapiens")
  expect_false("SRRT01" %in% done)                                   # no marker: retried on the next attempt
  expect_true(length(done) >= 1)                                     # jobs running at the failure were finished
  for (r in done) expect_true(file.exists(file.path(root, "proc", "trim", "SINGLE", paste0("collapsed_trimmed_", r, ".fasta.gz"))))
})
