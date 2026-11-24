# Unit tests for R/resources.R: the one thread formula, the step registry, the resource mode file
# (parse, expiry, atomic write, CLI) and the STAR sort threads / scratch estimate.

no_mode_file <- function(env = parent.frame())
  withr::local_options(list(massiveNGSpipe.resource_mode_file = file.path(tempdir(), "no_such_RESOURCE_MODE")), .local_envir = env)

snap <- function(cores = 48, load1 = 0, mem = 500, mode = NULL)
  list(cores = cores, load1 = load1, mem_available_gb = mem, reserve_gb = 80, mode = mode)

mode_lines <- function(mode = "conserve", set_at = "2026-11-24T10:00:00Z", expires_at = "2026-11-24T14:00:00Z", ...)
  c(paste0("mode: ", mode), paste0("set_at: ", set_at), paste0("expires_at: ", expires_at), ...)
now <- as.POSIXct("2026-11-24 12:00:00", tz = "UTC")

test_that("resource_mode_parse: valid, malformed, expired, future-dated, too long", {
  p <- resource_mode_parse(mode_lines(), now = now, max_hours = 6)
  expect_true(p$ok); expect_identical(p$mode$mode, "conserve"); expect_identical(p$mode$factor, 0.5)
  expect_true(is.na(p$mode$max_mem_gb))
  p <- resource_mode_parse(c(mode_lines("minimal"), "max_mem_gb: 40", "max_threads: 8", "reason: VM stall"), now = now, max_hours = 6)
  expect_true(p$ok); expect_identical(p$mode$max_mem_gb, 40); expect_identical(p$mode$max_threads, 8)
  expect_identical(p$mode$reason, "VM stall")
  expect_match(resource_mode_parse(mode_lines("wild"), now = now)$why, "malformed")
  expect_match(resource_mode_parse("garbage", now = now)$why, "malformed")
  expect_match(resource_mode_parse(character(), now = now)$why, "malformed")
  expect_match(resource_mode_parse("mode: conserve", now = now)$why, "malformed")         # no times
  expect_match(resource_mode_parse(c(mode_lines(), "max_threads: -1"), now = now)$why, "malformed")
  expect_identical(resource_mode_parse(mode_lines(expires_at = "2026-11-24T11:59:00Z"), now = now)$why, "expired")
  expect_match(resource_mode_parse(mode_lines(set_at = "2026-11-24T12:30:00Z", expires_at = "2026-11-24T13:00:00Z"),
                                   now = now, skew_min = 10)$why, "future")
  # within the clock-skew tolerance a slightly future set_at is accepted
  expect_true(resource_mode_parse(mode_lines(set_at = "2026-11-24T12:05:00Z"), now = now, skew_min = 10)$ok)
  expect_match(resource_mode_parse(mode_lines(set_at = "2026-11-24T10:00:00Z", expires_at = "2026-11-25T10:00:00Z"),
                                   now = now, max_hours = 6)$why, "more than 6 h")
})

test_that("resource_mode_read: missing / disabled = NULL silently; invalid = NULL with one warning", {
  withr::local_options(list(massiveNGSpipe.resource_mode_file = "none"))
  expect_null(resource_mode_file()); expect_null(resource_mode_read())
  f <- tempfile("RESOURCE_MODE_")
  expect_silent(expect_null(resource_mode_read(f)))
  writeLines(mode_lines(expires_at = "2026-11-24T11:00:00Z"), f)
  expect_warning(expect_null(resource_mode_read(f, now = now)), "expired")
  expect_silent(expect_null(resource_mode_read(f, now = now)))                     # once per file version
  writeLines("half a fi", f); Sys.setFileTime(f, Sys.time() + 5)
  expect_warning(expect_null(resource_mode_read(f, now = now)), "malformed")
  writeLines(mode_lines(), f); Sys.setFileTime(f, Sys.time() + 10)
  m <- suppressMessages(resource_mode_read(f, now = now))
  expect_identical(m$mode, "conserve")
  # optional rename of an expired file
  withr::local_options(list(massiveNGSpipe.resource_mode_rename_expired = TRUE))
  writeLines(mode_lines(expires_at = "2026-11-24T11:00:00Z"), f); Sys.setFileTime(f, Sys.time() + 15)
  suppressWarnings(resource_mode_read(f, now = now))
  expect_false(file.exists(f)); expect_true(file.exists(paste0(f, ".expired")))
})

test_that("resource_mode_set writes atomically, refresh / clear / show work, limits are enforced", {
  d <- withr::local_tempdir(); f <- file.path(d, "RESOURCE_MODE")
  t0 <- Sys.time()
  suppressMessages(resource_mode_set("minimal", hours = 1, reason = "test\nline", max_mem_gb = 30, set_by = "tester", path = f, now = t0))
  expect_identical(list.files(d, all.files = TRUE, no.. = TRUE), c("RESOURCE_MODE", "RESOURCE_MODE.log"))  # no temp left
  m <- suppressMessages(resource_mode_read(f, now = t0 + 60))
  expect_identical(m$mode, "minimal"); expect_identical(m$max_mem_gb, 30); expect_identical(m$reason, "test line")
  expect_null(suppressWarnings(resource_mode_read(f, now = t0 + 3601)))             # expired after 1 h
  suppressMessages(resource_mode_refresh(hours = 2, path = f, now = t0 + 1800))     # heartbeat
  expect_identical(suppressMessages(resource_mode_read(f, now = t0 + 3601))$mode, "minimal")
  expect_error(resource_mode_set("conserve", hours = 7, reason = "x", path = f), "hours")
  expect_error(resource_mode_set("conserve", hours = 1, reason = "", path = f), "reason")
  out <- capture.output(suppressMessages(resource_mode_show(f, now = t0 + 1900)))
  expect_true(any(grepl("minimal", out)))
  expect_false(suppressMessages(resource_mode_clear(only_set_by = "someone_else", path = f, now = t0 + 1900)))
  expect_true(file.exists(f))
  expect_true(suppressMessages(resource_mode_clear(only_set_by = "tester", path = f, now = t0 + 1900)))
  expect_false(file.exists(f))
  expect_error(resource_mode_set("conserve", hours = 1, reason = "x", path = file.path(d, "nodir", "RM")), "does not exist")
})

test_that("resource_threads_formula: load-based bounds 1..max, min floor on CPU only", {
  f <- resource_threads_formula
  expect_identical(f(16, min = 4, cores = 48, load1 = 0, mem_gb = 500), 16L)            # idle: max
  expect_identical(f(16, min = 4, cores = 48, load1 = 40, mem_gb = 500), 8L)            # 8 idle cores
  expect_identical(f(16, min = 4, cores = 48, load1 = 80, mem_gb = 500), 4L)            # overloaded: floor
  expect_identical(f(3, min = 4, cores = 48, load1 = 80, mem_gb = 500), 3L)             # floor never above max
  expect_identical(f(16, min = 4, cores = 48, load1 = 0, mem_gb = 90, per_worker_gb = 2), 5L)  # (90-80)/2
  expect_identical(f(16, min = 4, cores = 48, load1 = 0, mem_gb = 70, per_worker_gb = 2), 1L)  # below reserve: 1
  expect_identical(f(16, min = 4, cores = 48, load1 = 0, mem_gb = 500, n = 3), 3L)
  expect_identical(f(16, cores = NA, load1 = NA, fallback = 16), 16L)                   # no /proc: fallback
  expect_identical(f(16, cores = NA, load1 = NA, fallback = 1), 1L)
  expect_identical(f(16, cores_per_worker = 0, cores = 48, load1 = 1000), 16L)         # no CPU check
  for (mx in c(1, 2, 7, 46)) for (ld in c(0, 20, 100)) for (mem in c(50, 100, 600)) {
    x <- f(mx, min = 4, cores = 48, load1 = ld, mem_gb = mem, per_worker_gb = 4)
    expect_true(x >= 1L && x <= mx)
  }
})

test_that("mode caps never scale up: result with a mode <= result without, for every step", {
  modes <- list(list(mode = "normal", factor = 1, max_threads = NA, max_mem_gb = NA),
                list(mode = "normal", factor = 1, max_threads = 3, max_mem_gb = NA),
                list(mode = "conserve", factor = 0.5, max_threads = NA, max_mem_gb = NA),
                list(mode = "minimal", factor = 0.25, max_threads = NA, max_mem_gb = 10))
  for (st in resource_steps()$step) for (mx in c(1, 4, 16, 46)) for (ld in c(0, 30, 90)) {
    base <- resource_threads(st, mx, snapshot = snap(load1 = ld))
    expect_true(base >= 1 && base <= mx, info = st)
    for (m in modes) {
      capped <- resource_threads(st, mx, snapshot = snap(load1 = ld, mode = m))
      expect_true(capped >= 1 && capped <= base, info = paste(st, m$mode, mx, ld))
    }
  }
  m <- list(mode = "conserve", factor = 0.5, max_threads = NA, max_mem_gb = NA)
  expect_identical(resource_threads("star", 46, snapshot = snap(mode = m)), 23L)
  m$max_threads <- 10
  expect_identical(resource_threads("star", 46, snapshot = snap(mode = m)), 10L)
  m <- list(mode = "normal", factor = 1, max_threads = NA, max_mem_gb = 16)
  expect_identical(resource_threads("trim_workers", 8, snapshot = snap(mode = m)), 2L)   # 16 GB / 8 GB per worker
})

test_that("dynamic_thread_count = FALSE returns the configured value exactly, ignoring load and mode", {
  m <- list(mode = "minimal", factor = 0.25, max_threads = 1, max_mem_gb = 1)
  for (st in c(resource_steps()$step, "unregistered_step")) for (v in list(1L, 8, 46L, 32)) {
    expect_identical(resource_threads(st, v, config = list(dynamic_thread_count = FALSE), snapshot = snap(load1 = 500, mem = 1, mode = m)), v)
  }
  withr::local_options(list(massiveNGSpipe.dynamic_thread_count = FALSE))
  expect_identical(resource_threads("star", 46), 46)
  # every call site through bpparam_from_config: the configured count
  config <- fake_config(extra = list(dynamic_thread_count = FALSE, thread_type = BiocParallel::MulticoreParam,
                                     threads = list(main = 3, default = 7, trim = 5, pcounts = 2),
                                     parallel_conf = list(log = TRUE, logdir = tempfile("bplog_"), jobname = "x", stop.on.error = TRUE)))
  for (st in c("default", "trim", "pcounts")) expect_identical(BiocParallel::bpnworkers(bpparam_from_config(config, st)),
                                                               as.integer(config$threads[[st]]))
})

test_that("no mode file: every registered call site is load-based within 1..max; idle host = configured value", {
  no_mode_file()
  for (st in c(resource_steps()$step, "unregistered_step")) for (v in c(1L, 5L, 16L, 46L)) {
    x <- resource_threads(st, v)
    expect_true(x >= 1L && x <= v, info = st)
    # an idle host with plenty of memory gives the configured value (today's behaviour)
    expect_identical(resource_threads(st, v, snapshot = snap(cores = 512, load1 = 0, mem = 1e5)), v, info = st)
  }
  expect_identical(resource_threads("default", NULL), NULL)
})

test_that("legacy dynamic helpers keep their behaviour and gain only the mode cap", {
  no_mode_file()
  withr::local_options(list(massiveNGSpipe.trim_workers = NULL))
  expect_identical(trim_workers_dynamic(20, 48, 28, 500), 5L)
  withr::local_options(list(massiveNGSpipe.resource_mode_file = f <- tempfile("RM_")))
  suppressMessages(resource_mode_set("conserve", hours = 1, reason = "t", path = f))
  expect_identical(suppressMessages(trim_workers_dynamic(20, 48, 0, 500)), 4L)          # 8 * 0.5
  withr::local_options(list(massiveNGSpipe.trim_workers = 6))
  expect_identical(suppressMessages(trim_workers_dynamic(20, 48, 60, 500)), 3L)         # fixed 6, capped x0.5
  withr::local_options(list(ORFik.derive_threads = 10))
  expect_identical(suppressMessages(derive_threads_dynamic()), 5L)
})

test_that("STAR sort threads: 16 with headroom, 6 without, never above mapping threads; extra args", {
  no_mode_file()
  plan <- function(free, needed = 30, sorter = "STAR") list(sorter = sorter, decision = list(free_gb = free, needed_gb = needed))
  s <- snap(cores = 512)
  expect_identical(star_sort_threads(list(), plan(300), 46, snapshot = s), 16L)
  expect_identical(star_sort_threads(list(), plan(40), 46, snapshot = s), 6L)          # < 14 GB headroom
  expect_identical(star_sort_threads(list(), plan(300), 8, snapshot = s), 8L)
  expect_identical(star_sort_threads(list(star_sort_threads = 6), plan(300), 46, snapshot = s), 6L)
  expect_true(is.na(star_sort_threads(list(), plan(300, sorter = "Rsamtools"), 46, snapshot = s)))
  expect_identical(star_sort_threads(list(), list(sorter = "STAR", decision = list(decision = "proceed")), 46, snapshot = s), 6L)
  expect_identical(star_sort_threads(list(dynamic_thread_count = FALSE), plan(300), 46), 16L)
  expect_identical(star_extra_args("", 16L, 46L), "--outBAMsortingThreadN 16")
  expect_identical(star_extra_args("--limitOutSJcollapsed 5000000", 16L, 46L), "--limitOutSJcollapsed 5000000 --outBAMsortingThreadN 16")
  expect_identical(star_extra_args("--outBAMsortingThreadN 4", 16L, 46L), "--outBAMsortingThreadN 4")   # user wins
  expect_identical(star_extra_args("", 6L, 46L), "")                                    # STAR default: unchanged
  expect_identical(star_extra_args("", NA_integer_, 46L), "")
  expect_identical(star_extra_args("", 4L, 46L), "--outBAMsortingThreadN 4")
})

test_that("align scratch estimate: ~1 GB per million reads + input + margin", {
  cfg <- list(ssd_min_free_gb = 50)
  # benchmark library L: 0.48 GB gzipped, 67.9 M reads, needed 54 GB temp
  need <- align_scratch_need_gb("x.fasta.gz", cfg, size = 0.48e9)
  expect_true(need > 54 + 50 && need < 130)
  expect_equal(align_scratch_need_gb("x.fastq", list(ssd_min_free_gb = 0, scratch_gb_per_million_reads = 1), size = 3e9),
               3e9 / 1024^3 + 100)
  expect_equal(align_scratch_need_gb("x.fastq", cfg, size = NA), 50)
})

test_that("bam_sort_plan: an active max_mem_gb mode sends a large sort to the external sorter within the cap", {
  f <- tempfile("RM_")
  withr::local_options(list(massiveNGSpipe.resource_mode_file = f))
  usage <- list(Memory_Total_GB = 500, Memory_Available_GB = 400, Memory_Usage_GB = 100)
  cfg <- list(ram_reserve_gb = 80)
  expect_identical(bam_sort_plan(60e9, cfg, usage = usage, external_available = TRUE)$sorter, "STAR")
  suppressMessages(resource_mode_set("normal", hours = 1, reason = "t", max_mem_gb = 40, path = f))
  p <- suppressMessages(bam_sort_plan(60e9, cfg, usage = usage, external_available = TRUE))
  expect_identical(p$sorter, "Rsamtools"); expect_identical(p$limit_bytes, 16e9)
  suppressMessages(resource_mode_set("minimal", hours = 1, reason = "t", path = f))
  expect_identical(suppressMessages(bam_sort_plan(600e9, cfg, usage = usage, external_available = TRUE))$limit_bytes, 4e9)
})

test_that("mngsp resources set / show / clear", {
  f <- tempfile("RM_")
  expect_identical(suppressMessages(mngsp_cli(c("resources", "set", "conserve", "--reason", "VM stall", "--hours", "1",
                                                "--set-by", "cli-test", "--mode-file", f))), 0L)
  expect_identical(suppressMessages(resource_mode_read(f))$set_by, "cli-test")
  out <- capture.output(suppressMessages(mngsp_cli(c("resources", "show", "--mode-file", f))))
  expect_true(any(grepl("conserve", out)))
  expect_identical(suppressMessages(mngsp_cli(c("resources", "clear", "--mode-file", f))), 0L)
  expect_false(file.exists(f))
  expect_identical(suppressMessages(mngsp_cli(c("resources", "set", "conserve", "--hours", "9", "--reason", "x", "--mode-file", f))), 1L)
})
