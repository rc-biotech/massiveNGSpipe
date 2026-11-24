# Command line interface (R/cli.R, R/cli_study_list.R): parsing, settings merge, dry-run planning, reset safety.

cli_fixture <- function() {
  root <- tempfile("mngsp_cli_")
  data <- c(fastq = file.path(root, "fastq"), bam = file.path(root, "bam"), ref = file.path(root, "ref"),
            exp = file.path(root, "exp"))
  for (d in data) dir.create(d, recursive = TRUE)
  project <- file.path(root, "project"); dir.create(project)
  final <- data.table::data.table(
    study_accession = c("PRJNA1", "PRJNA1", "GSE2", "GSE2", "PRJNA3", "PRJNA1028780", "PRJNA4"),
    BioProject = c("PRJNA1", "PRJNA1", "PRJNA2", "PRJNA2", "PRJNA3", "PRJNA1028780", "PRJNA4"),
    Run = c("SRR1", "SRR2", "SRR3", "SRR4", "SRR5", "SRR6", "SRR7"),
    ScientificName = c("Homo sapiens", "Homo sapiens", "Homo sapiens", "Mus musculus", "Homo sapiens",
                       "Homo sapiens", "Homo sapiens"),
    LibraryLayout = "SINGLE")
  data.table::fwrite(final, file.path(project, "FINAL_LIST.csv"))
  # complete human reference, no mouse reference
  hs <- file.path(data["ref"], "homo_sapiens")
  dir.create(file.path(hs, "STAR_index", "genomeDir"), recursive = TRUE)
  saveRDS(TRUE, file.path(hs, "outputs.rds")); file.create(file.path(hs, "STAR_index", "genomeDir", "SA"))
  writeLines(c("# test list", "PRJNA1", "PRJNA2", "PRJNA3", "PRJNA1028780", "PRJNA4", "PRJNA_MISSING"),
             file.path(root, "studies.txt"))
  list(root = root, data = data, project = project, list = file.path(root, "studies.txt"))
}

local_cli_env <- function(fx, env = parent.frame()) {
  testthat::local_mocked_bindings(config = function(...) fx$data, .package = "ORFik", .env = env)
  testthat::local_mocked_bindings(blas_set_num_threads = function(...) NULL, omp_set_num_threads = function(...) NULL,
                                  .package = "RhpcBLASctl", .env = env)
  withr::local_options(massiveNGSpipe.barcode_detector = "guarded", .local_envir = env)
}

base_flags <- function(fx) c("--project-dir", fx$project, "--no-discord", "--no-preflight", "--google-sheet", "none",
                             "--ssd-scratch-dir", "none")

test_that("cli_parse handles subcommands, dashes/underscores, =values, booleans and threads per step", {
  p <- cli_parse(c("run", "list.txt", "--apply", "--manual_dir", "m", "--skip-failed=2", "--no-preflight",
                   "--reset", "TRUE", "--threads", "12", "--threads-trim", "4", "--accepted-lengths-rpf", "20,25:27"))
  expect_identical(p$command, "run")
  expect_identical(p$positional, "list.txt")
  expect_true(p$opts$apply); expect_true(p$opts$reset); expect_false(p$opts$preflight)
  expect_identical(p$opts$manual_dir, "m"); expect_identical(p$opts$skip_failed, "2")
  expect_identical(p$opts$threads, 12L); expect_identical(p$opts$threads_trim, 4L)
  expect_identical(p$opts$accepted_lengths_rpf, c(20L, 25L, 26L, 27L))
  expect_identical(cli_parse(character())$command, "help")
  expect_identical(cli_parse("--help")$command, "help")
  expect_true(cli_parse(c("run", "--help"))$help)
  expect_identical(cli_parse(c("retry", "A-homo_sapiens", "B-mus_musculus"))$positional, c("A-homo_sapiens", "B-mus_musculus"))
})

test_that("cli_parse gives helpful errors", {
  expect_error(cli_parse(c("rnu", "x")), "did you mean: run")
  expect_error(cli_parse(c("run", "x", "--manul-dir", "d")), "did you mean --manual-dir")
  expect_error(cli_parse(c("status", "--apply")), "does not apply to 'mngsp status'")
  expect_error(cli_parse(c("run", "x", "--keep")), "--keep needs a value")
  expect_error(cli_parse(c("run", "x", "--threads", "many")), "--threads expects a whole number")
  expect_error(cli_parse(c("run", "x", "--mode", "cloud")), "one of online, local")
  expect_error(cli_parse(c("--apply", "run")), "the command comes first")
})

test_that("settings merge: defaults < mngsp.conf < <list>.conf < '#:' lines < --config < flags, paths relative to their file", {
  d <- tempfile("mngsp_conf_"); dir.create(d)
  writeLines("PRJNA1", file.path(d, "keep.txt"))
  writeLines(c("# folder settings", "threads: 4", "max_unprocessed_downloads: 60", "keep: keep.txt",
               "split-unique-mappers: TRUE", "file_prefix: folder_"), file.path(d, "mngsp.conf"))
  writeLines(c("threads = 5", "wait: 30"), file.path(d, "wave.conf"))
  writeLines(c("#: threads: 6", "#: subset_suffix: _v4", "PRJNA1"), file.path(d, "wave.txt"))
  writeLines("threads: 8", file.path(d, "extra.conf"))
  r <- cli_resolve(cli_parse(c("run", file.path(d, "wave.txt"), "--config", file.path(d, "extra.conf"))), cwd = d)
  expect_identical(r$values$threads, 8L)
  expect_match(r$source[["threads"]], "extra.conf:1$")
  expect_identical(r$values$wait, 30)
  expect_identical(r$values$subset_suffix, "_v4")
  expect_identical(r$values$max_unprocessed_downloads, 60L)
  expect_identical(r$values$keep, normalizePath(file.path(d, "keep.txt")))
  expect_identical(r$values$file_prefix, "folder_")
  expect_identical(r$values$skip_dbgap, TRUE)
  expect_identical(r$source[["skip_dbgap"]], "default")
  r2 <- cli_resolve(cli_parse(c("run", "wave.txt", "--config", "extra.conf", "--threads", "12", "--file-prefix", "x_")), cwd = d)
  expect_identical(r2$values$threads, 12L)
  expect_identical(r2$source[["threads"]], "command line")
  expect_identical(r2$values$file_prefix, "x_")
  r3 <- cli_resolve(cli_parse(c("run", "wave.txt", "--no-config")), cwd = d)
  expect_null(r3$values$threads)
  writeLines("apply: TRUE", file.path(d, "bad.conf"))
  expect_error(cli_resolve(cli_parse(c("run", "wave.txt", "--config", "bad.conf")), cwd = d), "only be given on the command line")
  writeLines("thread: 3", file.path(d, "typo.conf"))
  expect_error(cli_resolve(cli_parse(c("run", "wave.txt", "--config", "typo.conf")), cwd = d), "typo.conf:1: unknown setting 'thread'")
  expect_error(cli_resolve(cli_parse(c("run", "nolist.txt")), cwd = d), "study list not found")
  expect_error(cli_resolve(cli_parse(c("run", "wave.txt", "--keep", "nofile")), cwd = d), "--keep: file not found")
  expect_error(cli_resolve(cli_parse(c("plan")), cwd = d), "needs a study list")
})

test_that("cli_build_config passes settings to pipeline_config and config fields; defaults stay pipeline_config's", {
  fx <- cli_fixture(); local_cli_env(fx)
  r <- cli_resolve(cli_parse(c("run", fx$list, base_flags(fx), "--threads", "6", "--threads-trim", "2",
                               "--blas-cap", "4", "--ram-reserve-gb", "80", "--bam-sort-memory-gb", "40",
                               "--split-unique-mappers", "--max-unprocessed-downloads", "60", "--skip-failed", "2026-11-20")))
  config <- cli_build_config(r$values)
  expect_identical(config$project, fx$project)
  expect_null(config$google_url); expect_null(config$discord_webhook); expect_null(config$ssd_scratch_dir)
  expect_identical(config$threads$default, 6L); expect_identical(config$threads$trim, 2L)
  expect_identical(config$threads$pshifted, 4L)
  expect_identical(config$ram_reserve_gb, 80); expect_identical(config$bam_sort_memory_gb, 40)
  expect_true(config$split_unique_mappers); expect_identical(config$max_unprocessed_downloads, 60L)
  expect_identical(config$skip_failed_sessions, "2026-11-20")
  expect_identical(config$barcode_detector, "guarded")
  d <- cli_build_config(cli_resolve(cli_parse(c("run", fx$list, base_flags(fx))))$values)
  pc <- formals(pipeline_config)
  expect_identical(d$split_unique_mappers, pc$split_unique_mappers)
  expect_identical(d$max_unprocessed_downloads, pc$max_unprocessed_downloads)
  expect_identical(d$reuse_shifts_if_existing, pc$reuse_shifts_if_existing)
  expect_identical(d$skip_failed_sessions, pc$skip_failed_sessions)
  expect_true(d$delete_raw_files)   # mode online
})

test_that("mngsp plan: dry-run planning of a tiny study list (keep, dbGaP, organisms, exclude-runs, missing)", {
  fx <- cli_fixture(); local_cli_env(fx)
  writeLines("PRJNA4", file.path(fx$root, "keep.txt"))
  data.table::fwrite(data.table::data.table(BioProject = c("PRJNA1", "PRJNA2", "PRJNA3"),
                                            ScientificName = c("Homo sapiens", "Mus musculus", "Homo sapiens")),
                     file.path(fx$root, "orgs.csv"))
  data.table::fwrite(data.table::data.table(Run = "SRR2", reason = "KEEP FALSE"), file.path(fx$root, "keep_false.csv"))
  dir.create(file.path(fx$project, "flags", "trim"), recursive = TRUE)
  saveRDS(TRUE, file.path(fx$project, "flags", "trim", "PRJNA1-homo_sapiens.rds"))
  expect_message(status <- mngsp_cli(c("plan", fx$list, base_flags(fx), "--keep", file.path(fx$root, "keep.txt"),
                                       "--organisms", file.path(fx$root, "orgs.csv"),
                                       "--exclude-runs", file.path(fx$root, "keep_false.csv")), cwd = fx$root),
                 "Dry run: nothing changed")
  expect_identical(status, 0L)
  plan <- data.table::fread(file.path(fx$project, "mngsp_plan_studies.csv"))
  expect_setequal(plan$experiment, c("PRJNA1-homo_sapiens", "PRJNA2-mus_musculus", "PRJNA3-homo_sapiens"))
  expect_identical(plan[experiment == "PRJNA1-homo_sapiens"]$runs, 1L)   # SRR2 excluded
  expect_identical(plan[experiment == "PRJNA1-homo_sapiens"]$flags_set, 1L)
  expect_identical(plan[experiment == "PRJNA3-homo_sapiens"]$next_step, "start")
  expect_false(plan[experiment == "PRJNA2-mus_musculus"]$reference_complete)
  expect_true(plan[experiment == "PRJNA3-homo_sapiens"]$reference_complete)
  p <- study_list_plan(fx$list, cli_build_config(cli_resolve(cli_parse(c("plan", fx$list, base_flags(fx))))$values))
  expect_identical(p$missing, "PRJNA_MISSING")
  expect_identical(unname(p$counts["dbgap"]), 1L)
  expect_false("PRJNA1028780" %in% p$plan$study)
})

test_that("--reset without --apply changes nothing; study_list_reset clears flags, folder and experiment CSV", {
  fx <- cli_fixture(); local_cli_env(fx)
  args <- c("run", fx$list, base_flags(fx), "--reset")
  mngsp_cli(c("plan", fx$list, base_flags(fx)), cwd = fx$root)       # creates the flag folders
  exp <- "PRJNA1-homo_sapiens"
  saveRDS(TRUE, file.path(fx$project, "flags", "trim", paste0(exp, ".rds")))
  dir.create(file.path(fx$data["bam"], exp, "trim"), recursive = TRUE)
  writeLines("x", file.path(fx$data["bam"], exp, "trim", "barcodes_manual.csv"))
  writeLines("x", file.path(fx$data["exp"], paste0(exp, ".csv")))
  snapshot <- function() {
    f <- list.files(fx$root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
    f <- f[!grepl("mngsp_plan_", f)]
    data.table::data.table(f = f, mtime = file.mtime(f), size = file.size(f))[order(f)]
  }
  before <- snapshot()
  expect_message(status <- mngsp_cli(args, cwd = fx$root), "--reset needs --apply")
  expect_identical(status, 0L)
  expect_identical(snapshot(), before)
  config <- cli_build_config(cli_resolve(cli_parse(args))$values)
  expect_message(study_list_reset(config, exp), "Reset 1 experiments")
  expect_false(file.exists(file.path(fx$project, "flags", "trim", paste0(exp, ".rds"))))
  expect_false(dir.exists(file.path(fx$data["bam"], exp)))
  expect_false(file.exists(file.path(fx$data["exp"], paste0(exp, ".csv"))))
})

test_that("mngsp init writes a settings file that parses back, and mngsp config reports sources", {
  fx <- cli_fixture(); local_cli_env(fx)
  d <- file.path(fx$root, "demo")
  expect_message(mngsp_cli(c("init", d)), "Created")
  expect_true(file.exists(file.path(d, "mngsp.conf"))); expect_true(file.exists(file.path(d, "studies.txt")))
  expect_identical(read_study_list(file.path(d, "studies.txt")), "PRJNA779948")
  expect_message(s <- mngsp_cli(c("init", d)), "exists \\(use --force")
  expect_identical(s, 1L)
  r <- cli_resolve(cli_parse(c("plan", file.path(d, "studies.txt"))), cwd = d)
  expect_identical(r$values$project_dir, normalizePath(file.path(d, "project")))
  expect_identical(r$values$preset, "Ribo-seq")
  out <- capture.output(expect_message(mngsp_cli(c("config", fx$list, base_flags(fx), "--threads", "3")), "Project:"))
  expect_true(any(grepl("threads +3 +command line", out)))
})

test_that("help texts list every command and option", {
  h <- cli_help()
  for (cmd in c("init", "curate", "plan", "run", "status", "errors", "retry", "preflight", "config")) expect_match(h, cmd)
  hr <- cli_help("run")
  expect_match(hr, "--manual-dir <dir>", fixed = TRUE)
  expect_match(hr, "--threads-<step> <n>", fixed = TRUE)
  expect_match(hr, "[default: Ribo-seq]", fixed = TRUE)
  expect_identical(capture.output(s <- mngsp_cli(c("errors", "--help")))[1], "Usage: mngsp errors [--session <n>] [--events <n>]")
  expect_identical(s, 0L)
})

test_that("shell wrapper: mngsp --help and mngsp run --help work", {
  skip_on_cran()
  skip_if(nzchar(Sys.getenv("MNGSP_SKIP_SHELL_TEST")), "MNGSP_SKIP_SHELL_TEST set")
  bin <- system.file("bin", "mngsp", package = "massiveNGSpipe")
  skip_if(!nzchar(bin) || !nzchar(Sys.which("bash")), "no bin/mngsp or bash")
  root <- dirname(dirname(dirname(bin)))
  env <- if (basename(dirname(dirname(bin))) == "inst") paste0("MNGSP_DEV=", root) else character()   # source checkout
  out <- suppressWarnings(system2("bash", c(shQuote(bin), "--help", "2>/dev/null"), stdout = TRUE, env = env))
  expect_null(attr(out, "status"))
  expect_true(any(grepl("^Usage: mngsp <command>", out)))
  out <- suppressWarnings(system2("bash", c(shQuote(bin), "run", "--help", "2>/dev/null"), stdout = TRUE, env = env))
  expect_null(attr(out, "status"))
  expect_true(any(grepl("^Usage: mngsp run <study list>", out)))
  expect_true(any(grepl("--organisms <file>", out, fixed = TRUE)))
  out <- suppressWarnings(system2("bash", c(shQuote(bin), "rnu", "2>&1"), stdout = TRUE, env = env))
  expect_identical(attr(out, "status"), 1L)
})
