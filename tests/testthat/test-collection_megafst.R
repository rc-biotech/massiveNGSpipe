# Unit tests for R/collection_megafst.R on small synthetic covRle libraries (no production data).
# Reference values come from the original page writer save_coverage_tiles().

fake_covrle_files <- function(dir, runs, lens = c(chrA = 5500L, chrB = 2300L), seed = 1) {
  set.seed(seed)
  vapply(runs, function(run) {
    mk <- function() IRanges::RleList(lapply(lens, function(L) {
      S4Vectors::Rle(as.numeric(stats::rpois(L, 0.3) * sample(c(0, 1, 3), L, replace = TRUE)))
    }), compress = FALSE)
    path <- file.path(dir, paste0(run, ".covqs"))
    qs2::qs_save(ORFik::covRle(mk(), mk()), path)
    path
  }, character(1))
}

fake_seqinfo <- function(lens = c(chrA = 5500L, chrB = 2300L)) GenomeInfoDb::Seqinfo(names(lens), lens)

read_pages <- function(dir) {
  idx <- fst::read_fst(file.path(dir, "coverage_index.fst"), as.data.table = TRUE)
  files <- file.path(dir, basename(c(idx$file_forward, idx$file_reverse)))
  stats::setNames(lapply(files, fst::read_fst, as.data.table = TRUE), basename(files))
}

legacy_pages <- function(cov, runs, si, dir, tile) {
  dir.create(dir, showWarnings = FALSE)
  covs <- lapply(cov, ORFik:::read_RDSQS)
  fp <- file.path(dir, c("fwd.qs", "rev.qs"))
  fl <- lapply(covs, f); attr(fl, "strand") <- "+"; qs2::qs_save(fl, fp[1])
  rl <- lapply(covs, r); attr(rl, "strand") <- "-"; qs2::qs_save(rl, fp[2])
  pt <- megafst_page_table(si, tile, dir); chrs <- unique(pt$chr)
  for (p in fp) save_coverage_tiles(p, chrs, GenomeInfoDb::seqlengths(si)[chrs], dir, tile, runs,
                                    verbose = FALSE, remake = TRUE, BPPARAM = NULL)
  fst::write_fst(pt, file.path(dir, "coverage_index.fst"))
  read_pages(dir)
}

test_that("megafst_update full build is cell-identical to the original page writer", {
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:6)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo()
  ref <- legacy_pages(cov, runs, si, file.path(dir, "legacy"), tile = 2000)
  megafst_update(runs, cov, si, file.path(dir, "new"), tile_size = 2000, batch_size = 4,
                 BPPARAM = BiocParallel::SerialParam())
  expect_identical(read_pages(file.path(dir, "new")), ref)
  expect_identical(megafst_update(runs, cov, si, file.path(dir, "new"), tile_size = 2000,
                                  BPPARAM = BiocParallel::SerialParam())$plan$mode, "none")
})

test_that("incremental add, replace and removal reproduce a full build", {
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:6)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo(); sp <- BiocParallel::SerialParam()
  ref <- legacy_pages(cov, runs, si, file.path(dir, "legacy"), tile = 2000)
  out <- file.path(dir, "inc")
  megafst_update(runs[1:4], cov[1:4], si, out, tile_size = 2000, BPPARAM = sp)
  plan <- megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = sp)$plan
  expect_identical(plan$mode, "incremental"); expect_setequal(plan$add, runs[5:6])
  expect_identical(read_pages(out), ref)
  Sys.setFileTime(cov[2], Sys.time() + 10)
  plan <- megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = sp)$plan
  expect_identical(plan$replace, runs[2])
  expect_identical(read_pages(out), ref)
  writeLines("x", file.path(out, "keep.txt"))
  plan <- megafst_update(runs[-1], cov[-1], si, out, tile_size = 2000, BPPARAM = sp)$plan
  expect_identical(plan$mode, "full")
  expect_identical(read_pages(out), lapply(ref, function(x) x[, -1]))
  expect_true(file.exists(file.path(out, "keep.txt")))
  expect_false(dir.exists(paste0(out, ".staging")))
})

test_that("an interrupted append resumes to the same result", {
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:5)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo(); sp <- BiocParallel::SerialParam()
  ref <- legacy_pages(cov, runs, si, file.path(dir, "legacy"), tile = 2000)
  out <- file.path(dir, "res")
  megafst_update(runs[1:3], cov[1:3], si, out, tile_size = 2000, BPPARAM = sp)
  pages <- megafst_page_table(si, 2000, out)
  covs <- stats::setNames(lapply(cov[4:5], ORFik:::read_RDSQS), runs[4:5])
  megafst_write_batch(covs, pages[1:2], "add", 50, sp)
  expect_identical(megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = sp)$plan$mode, "incremental")
  expect_identical(read_pages(out), ref)
})

test_that("plan triggers compaction and full rebuilds where required", {
  lib <- data.table::data.table(Run = c("a", "b"), cov_path = c("x", "y"), cov_mtime = 1, cov_size = 1)
  p <- list(format_version = 1L, tile_size = 2e6, compress = 50, seq_names = "1", seq_lengths = 10,
            replaced_since_full = 0, appends_since_full = 64)
  new <- rbind(lib, data.table::data.table(Run = "c", cov_path = "z", cov_mtime = 1, cov_size = 1))
  expect_identical(megafst_plan(lib, new, p, p)$mode, "full")              # > 64 append passes
  p$appends_since_full <- 0
  expect_identical(megafst_plan(lib, new, p, p)$mode, "incremental")
  expect_identical(megafst_plan(lib, lib, p, p)$mode, "none")
  expect_identical(megafst_plan(NULL, lib, NULL, p)$mode, "full")
  expect_identical(megafst_plan(lib, lib, p, modifyList(p, list(tile_size = 1e6)))$mode, "full")
})

test_that("wait_for_memory proceeds, waits and errors", {
  expect_true(wait_for_memory(10, "x", 5, 0.01, 0.01, function() 100))
  expect_error(wait_for_memory(90, "x", 32, timeout_minutes = 0.0005, poll_seconds = 0.01,
                               available_fun = function() 100), "Not enough memory")
  calls <- 0
  avail <- function() { calls <<- calls + 1; if (calls < 3) 20 else 200 }
  expect_true(suppressMessages(wait_for_memory(50, "x", 10, timeout_minutes = 1, poll_seconds = 0.01,
                                               available_fun = avail)))
  expect_gte(calls, 3)
})

test_that("megacovrle_update caches by Run", {
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:4)
  cov <- fake_covrle_files(dir, runs); mc <- file.path(dir, "mc", c("f.qs", "r.qs"))
  expect_identical(megacovrle_update(runs[1:3], cov[1:3], mc), "full")
  expect_identical(megacovrle_update(runs, cov, mc), "incremental")
  expect_identical(megacovrle_update(runs, cov, mc), "none")
  fw <- qs2::qs_read(mc[1])
  expect_identical(names(fw), runs)
  expect_identical(fw[["SRR2"]], f(ORFik:::read_RDSQS(cov[2])))
  expect_identical(megacovrle_update(runs[-1], cov[-1], mc), "incremental")
  expect_identical(names(qs2::qs_read(mc[2])), runs[-1])
})

test_that("megacovrle_update adopts pre-caching unnamed files, names them on next load", {
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:3)
  cov <- fake_covrle_files(dir, runs); mc <- file.path(dir, "mc", c("f.qs", "r.qs")); dir.create(dirname(mc[1]))
  csv <- file.path(dir, "collection.csv"); writeLines("x", csv)
  covs <- lapply(cov, ORFik:::read_RDSQS)
  Sys.sleep(1.1)
  qs2::qs_save(lapply(covs, f), mc[1]); qs2::qs_save(lapply(covs, r), mc[2])   # legacy: unnamed
  expect_identical(megacovrle_update(runs, cov, mc, collection_csv = csv), "adopted")
  expect_identical(megacovrle_update(runs, cov, mc, collection_csv = csv), "none")
  cov4 <- fake_covrle_files(dir, "SRR4", seed = 9)
  expect_identical(megacovrle_update(c(runs, "SRR4"), c(cov, cov4), mc, collection_csv = csv), "incremental")
  fw <- qs2::qs_read(mc[1])
  expect_identical(names(fw), c(runs, "SRR4"))
  expect_identical(unname(fw[1:3]), unname(lapply(covs, f)))
})

test_that("pre-caching pages holding a subset of libraries are adopted and extended", {
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:6)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo(); sp <- BiocParallel::SerialParam()
  ref <- legacy_pages(cov, runs, si, file.path(dir, "full"), tile = 2000)
  old <- file.path(dir, "old")
  legacy_pages(cov[1:4], runs[1:4], si, old, tile = 2000)          # legacy pages, no manifest, 4 libraries
  dir.create(paste0(old, ".staging")); writeLines("x", file.path(paste0(old, ".staging"), "stale.fst"))
  plan <- megafst_update(runs, cov, si, old, tile_size = 2000, BPPARAM = sp)$plan
  expect_identical(plan$mode, "incremental"); expect_setequal(plan$add, runs[5:6])
  expect_identical(read_pages(old), ref)
  expect_false(dir.exists(paste0(old, ".staging")))
  exact <- file.path(dir, "exact"); legacy_pages(cov, runs, si, exact, tile = 2000)
  expect_identical(megafst_update(runs, cov, si, exact, tile_size = 2000, BPPARAM = sp)$plan$mode, "none")
  expect_true(file.exists(file.path(exact, "megafst_manifest.fst")))   # adoption persisted
  # a library present in the pages but removed from the collection cannot be adopted
  rm_dir <- file.path(dir, "rm"); legacy_pages(cov, runs, si, rm_dir, tile = 2000)
  expect_identical(megafst_update(runs[-2], cov[-2], si, rm_dir, tile_size = 2000, BPPARAM = sp)$plan$mode, "full")
})

test_that("C++ run-length writer (when available) gives pages identical to the R path", {
  skip_if_not(megafst_rle_available(), "fst without write_rle_columns_fst")
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:6)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo(); sp <- BiocParallel::SerialParam()
  withr::with_options(list(massiveNGSpipe.megafst_rle = FALSE),
    megafst_update(runs[1:4], cov[1:4], si, file.path(dir, "r"), tile_size = 2000, BPPARAM = sp))
  megafst_update(runs[1:4], cov[1:4], si, file.path(dir, "c"), tile_size = 2000, BPPARAM = sp)
  expect_identical(read_pages(file.path(dir, "c")), read_pages(file.path(dir, "r")))
  withr::with_options(list(massiveNGSpipe.megafst_rle = FALSE),
    megafst_update(runs, cov, si, file.path(dir, "r"), tile_size = 2000, BPPARAM = sp))   # R append
  megafst_update(runs, cov, si, file.path(dir, "c"), tile_size = 2000, BPPARAM = sp)      # C++ append
  expect_identical(read_pages(file.path(dir, "c")), read_pages(file.path(dir, "r")))
  Sys.setFileTime(cov[3], Sys.time() + 20)                                                 # one library changed
  withr::with_options(list(massiveNGSpipe.megafst_rle = FALSE),
    expect_identical(megafst_update(runs, cov, si, file.path(dir, "r"), tile_size = 2000, BPPARAM = sp)$plan$replace, runs[3]))
  expect_identical(megafst_update(runs, cov, si, file.path(dir, "c"), tile_size = 2000, BPPARAM = sp)$plan$replace, runs[3])
  expect_identical(read_pages(file.path(dir, "c")), read_pages(file.path(dir, "r")))      # C++ replace == R replace
})

test_that("stale .tmp files from an interrupted run are removed on the next update", {
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:3)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo(); sp <- BiocParallel::SerialParam()
  out <- file.path(dir, "m"); megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = sp)
  writeLines("half", file.path(out, "coverage_chrA_part1_forward.fst.tmp"))
  dir.create(paste0(out, ".staging")); writeLines("half", file.path(paste0(out, ".staging"), "x.fst.tmp"))
  megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = sp)
  expect_length(list.files(out, "\\.tmp$"), 0)
  expect_false(dir.exists(paste0(out, ".staging")))
  mc <- file.path(dir, "mc", c("f.qs", "r.qs")); megacovrle_update(runs, cov, mc)
  writeLines("half", paste0(mc[1], ".tmp"))
  megacovrle_update(runs, cov, mc)
  expect_false(file.exists(paste0(mc[1], ".tmp")))
})

test_that("forked page workers do not deadlock (no OpenMP thread changes inside workers)", {
  skip_on_os("windows")
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:4)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo()
  fst::threads_fst(2, reset_after_fork = FALSE); fst::write_fst(data.frame(a = 1:1e5), file.path(dir, "warm.fst"))
  res <- callr::r(function(runs, cov, si, out) {
    suppressMessages(devtools::load_all("~/livemount/forks/massiveNGSpipe", quiet = TRUE))
    fst::threads_fst(2, reset_after_fork = FALSE); fst::write_fst(data.frame(a = 1:1e5), tempfile())
    megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = BiocParallel::MulticoreParam(2))$plan$mode
  }, args = list(runs, cov, si, file.path(dir, "m")), timeout = 120)
  expect_identical(res, "full")
})

test_that("removed libraries are tombstoned (zeroed, kept) instead of forcing a full rebuild; re-adding overwrites", {
  skip_if_not(megafst_rle_available())
  dir <- withr::local_tempdir(); runs <- paste0("SRR", 1:10)
  cov <- fake_covrle_files(dir, runs); si <- fake_seqinfo(); sp <- BiocParallel::SerialParam()
  ref <- legacy_pages(cov, runs, si, file.path(dir, "legacy"), tile = 2000)
  out <- file.path(dir, "tomb")
  megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = sp)
  res <- megafst_update(runs[-3], cov[-3], si, out, tile_size = 2000, BPPARAM = sp)
  expect_identical(res$plan$mode, "incremental"); expect_identical(res$plan$removed, "SRR3")
  pages <- read_pages(out)
  for (nm in names(pages)) {
    expect_true(all(pages[[nm]]$SRR3 == 0L))                               # tombstone zeroed
    expect_identical(pages[[nm]][, -"SRR3"], ref[[nm]][, -"SRR3"])        # all other columns untouched
  }
  expect_identical(readRDS(file.path(out, "megafst_params.rds"))$tombstones, "SRR3")
  expect_false("SRR3" %in% fst::read_fst(file.path(out, "megafst_manifest.fst"))$Run)
  expect_identical(megafst_update(runs[-3], cov[-3], si, out, tile_size = 2000, BPPARAM = sp)$plan$mode, "none")
  # the library comes back: its existing column is overwritten (replace), pages equal a full build again
  res <- megafst_update(runs, cov, si, out, tile_size = 2000, BPPARAM = sp)
  expect_identical(res$plan$replace, "SRR3"); expect_length(res$plan$add, 0)
  expect_identical(read_pages(out), ref)
  expect_length(readRDS(file.path(out, "megafst_params.rds"))$tombstones, 0)
})

test_that("tombstones count toward the compaction threshold", {
  lib <- data.table::data.table(Run = paste0("R", 1:8), cov_path = paste0("p", 1:8), cov_mtime = 1, cov_size = 1)
  p <- list(format_version = 1, tile_size = 1, compress = 50, seq_names = "a", seq_lengths = 1,
            replaced_since_full = 0, appends_since_full = 0, tombstones = character())
  expect_identical(megafst_plan(lib, lib[-1], p, p)$mode, "incremental")            # 1 of 7 < 25%
  p$tombstones <- c("X1", "X2")
  expect_identical(megafst_plan(lib, lib[-1], p, p)$mode, "full")                   # 3 tombstones > 25% of 7
})
