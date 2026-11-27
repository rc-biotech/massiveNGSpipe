frame_test_collection <- function() {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  exp_dir <- file.path(root, "experiments")
  aligned <- file.path(root, "study-homo_sapiens", "aligned")
  coll_dir <- file.path(root, "all_samples-Homo_sapiens")
  dir.create(exp_dir)
  dir.create(coll_dir)
  dir.create(aligned, recursive = TRUE)
  testthat::local_mocked_bindings(config = function(...) c(exp = exp_dir),
                                  .package = "ORFik", .env = parent.frame())

  runs <- c("SRR90001", "SRR90002")
  paths <- file.path(aligned, paste0(runs, ".bam"))
  writeLines("placeholder", paths[1])
  writeLines("placeholder", paths[2])
  template <- ORFik::ORFik.template.experiment()
  txdb <- template@txdb
  fa <- template@fafile
  cds <- GenomicRanges::GRangesList(
    plus = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(11, 80), c(40, 100)), "+"),
    minus = GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 180), "-"))
  testthat::local_mocked_bindings(loadRegion = function(...) cds,
                                  .package = "ORFik", .env = parent.frame())

  make_experiment <- function(exp, condition, result) {
    ORFik::create.experiment(aligned, exper = exp, txdb = txdb, fa = fa,
      organism = "Homo sapiens", libtype = "RFP", stage = "test", condition = condition,
      rep = 1:2, fraction = "test", files = paths, saveDir = exp_dir, result_folder = result)
  }
  make_experiment("study-homo_sapiens", c("WT", "mutant"), aligned)
  make_experiment("all_samples-Homo_sapiens", c("irrelevantA", "irrelevantB"), coll_dir)
  branches <- c(all = aligned, unique = file.path(aligned, "unique_mappers"))
  for (branch in branches) {
    dir.create(file.path(branch, "ofst"), recursive = TRUE)
    dir.create(file.path(branch, "pshifted"))
  }
  branch <- branches[["all"]]
  raw <- file.path(branch, "ofst", paste0(runs, ".ofst"))
  # +: 5' ends at 10 (outside CDS) and 39 (frame 1); the second P-site
  # crosses a junction to 80 (frame 0). -: 5' end at 180 (frame 0).
  # The 4-nt read is outside the accepted lengths.
  reads <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(rep("chr1", 4)),
    pos = c(10L, 39L, 176L, 25L), cigar = c("5M", "2M39N3M", "5M", "4M"),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("+", "+", "-", "+"))))
  S4Vectors::mcols(reads)$score <- c(2, 5, 7, 100)
  shifts <- data.frame(fraction = 5L, offsets_start = -2L)
  for (i in seq_along(raw)) {
    ORFik::export.ofst(reads, raw[i])
    shifted <- suppressMessages(ORFik::shiftFootprints(reads, shifts))
    ORFik::export.ofst(shifted, file.path(branch, "pshifted", paste0(runs[i], "_pshifted.ofst")))
    S4Vectors::mcols(reads)$score <- S4Vectors::mcols(reads)$score * 2
  }
  table_path <- file.path(branch, "pshifted", "shifting_table.rds")
  offsets <- setNames(list(shifts, shifts), raw)
  saveRDS(offsets, table_path)
  unique_raw <- file.path(branches[["unique"]], "ofst", basename(raw[1]))
  file.copy(raw[1], unique_raw)
  file.copy(file.path(branch, "pshifted", paste0(runs[1], "_pshifted.ofst")),
    file.path(branches[["unique"]], "pshifted"))
  saveRDS(setNames(list(shifts), unique_raw), file.path(branches[["unique"]], "pshifted", "shifting_table.rds"))
  dir.create(file.path(aligned, "QC_STATS"))
  list(collection = "all_samples-Homo_sapiens", coll_dir = coll_dir, runs = runs,
       qc = file.path(aligned, "QC_STATS", "Ribo_frames_all.csv"),
       table_path = table_path, offsets = offsets, make_experiment = make_experiment,
       aligned = aligned, txdb = txdb, exp_file = file.path(exp_dir, "study-homo_sapiens.csv"))
}

test_that("without stored QC, P-sites are counted from the shifted ofst and 5' frames are rotated by the offset", {
  w <- frame_test_collection()
  data.table::fwrite(data.table::data.table(fraction = "stale_label", length = 5L,
                                           frame = 0L, score = 999), w$qc)
  expect_message(out <- collection_frames_update(w$collection),
                 "2 Run\\(s\\) had no usable stored frame QC")
  expect_identical(data.table::key(out), c("Run", "mapper", "read_length"))
  expect_identical(names(out), names(collection_frames_empty()))
  expect_identical(out$Run, w$runs[c(1, 1, 2)])
  expect_identical(out$mapper, c("all", "unique", "all"))
  expect_true(all(is.na(out$library_key)))
  expect_identical(out$experiment, rep("study-homo_sapiens", 3))
  expect_identical(out$read_length, rep(5L, 3))
  expect_identical(out$offset, rep(-2L, 3))
  # P-site frames 0/1/2 = 5/2/7; pre = frame (f - 2) mod 3. The score-2 read has its 5' end
  # outside the CDS (position 10), so rotation counts it in pre2 where the raw 5' end would not.
  expect_equal(as.numeric(unlist(out[1, c("frame0", "frame1", "frame2", "pre0", "pre1", "pre2")])),
               c(5, 2, 7, 7, 5, 2))
  expect_equal(out$frame0, c(5, 5, 10))
  path <- file.path(w$coll_dir, "QC_STATS", "frames_per_library_length.fst")
  expect_equal(fst::read_fst(path), as.data.frame(out))
})

test_that("rotation equals counting raw 5' ends for reads inside the CDS, including across a junction", {
  dir <- withr::local_tempdir()
  exons <- qc_cds_exon_table(GenomicRanges::GRangesList(
    plus = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(11, 80), c(40, 100)), "+"),
    minus = GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 180), "-")))
  reads <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(rep("chr1", 5)),
    pos = c(39L, 176L, 16L, 22L, 150L), cigar = c("2M39N3M", "5M", "5M", "6M", "6M"),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("+", "-", "+", "+", "-"))))
  S4Vectors::mcols(reads)$score <- c(5, 7, 3, 11, 13)
  offsets <- data.frame(fraction = 5:6, offsets_start = c(-2L, -4L))
  raw <- file.path(dir, "raw.ofst")
  shifted <- file.path(dir, "shifted.ofst")
  ORFik::export.ofst(reads, raw)
  ORFik::export.ofst(suppressMessages(ORFik::shiftFootprints(reads, offsets)), shifted)
  counts <- collection_frames_count(qc_frame_sums(qc_read_ofst(shifted), exons), offsets)
  # Previous implementation: frames of the raw 5' ends.
  r <- ORFik::import.ofst(raw)
  strand <- as.character(BiocGenerics::strand(r))
  pre <- qc_frame_sums(data.table::data.table(seqnames = as.character(GenomicRanges::seqnames(r)),
    start = ifelse(strand == "-", BiocGenerics::end(r), BiocGenerics::start(r)), strand = strand,
    size = ORFik::readWidths(r, along.reference = TRUE), score = S4Vectors::mcols(r)$score), exons)
  ref <- matrix(0, 2, 3)
  ref[cbind(pre$fraction - 4L, pre$frame + 1L)] <- pre$score
  expect_true(any(ref[, 1] > 0) && any(ref[, 2] > 0) && any(ref[, 3] > 0))
  expect_equal(as.matrix(counts[, .(pre0, pre1, pre2)]), ref, ignore_attr = TRUE)
  expect_equal(sum(counts[, .(pre0, pre1, pre2)]), sum(counts[, .(frame0, frame1, frame2)]))
})

test_that("stored QC is used before the ofst: per-sample cache by Run, then the study table by label", {
  w <- frame_test_collection()
  df <- ORFik::read.experiment(w$exp_file, validate = FALSE, output.env = new.env())
  labels <- sub("^RFP_", "", shift_qc_sample_names(df))
  sample_dir <- file.path(w$aligned, "QC_STATS", "per_sample_qc", "v2_complete_TIS_windows", "all")
  dir.create(sample_dir, recursive = TRUE)
  data.table::fwrite(data.table::data.table(fraction = labels[1], frame = 0:1, score = c(50, 60), length = 5L),
                     file.path(sample_dir, paste0(w$runs[1], "_frames.csv")))
  # libtype-prefixed labels are an accepted variant
  data.table::fwrite(data.table::data.table(fraction = paste0("RFP_", labels), length = 5L,
                                           frame = 2L, score = c(1, 30)), w$qc)
  read_ofst <- qc_read_ofst
  reads <- character()
  testthat::local_mocked_bindings(qc_read_ofst = function(path) {
    reads <<- c(reads, path)
    read_ofst(path)
  })
  expect_message(out <- collection_frames_update(w$collection), "1 Run\\(s\\) had no usable stored frame QC")
  expect_length(reads, 1L)
  expect_match(reads, "unique_mappers")
  expect_equal(as.numeric(unlist(out[1, c("frame0", "frame1", "frame2", "pre0", "pre1", "pre2")])),
               c(50, 60, 0, 0, 50, 60))
  expect_equal(as.numeric(unlist(out[3, c("frame0", "frame1", "frame2", "pre0", "pre1", "pre2")])),
               c(0, 0, 30, 30, 0, 0))
  expect_equal(out$frame0[2], 5)
})

test_that("read lengths without an offset are dropped with a warning, and a failing study is skipped", {
  counts <- collection_frames_count(data.table::data.table(fraction = c(5L, 7L), frame = c(1, 0), score = c(3, 9)),
                                    data.frame(fraction = 5L, offsets_start = -2L))
  expect_identical(attr(counts, "dropped_lengths"), 7L)
  expect_equal(as.numeric(unlist(counts[, .(frame0, frame1, frame2, pre0, pre1, pre2)])), c(0, 3, 0, 0, 0, 3))

  w <- frame_test_collection()
  sample_dir <- file.path(w$aligned, "QC_STATS", "per_sample_qc", "v2_complete_TIS_windows", "all")
  dir.create(sample_dir, recursive = TRUE)
  data.table::fwrite(data.table::data.table(fraction = "x", frame = 0:1, score = c(50, 60), length = c(5L, 9L)),
                     file.path(sample_dir, paste0(w$runs[1], "_frames.csv")))
  expect_warning(out <- collection_frames_update(w$collection), "1 read length\\(s\\) without a P-site offset dropped \\(9\\)")
  expect_equal(out$frame0[1], 50)
  expect_equal(out$frame1[1], 0)

  testthat::local_mocked_bindings(collection_frames_study = function(libs, ...) stop("broken study"))
  expect_warning(expect_warning(out <- collection_frames_update(w$collection, force = TRUE), "1 study\\(ies\\) skipped.*broken study"), "have no frame rows")
  expect_equal(nrow(out), 0L)
})

test_that("a fresh study report maps exact labels and rejects ambiguous or stale labels", {
  w <- frame_test_collection()
  df <- ORFik::read.experiment(w$exp_file, validate = FALSE, output.env = new.env())
  labels <- sub("^RFP_", "", shift_qc_sample_names(df))
  data.table::fwrite(data.table::data.table(fraction = rev(labels), length = 5L,
                                           frame = 0L, score = c(10, 5)), w$qc)
  mapped <- collection_frames_study_qc(w$qc, w$exp_file, w$txdb)
  expect_equal(mapped[[w$runs[1]]]$score, 5)
  expect_equal(mapped[[w$runs[2]]]$score, 10)
  Sys.setFileTime(w$qc, file.mtime(w$exp_file) - 1)
  expect_null(collection_frames_study_qc(w$qc, w$exp_file, w$txdb))
})

test_that("relabeling reuses the cache and changing one offset only invalidates its Run and mapper", {
  w <- frame_test_collection()
  count <- collection_frames_count
  calls <- 0L
  testthat::local_mocked_bindings(collection_frames_count = function(...) {
    calls <<- calls + 1L
    count(...)
  })
  out <- collection_frames_update(w$collection)
  expect_equal(calls, 3L)
  w$make_experiment(w$collection, c("newA", "newB"), w$coll_dir)
  expect_equal(collection_frames_update(w$collection), out)
  expect_equal(calls, 3L)

  w$offsets[[1]]$offsets_start <- -1L
  saveRDS(w$offsets, w$table_path)
  updated <- collection_frames_update(w$collection)
  expect_equal(calls, 4L)
  expect_identical(updated$offset, c(-1L, -2L, -2L))
})

test_that("the RFP collection driver runs frame QC by default and allows opting out", {
  w <- frame_test_collection()
  called <- character()
  testthat::local_mocked_bindings(collection_frames_update = function(collection, BPPARAM) {
    called <<- c(called, ORFik::name(collection))
  })
  bp <- BiocParallel::SerialParam()
  suppressMessages(fst_megaformat_make(organisms = w$collection, fst_threads = 1, BPPARAM = bp))
  expect_identical(called, w$collection)
  suppressMessages(fst_megaformat_make(organisms = w$collection, fst_threads = 1, BPPARAM = bp,
                                      make_frames = FALSE))
  expect_identical(called, w$collection)
})

test_that("a study report with a Run column maps by Run, even with stale labels", {
  w <- frame_test_collection()
  data.table::fwrite(data.table::data.table(fraction = c("stale_label_1", "stale_label_2"), length = 5L,
                                           frame = 0L, score = c(10, 5), Run = rev(w$runs)), w$qc)
  mapped <- collection_frames_study_qc(w$qc, w$exp_file, w$txdb)
  expect_equal(mapped[[w$runs[1]]]$score, 5)
  expect_equal(mapped[[w$runs[2]]]$score, 10)
})

test_that("a legacy study report where one label stands for two Runs is not mapped (no label collision)", {
  w <- frame_test_collection()
  df <- ORFik::read.experiment(w$exp_file, validate = FALSE, output.env = new.env())
  labels <- sub("^RFP_", "", shift_qc_sample_names(df))
  # a stale report that lists only one label twice cannot be attributed to both Runs
  data.table::fwrite(data.table::data.table(fraction = c(labels[1], labels[1]), length = 5L,
                                           frame = 0L, score = c(10, 5)), w$qc)
  expect_null(collection_frames_study_qc(w$qc, w$exp_file, w$txdb))
})
