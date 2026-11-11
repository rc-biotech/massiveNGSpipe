# Organism masters (pipeline_merge_org) and collections (pipeline_collection_org) on a two-study fake organism built
# from ORFik's template genome/annotation. ORFik's config is mocked so experiments and tracks live in a temp folder.
skip_if_not_installed("rtracklayer")
skip_on_os("windows")

fake_organism <- function(env = parent.frame()) {
  root <- tempfile("mergeworld"); dir.create(root)
  withr::defer(unlink(root, recursive = TRUE), envir = env)
  bam <- file.path(root, "bam"); exp_dir <- file.path(root, "exp"); dir.create(bam); dir.create(exp_dir)
  testthat::local_mocked_bindings(config = function(...) c(exp = exp_dir, bam = bam, ref = root, fastq = root),
                                  .package = "ORFik", .env = env)
  txdb <- system.file("extdata/references/homo_sapiens", "Homo_sapiens_dummy.gtf.db", package = "ORFik")
  fa <- system.file("extdata/references/homo_sapiens", "Homo_sapiens_dummy.fasta", package = "ORFik")
  tpl <- function(x) system.file("extdata/Homo_sapiens_sample", paste0("RFP_", x, ".ofst"), package = "ORFik")
  studies <- c(PRJTESTA = "WT", PRJTESTB = "Mutant")
  for (s in names(studies)) {
    d <- file.path(bam, paste0(s, "-homo_sapiens"), "aligned", "pshifted")
    dir.create(file.path(d, "pshifted_merged"), recursive = TRUE)
    dir.create(file.path(d, "unique_mappers", "pshifted_merged"), recursive = TRUE)
    runs <- paste0("SRR", c(PRJTESTA = 1000, PRJTESTB = 2000)[[s]] + 1:2)
    file.copy(tpl(paste0(studies[[s]], "_rep", 1:2)), file.path(d, paste0(runs, ".ofst")))
    # study pools: all mappers = rep1 + rep2 rows, unique mappers = rep1 only (both point P-site tables)
    a <- rbind(fst::read_fst(tpl(paste0(studies[[s]], "_rep1"))), fst::read_fst(tpl(paste0(studies[[s]], "_rep2"))))
    fst::write_fst(a, file.path(d, "pshifted_merged", "RFP.ofst"))
    file.copy(tpl(paste0(studies[[s]], "_rep1")), file.path(d, "unique_mappers", "pshifted_merged", "RFP.ofst"))
    ORFik::create.experiment(d, exper = paste0(s, "-homo_sapiens"), txdb = txdb, fa = fa, organism = "Homo sapiens",
                             libtype = "RFP", rep = 1:2, saveDir = exp_dir)
    df <- ORFik::read.experiment(paste0(s, "-homo_sapiens"), validate = FALSE, output.env = new.env())
    suppressMessages(ORFik::countTable_regions(df, lib.type = "ofst", forceRemake = TRUE, BPPARAM = BiocParallel::SerialParam()))
  }
  config <- list(all_mappers = TRUE, split_unique_mappers = TRUE, preset = "Ribo-seq", mode = "online",
                 config = c(bam = bam, exp = exp_dir))
  list(root = root, bam = bam, exp = exp_dir, config = config, txdb = txdb,
       done = c("Homo sapiens" = "PRJTESTA-homo_sapiens", "Homo sapiens" = "PRJTESTB-homo_sapiens"))
}

pools <- function(w, mapper) {
  sub <- if (mapper == "all") "" else "unique_mappers"
  data.table::rbindlist(lapply(c("PRJTESTA", "PRJTESTB"), function(s)
    fst::read_fst(file.path(w$bam, paste0(s, "-homo_sapiens"), "aligned", "pshifted", sub, "pshifted_merged", "RFP.ofst"),
                  as.data.table = TRUE)), fill = TRUE)
}

master_files <- function(dir) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE)
  f <- f[!grepl("merge_input_files.csv$|removal_summary", f)]
  stats::setNames(unname(tools::md5sum(f)), sub(paste0("^", dir, "/"), "", f))
}

test_that("organism master: totals, per-position sums, covRle, bigWig and count tables are correct for both mappers", {
  w <- fake_organism()
  suppressMessages(pipeline_merge_org(w$config, pipelines = list(), done_experiments = w$done, suffix = "_test",
                                      parallel_mappers = TRUE))
  out <- file.path(w$bam, "all_merged-Homo_sapiens_test")
  inv <- data.table::fread(file.path(out, "merge_input_files.csv"))
  expect_equal(nrow(inv), 4L); expect_setequal(inv$mapper, c("all", "unique"))
  df <- ORFik::read.experiment("all_merged-Homo_sapiens_test", validate = FALSE, output.env = new.env())
  expect_equal(nrow(df), 1L)
  for (mapper in c("all", "unique")) {
    d <- df; if (mapper == "unique") ORFik::uniqueMappers(d) <- TRUE
    # read each branch's own files from its folder (filepath(d, "default") resolves to the all-mapper track)
    master <- file.path(ORFik::libFolder(d), "all.ofst")
    m <- fst::read_fst(master, as.data.table = TRUE)
    p <- pools(w, mapper)
    expect_equal(sum(as.numeric(m$score)), sum(as.numeric(p$score)))
    key <- c("seqnames", "start", "strand", "size")
    exp_sum <- p[, .(score = sum(score)), by = key][order(seqnames, start, strand, size)]
    got <- m[, .(score = sum(score)), by = key][order(seqnames, start, strand, size)]
    got[, (key[c(1, 3)]) := lapply(.SD, as.character), .SDcols = key[c(1, 3)]]
    exp_sum[, (key[c(1, 3)]) := lapply(.SD, as.character), .SDcols = key[c(1, 3)]]
    expect_equal(got, exp_sum, ignore_attr = TRUE)
    # covRle = summed P-sites per strand and position
    cov <- ORFik::fimport(file.path(ORFik::libFolder(d), "cov_RLE", "all.covqs"))
    expect_equal(sum(sum(ORFik::f(cov))) + sum(sum(ORFik::r(cov))), sum(as.numeric(p$score)))
    fwd <- p[strand == "+", .(s = sum(score)), by = .(seqnames = as.character(seqnames), start)]
    expect_true(all(mapply(function(ch, st, s) as.numeric(ORFik::f(cov)[[ch]][st]) == s, fwd$seqnames, fwd$start, fwd$s)))
    # bigWig = covRle, per strand
    bw <- list.files(file.path(ORFik::libFolder(d), "bigwig"), "\\.bigWig$", full.names = TRUE)
    expect_setequal(basename(bw), c("all_forward.bigWig", "all_reverse.bigWig"))
    imp <- rtracklayer::import.bw(grep("forward", bw, value = TRUE), as = "RleList")
    for (ch in names(ORFik::f(cov))) expect_equal(as.numeric(imp[[ch]]), as.numeric(ORFik::f(cov)[[ch]]))
    # count tables = direct weighted overlaps of the P-sites
    ct <- qs2::qs_read(file.path(ORFik::libFolder(d), "QC_STATS", "countTable_cds.qs"))
    gr <- ORFik::fimport(master)
    cds <- ORFik::loadRegion(w$txdb, "cds")
    expect_equal(unname(SummarizedExperiment::assay(ct)[, 1]), unname(ORFik::countOverlapsW(cds[rownames(ct)], gr, "score")))
  }
})

test_that("parallel and sequential merges give byte-identical masters", {
  w <- fake_organism()
  suppressMessages(pipeline_merge_org(w$config, pipelines = list(), done_experiments = w$done, suffix = "_par",
                                      parallel_mappers = TRUE))
  suppressMessages(pipeline_merge_org(w$config, pipelines = list(), done_experiments = w$done, suffix = "_seq",
                                      parallel_mappers = FALSE))
  a <- master_files(file.path(w$bam, "all_merged-Homo_sapiens_par"))
  b <- master_files(file.path(w$bam, "all_merged-Homo_sapiens_seq"))
  expect_setequal(names(a), names(b))
  expect_true(any(grepl("bigWig$", names(a)))); expect_true(any(grepl("unique_mappers/", names(a))))
  # qs count tables embed the experiment name in colData; compare their assays instead of bytes
  qs <- grepl("countTable_.*\\.qs$", names(a))
  expect_identical(a[!qs], b[names(a)[!qs]])
  for (f in names(a)[qs]) expect_identical(SummarizedExperiment::assay(qs2::qs_read(file.path(w$bam, "all_merged-Homo_sapiens_par", f))),
                                           SummarizedExperiment::assay(qs2::qs_read(file.path(w$bam, "all_merged-Homo_sapiens_seq", f))))
})

test_that("collection: all samples kept, count matrix is the column-bound study tables with Run and study keys", {
  w <- fake_organism()
  suppressMessages(pipeline_collection_org(w$config, pipelines = list(), done_experiments = w$done, path_suffix = "_test",
                                           BPPARAM = BiocParallel::SerialParam()))
  df <- ORFik::read.experiment("all_samples-Homo_sapiens_test", validate = FALSE, output.env = new.env())
  expect_equal(nrow(df), 4L)
  studies <- lapply(w$done, ORFik::read.experiment, validate = FALSE, output.env = new.env())
  expect_equal(ORFik::filepath(df, "default"), unlist(lapply(studies, ORFik::filepath, type = "default"), use.names = FALSE))
  for (region in c("mrna", "cds", "leaders", "trailers")) {
    ct <- qs2::qs_read(file.path(ORFik::QCfolder(df), paste0("countTable_", region, ".qs")))
    parts <- lapply(studies, function(s) SummarizedExperiment::assay(ORFik::countTable(s, region, type = "summarized")))
    expect_equal(unname(SummarizedExperiment::assay(ct)), unname(do.call(cbind, parts)))
    expect_equal(ct$Run, unlist(lapply(studies, ORFik::runIDs), use.names = FALSE))
    expect_equal(as.character(ct$study), rep(c("PRJTESTA", "PRJTESTB"), each = 2), info = paste(ct$study, collapse = ","))
    expect_equal(unname(readRDS(file.path(ORFik::QCfolder(df), paste0("totalCounts_", region, ".rds")))),
                 unname(colSums(SummarizedExperiment::assay(ct))))
  }
})

test_that("collection refuses study count tables whose labels do not match the experiment", {
  w <- fake_organism()
  s <- ORFik::read.experiment(w$done[[1]], validate = FALSE, output.env = new.env())
  f <- file.path(ORFik::QCfolder(s), "countTable_mrna.qs")
  ct <- qs2::qs_read(f); colnames(ct) <- rev(colnames(ct)); qs2::qs_save(ct, f)
  expect_error(suppressMessages(pipeline_collection_org(w$config, pipelines = list(), done_experiments = w$done,
                                                        path_suffix = "_test", BPPARAM = BiocParallel::SerialParam())),
               "do not match the experiment labels")
})

test_that("run_forked returns results in order, keeps NULL results and propagates errors", {
  r <- run_forked(list(a = function() 1, b = function() NULL, c = function() "x"), parallel = TRUE)
  expect_identical(r, list(a = 1, b = NULL, c = "x"))
  expect_identical(run_forked(list(a = function() 1, b = function() NULL), parallel = FALSE), list(a = 1, b = NULL))
  expect_error(run_forked(list(ok = function() 1, bad = function() stop("boom")), parallel = TRUE), "bad.*boom")
})
