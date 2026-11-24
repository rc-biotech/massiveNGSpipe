derive_available <- function() {
  !is.null(tryCatch(utils::getFromNamespace("ofst_derive_experiment", "ORFik"), error = function(e) NULL))
}

test_that("derive_master_formats: option off, several libraries, FALSE and errors all fall back", {
  skip_if_not(derive_available(), "ORFik without ofst_derive_experiment")
  calls <- 0L
  result <- FALSE
  local_mocked_bindings(ofst_derive_experiment = function(df, ...) {
    calls <<- calls + 1L
    if (identical(result, "error")) stop("boom")
    result
  }, .package = "ORFik")
  one <- data.frame(x = 1)

  withr::local_options(massiveNGSpipe.ofst_derive = FALSE)
  expect_false(derive_master_formats(one))
  expect_equal(calls, 0L)

  withr::local_options(massiveNGSpipe.ofst_derive = TRUE)
  expect_false(derive_master_formats(data.frame(x = 1:2)))   # group tracks: separate steps
  expect_equal(calls, 0L)
  expect_false(derive_master_formats(one))                   # not sorted point data
  expect_equal(calls, 1L)
  result <- "error"
  expect_warning(expect_false(derive_master_formats(one)), "using the separate steps")
  result <- invisible(list(counts = list()))
  expect_true(derive_master_formats(one))
})

test_that("make_additional_formats_internal skips the separate steps after a derive", {
  steps <- character()
  local_mocked_bindings(
    derive_master_formats = function(df) { steps <<- c(steps, "derive"); TRUE },
    make_point_covrle = function(df) { steps <<- c(steps, "covrle"); TRUE },
    master_bigwig = function(df, ...) { steps <<- c(steps, "bigwig"); invisible(NULL) },
    remove.experiments = function(df) invisible(NULL))
  make_additional_formats_internal(data.frame(), libtype_df_ref = TRUE, parallel = FALSE)
  expect_equal(steps, "derive")
})

# Real fixture: a sorted point OFST master built from ORFik's template RFP library. Old and new path must give the
# same files: covRle/covRleList/count tables byte-identical, bigWig files the same imported coverage (the ORFik
# writer's zoom summaries differ from rtracklayer's, full-resolution values do not).
test_that("ofst_derive path gives the same outputs as the separate steps (template master)", {
  skip_if_not(derive_available(), "ORFik without ofst_derive_experiment")
  tmpl <- ORFik::ORFik.template.experiment()
  rfp <- tmpl[tmpl$libtype == "RFP", ][1, ]
  x <- fst::read_fst(ORFik::filepath(rfp, "default"), as.data.table = TRUE)
  data.table::setorder(x, seqnames, start, strand, size)
  build <- function(derive) {
    d <- withr::local_tempdir(.local_envir = parent.frame(2))
    fst::write_fst(x, file.path(d, "all.ofst"))
    suppressMessages(ORFik::create.experiment(d, exper = "fx", txdb = tmpl@txdb, libtype = "RFP",
                                              fa = tmpl@fafile, organism = "Homo sapiens", saveDir = d))
    df <- ORFik::read.experiment(file.path(d, "fx.csv"), output.env = new.env())
    withr::with_options(list(massiveNGSpipe.ofst_derive = derive),
                        suppressMessages(make_additional_formats_internal(df, parallel = FALSE)))
    d
  }
  old <- build(FALSE); new <- build(TRUE)
  files <- setdiff(list.files(old, recursive = TRUE), "fx.csv")
  expect_setequal(setdiff(list.files(new, recursive = TRUE), "fx.csv"), files)
  expect_true(all(c("cov_RLE/all.covqs", "cov_RLE_List/all.covqs", "bigwig/all_forward.bigWig",
                    "QC_STATS/countTable_cds.qs") %in% files))
  for (f in grep("\\.bigWig$", files, value = TRUE, invert = TRUE))
    expect_equal(unname(tools::md5sum(file.path(new, f))), unname(tools::md5sum(file.path(old, f))), info = f)
  for (f in grep("\\.bigWig$", files, value = TRUE)) {
    a <- rtracklayer::import.bw(file.path(old, f), as = "RleList")
    b <- rtracklayer::import.bw(file.path(new, f), as = "RleList")
    expect_identical(b, a, info = f)
  }
})

test_that("derive_threads_dynamic respects a fixed option and stays within 1..16", {
  old <- options(ORFik.derive_threads = 3); on.exit(options(old))
  expect_identical(massiveNGSpipe:::derive_threads_dynamic(), 3L)
  options(ORFik.derive_threads = NULL)
  th <- massiveNGSpipe:::derive_threads_dynamic()
  expect_true(is.integer(th) && th >= 1L && th <= 16L)
  options(massiveNGSpipe.mem_reserve_gb = 1e9)  # no free memory above the reserve -> 1 thread
  expect_identical(massiveNGSpipe:::derive_threads_dynamic(), 1L)
  options(massiveNGSpipe.mem_reserve_gb = NULL)
})
