test_that("make_point_covrle builds every library of a multi-library (group) experiment", {
  skip_if(is.null(tryCatch(utils::getFromNamespace(".convert_point_ofst_cov", "ORFik"), error = function(e) NULL)),
          "ORFik without the point coverage builder")
  out <- withr::local_tempdir()
  calls <- character()
  local_mocked_bindings(
    filepath = function(df, type, ...) c("/x/pshifted/RFP_g1.ofst", "/x/pshifted/RFP_g2.ofst"),
    libFolder = function(df, ...) c(out, out),
    seqinfo = function(x) NULL)
  local_mocked_bindings(.convert_point_ofst_cov = function(in_file, out_list, out_merged, seqinfo, format) {
    calls <<- c(calls, basename(out_merged)); TRUE }, .package = "ORFik")
  expect_true(make_point_covrle(data.frame()))
  expect_equal(calls, c("RFP_g1", "RFP_g2"))
  expect_true(dir.exists(file.path(out, "cov_RLE")))
})

test_that("make_point_covrle reports FALSE when any library fails (caller falls back)", {
  skip_if(is.null(tryCatch(utils::getFromNamespace(".convert_point_ofst_cov", "ORFik"), error = function(e) NULL)),
          "ORFik without the point coverage builder")
  out <- withr::local_tempdir()
  local_mocked_bindings(
    filepath = function(df, type, ...) c("/x/RFP_g1.ofst", "/x/RFP_g2.ofst"),
    libFolder = function(df, ...) c(out, out),
    seqinfo = function(x) NULL)
  local_mocked_bindings(.convert_point_ofst_cov = function(in_file, ...) grepl("g1", in_file), .package = "ORFik")
  expect_false(make_point_covrle(data.frame()))
})
