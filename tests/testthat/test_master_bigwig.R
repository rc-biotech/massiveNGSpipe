test_that("master_bigwig sends multi-library experiments (group tracks) to convert_to_bigWig", {
  called <- NULL
  local_mocked_bindings(
    filepath = function(df, type, ...) c("/x/RFP_elongation.covqs", "/x/RFP_initiation.covqs"),
    convert_to_bigWig = function(df, in_files, ...) { called <<- in_files; invisible(NULL) },
    fimport = function(...) stop("fimport must not be called with several files"))
  expect_no_error(master_bigwig(data.frame(), parallel = FALSE))
  expect_equal(called, c("/x/RFP_elongation.covqs", "/x/RFP_initiation.covqs"))
})
