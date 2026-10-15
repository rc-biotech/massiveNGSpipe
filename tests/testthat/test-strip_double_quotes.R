test_that("strip_double_quotes removes quote runs so fread/fwrite round trips are stable", {
  dt <- data.table(t = c("study \"GTEx\"", "a \"\"\"\" b", "plain"), n = 1:3)
  expect_warning(strip_double_quotes(dt), "3 values|2 values")
  expect_identical(dt$t, c("study 'GTEx'", "a ' b", "plain"))
  f <- tempfile(fileext = ".csv")
  for (i in 1:3) { fwrite(dt, f); dt <- fread(f, colClasses = c(t = "character")); suppressWarnings(strip_double_quotes(dt)) }
  expect_identical(dt$t, c("study 'GTEx'", "a ' b", "plain"))
})
