test_that("remove_flag_all_exp(samples = TRUE) also clears per-sample markers incl. _unique; default keeps them", {
  d <- tempfile("proj"); dir.create(d)
  cfg <- list(project = d)
  for (s in c("trim", "ofst_unique", "aligned")) for (r in c("R1", "R2")) set_sample_flag(cfg, s, "E1", r)
  set_sample_flag(cfg, "trim", "E2", "R9")
  local_mocked_bindings(is_config = function(x) TRUE, remove_flag = function(...) NULL)
  cfg$flag <- list(trim = 1, ofst = 2, aligned = 3)
  suppressMessages(remove_flag_all_exp(cfg, c("trim", "ofst"), "E1"))
  expect_true(dir.exists(sample_flag_dir(cfg, "trim", "E1")))
  suppressMessages(remove_flag_all_exp(cfg, c("trim", "ofst"), "E1", samples = TRUE))
  expect_false(dir.exists(sample_flag_dir(cfg, "trim", "E1"))); expect_false(dir.exists(sample_flag_dir(cfg, "ofst_unique", "E1")))
  expect_true(dir.exists(sample_flag_dir(cfg, "aligned", "E1")))          # step not reset
  expect_true(dir.exists(sample_flag_dir(cfg, "trim", "E2")))             # other experiment untouched
  suppressMessages(remove_flag_all_exp_from(cfg, "ofst", "E1", samples = TRUE))
  expect_false(dir.exists(sample_flag_dir(cfg, "aligned", "E1")))
})
