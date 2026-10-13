test_that("enabling contam on an existing project: legacy aligned studies skip it, new studies require it", {
  d <- tempfile("flags")
  steps <- c("collapsed", "contam", "aligned", "cleanbam")
  flag <- setNames(file.path(d, steps), steps); for (f in flag) dir.create(f, recursive = TRUE)
  config <- list(flag = flag)
  for (s in c("collapsed", "aligned")) saveRDS(TRUE, file.path(flag[s], "old.rds"))   # aligned before contam existed
  saveRDS(TRUE, file.path(flag["collapsed"], "new.rds"))
  expect_false(step_is_next_not_done(config, "contam", "old"))
  expect_true(step_is_next_not_done(config, "cleanbam", "old"))
  expect_true(step_is_next_not_done(config, "contam", "new"))
  expect_false(step_is_next_not_done(config, "aligned", "new"))                        # needs contam first
  saveRDS(TRUE, file.path(flag["contam"], "new.rds"))
  expect_true(step_is_next_not_done(config, "aligned", "new"))
})
