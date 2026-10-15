test_that("pipeline_pshift raises a shifting error instead of returning silently (so it is recorded, not retried forever)", {
  flags <- character()
  local_mocked_bindings(step_is_next_not_done = function(...) TRUE,
                        pshift_needed_subset = function(...) simpleError("Library contained no periodic or accepted read-lengths"),
                        pshift_extra_lengths_hook = function(df, res, ...) res,
                        set_flag = function(config, step, e) flags <<- c(flags, e),
                        name = function(df) "E1")
  cfg <- list(all_mappers = TRUE, split_unique_mappers = TRUE, accepted_lengths_rpf = 26:30, reuse_shifts_if_existing = FALSE,
              pshift_extra_lengths = NULL, max_no_adapter_removed_pct = 80)
  expect_error(pipeline_pshift(list(list()), cfg, BPPARAM = BiocParallel::SerialParam()), "P-shifting failed: Library contained no periodic")
  expect_length(flags, 0)
})
