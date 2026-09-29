# Unit tests for R/pipeline_names.R -- tiny pure string/list helpers.

test_that("get_experiment_names extracts conf['exp'] from every organism of every pipeline", {
  pipelines <- list(
    study1 = list(organisms = list(
      "Homo sapiens" = list(conf = c(exp = "study1-homo_sapiens")),
      "Mus musculus" = list(conf = c(exp = "study1-mus_musculus"))
    )),
    study2 = list(organisms = list(
      "Homo sapiens" = list(conf = c(exp = "study2-homo_sapiens"))
    ))
  )
  names <- unlist(get_experiment_names(pipelines), use.names = FALSE)
  expect_setequal(names, c("study1-homo_sapiens", "study1-mus_musculus", "study2-homo_sapiens"))
})

test_that("organism_merged_exp_name / organism_collection_exp_name replace spaces with underscores", {
  expect_identical(organism_merged_exp_name("Homo sapiens"), "all_merged-Homo_sapiens")
  expect_identical(organism_collection_exp_name("Homo sapiens"), "all_samples-Homo_sapiens")
  expect_identical(organism_merged_exp_name(c("Homo sapiens", "Mus musculus")),
                   c("all_merged-Homo_sapiens", "all_merged-Mus_musculus"))
})

test_that("organism_merged_modalities_exp_name appends _modalities to the merged name", {
  expect_identical(organism_merged_modalities_exp_name("Homo sapiens"),
                   "all_merged-Homo_sapiens_modalities")
})

test_that("reference_folder_name lowercases, trims, and underscores the organism name", {
  expect_identical(reference_folder_name("Homo sapiens"), "homo_sapiens")
  expect_identical(reference_folder_name("  Homo Sapiens  "), "homo_sapiens")
})

test_that("reference_folder_name: full=TRUE joins onto full_path without touching ORFik::config()", {
  # full_path is passed explicitly, so ORFik::config() (the default) is
  # never evaluated -- this must work with no ORFik config set up at all.
  expect_identical(reference_folder_name("Homo sapiens", full = TRUE, full_path = "/refs"),
                   file.path("/refs", "homo_sapiens"))
})
