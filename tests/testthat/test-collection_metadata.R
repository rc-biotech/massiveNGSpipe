metadata_fixture_config <- function() {
  project <- tempfile("project"); dir.create(project)
  complete <- data.table(Run = paste0("R", 1:6), BioProject = c("P1", "P1", "P2", "P2", "P3", "P3"),
                         ScientificName = c(rep("Homo sapiens", 4), "Mus musculus", "Mus musculus"),
                         CELL_LINE = c("HeLa", "HeLa", "", "K562", "", ""), TISSUE = c("cervix", "cervix", "liver", "blood", "", ""),
                         sample_title = c("a", "b", "tumor x", "c", "d", "e"), sample_source = "")
  fwrite(complete, file.path(project, "complete.csv"))
  list(project = project, complete_metadata = file.path(project, "complete.csv"), mode = "online")
}

cell_line_fixture <- function() data.table(
  Tissue...Cell.line = c("HeLa | cervix", "K562 | blood"), Status..y.n. = c("y", "y"),
  Cancer.type = c("cervical", "leukemia"), Biological.class = c("cell line", "cell line"),
  Cell_type = c("epithelial", "lymphoblast"), Organ_system = c("reproductive", "blood"),
  Sex = c("female", "female"), Life_stage = c("adult", "adult"), CONDITION = "", TIMEPOINT = "")

test_that("extended metadata has exactly the done samples, in order, with cell line annotation", {
  config <- metadata_fixture_config()
  done <- fread(config$complete_metadata)[c(4, 1, 3, 5)]
  fwrite(done, file.path(config$project, "metadata_done_samples.csv"))
  expect_message(ext <- extend_done_samples_metadata(config, cell_lines = cell_line_fixture()),
                 "missing in the cell line sheet")
  expect_identical(ext$Run, c("R4", "R1", "R3", "R5"))
  expect_identical(ext$Cancer_type, c("leukemia", "cervical", NA, NA))
  expect_identical(ext$Cell_model[3], "cancer")          # "tumor" in sample title
  written <- fread(file.path(config$project, "metadata_done_samples_extended.csv"))
  expect_identical(written$Run, ext$Run)
  cancer <- fread(file.path(config$project, "complete_metadata_with_cancer_status.csv"))
  expect_setequal(cancer$Run, paste0("R", 1:4))
  # rewriting keeps a backup of the previous version
  Sys.sleep(1.1)
  extend_done_samples_metadata(config, cell_lines = cell_line_fixture())
  expect_length(list.files(file.path(config$project, "metadata_backups"), "extended"), 1)
  expect_false(file.exists(file.path(config$project, "metadata_done_samples_extended.csv.tmp")))
})

test_that("the cell line sheet is never written unless asked", {
  config <- metadata_fixture_config()
  fwrite(fread(config$complete_metadata), file.path(config$project, "metadata_done_samples.csv"))
  wrote <- FALSE
  local_mocked_bindings(write_sheet_safe = function(...) wrote <<- TRUE)
  suppressMessages(extend_done_samples_metadata(config, cell_lines = cell_line_fixture()))
  expect_false(wrote)
  suppressMessages(extend_done_samples_metadata(config, cell_lines = cell_line_fixture(), write_cell_line_sheet = TRUE))
  expect_true(wrote)
})

test_that("done samples metadata cannot silently shrink and must cover collections", {
  config <- metadata_fixture_config()
  studies <- list(P1 = c("R1", "R2"), P2 = c("R3", "R4"), P3 = c("R5", "R6"))
  collection_runs <- c("R1", "R3")
  fake_df <- function(e) structure(list(runs = if (e %in% names(studies)) studies[[e]] else collection_runs),
                                   class = "fake_exp")
  local_mocked_bindings(
    read.experiment = function(x, ...) fake_df(x),
    libraryTypes = function(df, ...) "RNA",
    list.experiments = function(...) data.table(name = c("all_samples-Homo_sapiens_RNA", "all_samples-Homo_sapiens_1_9_2025")),
    .package = "massiveNGSpipe")
  local_mocked_bindings(runIDs = function(df) df$runs, .package = "ORFik")
  run <- function(done, ...) suppressMessages(update_done_samples_metadata(
    config, pipelines = NULL, done_experiments = setNames(done, c("Homo sapiens", "Homo sapiens", "Mus musculus")[seq_along(done)]), ...))
  path <- file.path(config$project, "metadata_done_samples.csv")
  dt <- run(c("P1", "P2", "P3"), check_collections = FALSE)
  expect_identical(fread(path)$Run, paste0("R", 1:6))
  # a subset run must not shrink the file
  expect_error(run(c("P1", "P2"), check_collections = FALSE), "Refusing to remove 2 samples")
  expect_identical(fread(path)$Run, paste0("R", 1:6))
  run(c("P1", "P2"), check_collections = FALSE, force = TRUE)
  expect_identical(fread(path)$Run, paste0("R", 1:4))
  expect_true(length(list.files(file.path(config$project, "metadata_backups"))) >= 1)
  # current collection covered (dated snapshot ignored); a missing run stops
  expect_silent(suppressMessages(run(c("P1", "P2"))))
  collection_runs <- c("R1", "R7")
  expect_error(run(c("P1", "P2")), "collection runs lack done-samples metadata")
})
