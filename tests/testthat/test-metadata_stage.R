# Staging helpers for new runs (R/metadata_stage.R)

test_that("metadata_make_unique numbers replicates, then adds auto_ placeholders", {
  dt <- data.table::data.table(study_accession = "P1", ScientificName = "Homo sapiens", Run = paste0("R", 1:5),
                               LIBRARYTYPE = "RFP", CONDITION = c("WT", "WT", "KO", "KO", "KO"),
                               REPLICATE = c("", "", "1", "1", ""))
  u <- metadata_make_unique(dt)
  expect_identical(u$REPLICATE[1:2], c("1", "2"))
  # KO: two rows share REPLICATE 1 -> auto placeholders; the KO row with empty REPLICATE is already unique
  expect_identical(u$FRACTION[3:4], c("auto_1", "auto_2"))
  expect_identical(c(u$REPLICATE[5], u$FRACTION[5]), c("", ""))
  expect_identical(u$FRACTION[1:2], c("", ""))
  # other studies / organisms are independent
  dt2 <- rbind(dt[1], data.table::copy(dt[1])[, `:=`(study_accession = "P2", Run = "R9")])
  expect_identical(metadata_make_unique(dt2)$REPLICATE, c("", ""))
})

test_that("metadata_biosample_duplicates finds re-deposits in another BioProject only", {
  new <- data.table::data.table(Run = c("N1", "N2", "N3"), BioProject = c("PB", "PB", "PB"),
                                BioSample = c("S1", "S2", ""))
  have <- data.table::data.table(Run = c("O1", "O2"), BioProject = c("PA", "PB"), BioSample = c("S1", "S2"))
  d <- metadata_biosample_duplicates(new, have)
  expect_identical(d$Run, "N1")
  expect_identical(d$other_bioprojects, "PA")
  expect_identical(nrow(metadata_biosample_duplicates(new[3], have)), 0L)
})

test_that("metadata_title_assay flags disome / SSU / RNA / IP titles", {
  expect_identical(metadata_title_assay(c("Disome 60mers, WT, rep1", "MG1655 cells, disome fraction", "40S TCP-seq",
                                          "input RNA rep1", "MC4100 lacZ, IP", "WT_Ribo_rep1", NA, "ribosome footprints chip-free")),
                   c("disome", "disome", "SSU", "RNA", "IP", "", "", ""))
  # footprint libraries whose titles mention input/total are not RNA
  expect_identical(metadata_title_assay(c("Ribosome profiling of input(total translatome) for NAC-IP replicate 1",
                                          "Bacterial cells, Total translatome, WT, rep1", "total RNA input")),
                   c("", "", "RNA"))
})

test_that("metadata_title_labels moves replicate numbers and keeps the varying title part as CONDITION", {
  dt <- data.table::data.table(study_accession = "P", ScientificName = "Danio rerio",
                               sample_title = c("w-ct4-1", "w-ct4-2", "bmal2-ct4-1", "bmal2-ct4-2"))
  u <- metadata_title_labels(dt)
  expect_identical(u$CONDITION, c("w", "w", "bmal2", "bmal2"))
  expect_identical(u$REPLICATE, c("1", "2", "1", "2"))
  dt2 <- data.table::data.table(study_accession = "Q", ScientificName = "Mesocricetus auratus",
                                sample_title = c("Wild-type golden hamster MII oocyte Rep 1", "Piwil1-/- golden hamster MII oocyte Rep 1"),
                                CONDITION = c("", ""))
  u2 <- metadata_title_labels(dt2)
  expect_identical(u2$REPLICATE, c("1", "1"))
  expect_false(u2$CONDITION[1] == u2$CONDITION[2])
  # existing CONDITION is kept
  dt3 <- data.table::copy(dt)[, CONDITION := c("WT", "WT", "KO", "KO")]
  expect_identical(metadata_title_labels(dt3)$CONDITION, c("WT", "WT", "KO", "KO"))
})

test_that("metadata_title_labels overrides a lossy CONDITION and ignores replicate tokens", {
  dt <- data.table::data.table(study_accession = "P", ScientificName = "Chlorocebus sabaeus",
                               sample_title = c("noCHX_Ribo_9hpi_KO2_2", "noCHX_Ribo_9hpi_mock_2", "noCHX_Ribo_9hpi_WT_2"),
                               CONDITION = c("", "WT", "WT"))
  expect_identical(metadata_title_labels(dt)$CONDITION, c("ko2", "mock", "wt"))
  dt2 <- data.table::data.table(study_accession = "G", ScientificName = "Gallus gallus",
                                sample_title = c("fowl.testis.rep4.x", "fowl.testis.rep5.x"), CONDITION = "")
  expect_identical(metadata_title_labels(dt2)$CONDITION, c("", ""))
})
