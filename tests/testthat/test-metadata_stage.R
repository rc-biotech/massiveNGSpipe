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

test_that("metadata_biosample_labels fills only empty fields, ignores junk and canonicalises", {
  dt <- data.table::data.table(BioSample = c("SAM1", "SAM2", "SAM3"), CELL_LINE = c("", "HeLa", ""), TISSUE = "")
  attrs <- data.table::data.table(BioSample = c("SAM1", "SAM1", "SAM2", "SAM3"), attribute = c("cell_line", "tissue", "cell_line", "cell_line"),
                                  value = c("HEK-293", "kidney", "HepG2", "not applicable"))
  vocab <- data.table::data.table(Category = "Cell Lines", Column = "CELL_LINE", Main = "HEK293", alias = c("HEK293", "HEK-293"))
  u <- metadata_biosample_labels(dt, attrs, vocab)
  expect_identical(u$CELL_LINE, c("HEK293", "HeLa", ""))
  expect_identical(u$TISSUE, c("kidney", "", ""))
  expect_match(u$COMMENT[1], "author-declared")
  expect_identical(u$COMMENT[2:3], c("", ""))
})

test_that("metadata_biosample_attributes reads the cache without network when all records are cached", {
  d <- withr::local_tempdir()
  saveRDS(data.table::data.table(BioSample = "SAMX1", attribute = "tissue", value = "leaf"), file.path(d, "SAMX1.rds"))
  a <- metadata_biosample_attributes("SAMX1", d)
  expect_identical(a$value, "leaf")
  expect_identical(nrow(metadata_biosample_attributes(character(), d)), 0L)
})

test_that("metadata_cell_line_tokens extends truncated names and clears foreign lines", {
  dt <- data.table::data.table(ScientificName = "Homo sapiens", CELL_LINE = c("S2", "MDA", "U2OS", "S2", "HeLa"),
                               sample_title = c("AK_Riboseq_hESC_Rep1", "MDA-MB-231 WT, DMSO", "U2OS WT", "S2_1_1.fq.gz", "HeLa-S3 cells"),
                               sample_source = c("RUES2", "Breast", "Osteosarcoma", "", ""))
  u <- metadata_cell_line_tokens(dt)
  # S2 on a human sample is cleared before any token logic (Drosophila line), RUES2 then comes from BioSample/text later
  expect_identical(u$CELL_LINE, c("", "MDA-MB-231", "U2OS", "", "HeLa-S3"))
  expect_match(u$COMMENT[2], "extended")
  expect_match(u$COMMENT[1], "other species")
  # non-human species keep S2
  expect_identical(metadata_cell_line_tokens(data.table::data.table(ScientificName = "Drosophila melanogaster",
                   CELL_LINE = "S2", sample_title = "S2 cells"))$CELL_LINE, "S2")
})
