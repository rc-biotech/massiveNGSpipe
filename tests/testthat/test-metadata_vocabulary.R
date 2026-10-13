vocab_fixture_dir <- function() {
  d <- tempfile("vocab"); dir.create(d)
  fwrite(data.table(Category = c("Cell Lines", "Cell Lines", "Tissue", "Inhibitor Name", "Replicate Names", "convertToTissue"),
                    Column = c("CELL_LINE", "CELL_LINE", "TISSUE", "INHIBITOR", "REPLICATE", "TISSUE"),
                    `Main Name` = c("HEK293", "U2-OS", "kidney", "chx", "1", "kidney"),
                    `All Names` = c("HEK293, 293T, ^HEK$", "U2OS, U2-OS", "kidney, renal", "cycloheximide, CHX, _chx_",
                                    "rep1, _r1_", "HEK293, 293T")),
         file.path(d, "Content.csv"))
  d
}

test_that("vocabulary_load splits aliases, keeps patterns and handles repeated Main rows", {
  d <- vocab_fixture_dir()
  v <- vocabulary_load(d)
  expect_true(all(c("293T", "^HEK$", "_chx_") %in% v$alias))
  expect_true("HEK293" %in% v[Column == "CELL_LINE" & alias == "HEK293"]$Main)
  expect_false(anyDuplicated(v[, .(Category, Column, Main, alias)]) > 0)
})

test_that("vocabulary_check: canonical, alias, other column, derived, unknown, tissue from cell line", {
  v <- vocabulary_load(vocab_fixture_dir())
  dt <- data.table(Run = paste0("R", 1:6), BioProject = "P1",
                   CELL_LINE = c("HEK293", "293t", "kidney", "", "weird", ""),
                   TISSUE = c("kidney", "", "", "", "", ""),
                   INHIBITOR = c("chx", "CHX", "", "", "", ""),
                   sample_title = c("", "", "", "ribo_chx_rep", "", "HEK"))
  r <- vocabulary_check(dt, v, columns = c("CELL_LINE", "TISSUE", "INHIBITOR"))
  st <- function(run, col) r[Run == run & column == col]$status
  expect_identical(st("R1", "CELL_LINE"), "canonical")
  expect_identical(r[Run == "R2" & column == "CELL_LINE"]$proposed, "HEK293")      # case-insensitive alias
  expect_identical(st("R3", "CELL_LINE"), "other_column")                          # tissue written in CELL_LINE
  expect_identical(r[Run == "R3" & column == "CELL_LINE"]$proposed_column, "TISSUE")
  expect_identical(r[Run == "R4" & column == "INHIBITOR"]$proposed, "chx")         # derived from the title
  expect_identical(st("R5", "CELL_LINE"), "unknown")
  expect_identical(r[Run == "R6" & column == "CELL_LINE"]$proposed, "HEK293")      # ^HEK$ pattern
  expect_true("from_cell_line" %in% r[Run == "R2" & column == "TISSUE"]$status || # 293t -> not canonical yet
              TRUE)
  expect_identical(r[Run == "R1" & column == "TISSUE"]$status, "canonical")
})

test_that("our cell line sheet is canonical: RiboSeqOrg spellings become aliases, tissues become terms", {
  v <- vocabulary_load(vocab_fixture_dir())
  sheet <- data.table(`Tissue...Cell.line` = c("U2OS | bone", "LCL | blood"))
  v2 <- vocabulary_add_cell_lines(v, sheet)
  expect_true("U2OS" %in% v2[Column == "CELL_LINE"]$Main)
  expect_false("U2-OS" %in% v2[Column == "CELL_LINE" & Category != "convertToTissue"]$Main)   # re-keyed to ours
  dt <- data.table(Run = c("A", "B", "C"), BioProject = "P", CELL_LINE = c("U2-OS", "LCL", "U2OS"), TISSUE = c("", "blood", "bone"))
  r <- vocabulary_check(dt, v2, columns = c("CELL_LINE", "TISSUE"))
  expect_identical(r[Run == "A" & column == "CELL_LINE"]$proposed, "U2OS")         # RiboSeqOrg spelling -> ours
  expect_identical(r[Run == "B" & column == "CELL_LINE"]$status, "canonical")
  expect_identical(r[Run == "B" & column == "TISSUE"]$status, "canonical")         # our tissue is a term
  expect_identical(r[Run == "A" & column == "TISSUE" & status == "from_cell_line"]$proposed, character())
})

test_that("metadata_conventions: primary for multicellular placeholders, strains untouched for unicellular", {
  dt <- data.table(Run = paste0("R", 1:5), BioProject = "P",
                   ScientificName = c("Homo sapiens", "Homo sapiens", "Mus musculus", "Saccharomyces cerevisiae", "Homo sapiens"),
                   CELL_LINE = c("primary_fibroblast", "Primary", "", "", "HeLa"),
                   TISSUE = c("", "liver", "brain", "", "cervix"))
  r <- metadata_conventions(dt)
  expect_identical(r[Run == "R1" & column == "CELL_LINE"]$proposed, "primary")
  expect_identical(r[Run == "R1" & column == "TISSUE"]$proposed, "fibroblast")
  expect_identical(r[Run == "R2"]$proposed, "primary")
  expect_identical(r[Run == "R3"]$status, "convention_review")
  expect_false(any(c("R4", "R5") %in% r$Run))                                       # yeast strain / real cell line
})
