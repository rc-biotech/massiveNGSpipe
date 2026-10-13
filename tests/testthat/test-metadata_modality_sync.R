modality_fixture <- function() {
  base <- tempfile("modality"); dir.create(base)
  mk <- function(project, rows) {
    dir.create(file.path(base, project))
    fwrite(rows, file.path(base, project, "next_round_manual.csv"))
  }
  row <- function(run, lt, cond = "", rep = 1L) data.table(Run = run, BioProject = "P1", LIBRARYTYPE = lt,
                                                          CONDITION = cond, REPLICATE = rep, KEEP = TRUE)
  mk("NGS_pipeline", rbindlist(list(row("R1", "RFP"), row("R2", "RNA"), row("R3", "RNA"), row("R4", "RNA", "x"),
                                    row("R5", "SLAM"), row("R6", "disome"), row("R7", "SSU"))))
  mk("NGS_pipeline_RNA", rbindlist(list(row("R3", "RNA"), row("R4", "RNA", "y"), row("R9", "RNA"))))
  mk("NGS_pipeline_disome", rbindlist(list(row("D1", "disome"), row("R1b", "RFP"))))
  mk("NGS_pipeline_SSU", rbindlist(list(row("R7", "RNA"))))          # same run, other modality label
  base
}

test_that("plan: move, duplicate drop, conflicts and parking", {
  base <- modality_fixture()
  plan <- suppressMessages(sync_modalities("online", base_dir = base, use_google = FALSE, dry_run = TRUE))
  a <- setNames(plan$action, paste(plan$Run, plan$from))
  expect_identical(a[["R2 NGS_pipeline"]], "move")             # RNA not in RNA sheet
  expect_identical(a[["R3 NGS_pipeline"]], "drop_duplicate")   # same annotation already there
  expect_identical(a[["R4 NGS_pipeline"]], "conflict")         # different CONDITION
  expect_identical(a[["R5 NGS_pipeline"]], "park")             # no SLAM pipeline
  expect_identical(a[["R6 NGS_pipeline"]], "move")
  expect_identical(a[["R1b NGS_pipeline_disome"]], "move")     # RFP row in the disome sheet
  expect_identical(a[["R7 NGS_pipeline"]], "conflict")         # labelled SSU here, RNA in the SSU sheet
  expect_false("R1 NGS_pipeline" %in% names(a))                # own modality stays
  # dry run writes nothing
  expect_identical(nrow(fread(file.path(base, "NGS_pipeline", "next_round_manual.csv"))), 7L)
})

test_that("apply moves rows once, removes only confirmed rows, keeps conflicts", {
  base <- modality_fixture()
  suppressMessages(sync_modalities("online", base_dir = base, use_google = FALSE, dry_run = FALSE))
  rfp <- fread(file.path(base, "NGS_pipeline", "next_round_manual.csv"))
  rna <- fread(file.path(base, "NGS_pipeline_RNA", "next_round_manual.csv"))
  dis <- fread(file.path(base, "NGS_pipeline_disome", "next_round_manual.csv"))
  oth <- fread(file.path(base, "NGS_pipeline_other", "next_round_manual.csv"))
  expect_setequal(rfp$Run, c("R1", "R4", "R7", "R1b"))         # conflicts stay, R1b arrived
  expect_setequal(rna$Run, c("R3", "R4", "R9", "R2"))          # R2 added once, R4 untouched
  expect_identical(rna[Run == "R4"]$CONDITION, "y")
  expect_setequal(dis$Run, c("D1", "R6"))
  expect_identical(oth$Run, "R5")
  expect_identical(oth$LIBRARYTYPE, "SLAM")
  # idempotent: second run has nothing left to move
  plan2 <- suppressMessages(sync_modalities("online", base_dir = base, use_google = FALSE, dry_run = TRUE))
  expect_true(all(plan2$action == "conflict"))
})

test_that("on_conflict = archive: home sheet wins, stray copy archived, homeless runs kept", {
  base <- modality_fixture()
  # R8: labelled CAGE here and RNA only in the disome sheet (no home-sheet copy) -> must be kept
  fwrite(rbind(fread(file.path(base, "NGS_pipeline", "next_round_manual.csv")),
               data.table(Run = "R8", BioProject = "P1", LIBRARYTYPE = "CAGE", CONDITION = "", REPLICATE = 1L, KEEP = TRUE)),
         file.path(base, "NGS_pipeline", "next_round_manual.csv"))
  fwrite(rbind(fread(file.path(base, "NGS_pipeline_disome", "next_round_manual.csv")),
               data.table(Run = "R8", BioProject = "P1", LIBRARYTYPE = "RNA", CONDITION = "", REPLICATE = 1L, KEEP = TRUE)),
         file.path(base, "NGS_pipeline_disome", "next_round_manual.csv"))
  plan <- suppressMessages(sync_modalities("online", base_dir = base, use_google = FALSE, on_conflict = "archive", dry_run = FALSE))
  rfp <- fread(file.path(base, "NGS_pipeline", "next_round_manual.csv"))
  arch <- fread(file.path(base, "NGS_pipeline_other", "conflict_archive.csv"))
  rna <- fread(file.path(base, "NGS_pipeline_RNA", "next_round_manual.csv"))
  expect_identical(arch$Run, "R4")                             # annotation conflict: home sheet wins
  expect_identical(arch[Run == "R4"]$CONDITION, "x")           # stray annotation preserved in the archive
  expect_identical(rna[Run == "R4"]$CONDITION, "y")            # home sheet unchanged
  expect_false("R4" %in% rfp$Run)
  # R7 (SSU here, RNA-labelled in the SSU sheet) and R8 (CAGE here, RNA in the disome sheet) have no
  # copy in their label's home sheet: kept, not archived
  expect_true(all(c("R7", "R8") %in% rfp$Run))
  expect_identical(plan[Run == "R8" & from == "NGS_pipeline"]$action, "conflict")
})

test_that("online and local families never mix", {
  expect_true(all(grepl("^NGS_pipeline", modality_homes("online"))))
  expect_true(all(grepl("^local_pipeline", modality_homes("local"))))
  expect_identical(modality_other_project("local"), "local_pipeline_other")
})
