test_that("metadata_perturbation_labels finds KO/KD/OE/mutant/degron and controls", {
  sym <- list(`Homo sapiens` = c("upf1", "mettl3", "eif2a", "rpl12", "pabpc1", "xrn1", "hts1", "zfp36", "zfp36l1", "rps16b"),
              `Saccharomyces cerevisiae` = c("dhh1", "rps16b", "ski2", "dom34"))
  tt <- c("HEK293 siUPF1 rep1", "HEK293 siCtrl rep1", "METTL3 KO ribo", "eIF2A S51A rep2", "HCT116_PABPC1_AID_plus_IAA",
          "XRN1-null replicate 2", "untreated cells (no HTS1 knock-down)", "HTS1 knock-down", "Zfp36/Zfp36l1 dKO 2",
          "dhh1Δ ribo YPD rep1", "rps16b-Null YPAD RNA-Seq rep1, sigma1278b deletion collection",
          "ski2Δdom34Δ ribo", "Empty vector control rep 3", "HEK293 cold-shock 2h", "HeLa siSNORD rep1", "WT rep1")
  org <- c(rep("Homo sapiens", 9), rep("Saccharomyces cerevisiae", 3), "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens")
  r <- metadata_perturbation_labels(tt, org, sym)
  expect_equal(r$CONDITION, c("KD", "control", "KO", "mutant", "degron", "KO", "control", "KD", "DKO", "KO", "KO", "DKO", "EV", "", "", ""))
  expect_equal(r$GENE[c(1, 3, 4, 5, 6, 8, 9, 10, 11, 12)],
               c("UPF1", "METTL3", "eIF2A_S51A", "PABPC1", "XRN1", "HTS1", "Zfp36_Zfp36l1", "dhh1", "rps16b", "ski2_dom34"))
})

test_that("metadata_fill_perturbations only fills empty fields and repairs old KO-rule genes", {
  dt <- data.table(sample_title = c("METTL3 KO ribo", "METTL3 KO ribo 2", "WT ribo", "MV4-11_RPL22KO_rep1", "ski2KO dom34KO short"),
                   ScientificName = "Homo sapiens", CONDITION = c("", "KO", "WT", "KO", "KO"),
                   GENE = c("", "", "", "MV4-11_RPL22", "ski2_dom34KO"))
  sym <- list(`Homo sapiens` = c("mettl3", "rpl22", "ski2", "dom34"))
  r <- metadata_fill_perturbations(dt, sym)
  expect_equal(r$dt$CONDITION, c("KO", "KO", "WT", "KO", "KO"))
  expect_equal(r$dt$GENE, c("METTL3", "METTL3", "", "RPL22", "ski2_dom34KO"))   # double-gene value kept
  expect_true(all(r$changes$old[r$changes$column == "CONDITION"] == ""))
})
