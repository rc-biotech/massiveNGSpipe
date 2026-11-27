# fastqc_adapters_info(): adapters_manual.csv rows are found by the file's run id, also when the file's folder
# holds only that run (paired footprint mates trimmed from <trim>/.mate_input/<Run>/, PRJNA951361 2026-11-27).

test_that("manual adapter is used for a mate_input link whose folder holds only that run", {
  trim <- withr::local_tempdir()
  data.table::fwrite(data.table::data.table(Run = c("SRR24042802", "SRR24042795", "SRR24042796"),
                                            adapter = c("AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG")),
                     file.path(trim, "adapters_manual.csv"))
  mate_dir <- file.path(trim, ".mate_input", "SRR24042795"); dir.create(mate_dir, recursive = TRUE)
  f <- file.path(mate_dir, "SRR24042795.fastq"); writeLines(c("@r1", "ACGT", "+", "IIII"), f)
  expect_identical(fastqc_adapters_info(f, adapter_manual_dir = trim), "CCCCCCCCCC")
  expect_identical(detect_adapter(f, adapter_manual_dir = trim), "CCCCCCCCCC")
})

test_that("manual adapter lookup by name handles mate suffixes, prefixes, .gz and prefix-sharing run ids", {
  runs <- c("SRR1", "SRR10", "SRR100")
  expect_identical(manual_adapter_row_for_file(runs, "/x/SRR10_1.fastq"), 2L)
  expect_identical(manual_adapter_row_for_file(runs, "/x/SRR10_R2_001.fastq.gz"), 2L)
  expect_identical(manual_adapter_row_for_file(runs, "/x/trimmed_SRR100.fastq"), 3L)
  expect_identical(manual_adapter_row_for_file(runs, "/x/collapsed_trimmed2_SRR1.fasta"), 1L)
  expect_identical(manual_adapter_row_for_file(runs, "/x/SRR1000.fastq"), integer())
  expect_identical(manual_adapter_row_for_file(runs, "/x/SRR1.fastq.tmp"), integer())
})

test_that("raw folder layout (all runs present, paired) still resolves the manual adapter", {
  d <- withr::local_tempdir()
  data.table::fwrite(data.table::data.table(Run = c("SRR7", "SRR8"), adapter = c("AAAAAAAAAA", "TTTTTTTTTT")),
                     file.path(d, "adapters_manual.csv"))
  for (r in c("SRR7", "SRR8")) for (m in 1:2) writeLines(c("@r", "ACGT", "+", "IIII"), file.path(d, sprintf("%s_%d.fastq", r, m)))
  expect_identical(fastqc_adapters_info(file.path(d, "SRR8_1.fastq")), "TTTTTTTTTT")
})
