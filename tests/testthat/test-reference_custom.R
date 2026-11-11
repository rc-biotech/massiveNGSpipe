test_that("add_custom_genome stages local and gzipped files, records provenance and refuses to overwrite", {
  ref <- tempfile("refs"); dir.create(ref); src <- tempfile("src"); dir.create(src)
  fa <- file.path(src, "g.fa"); Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chr1 = strrep("ACGT", 100))), fa)
  gtf <- file.path(src, "a.gtf"); writeLines('chr1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";', gtf)
  R.utils::gzip(gtf, paste0(gtf, ".gz"), remove = FALSE)
  seen <- NULL
  local_mocked_bindings(get_annotation = function(organism, output.dir, genome, GTF, db, gene_symbols, assembly) {
    seen <<- list(genome = genome, GTF = GTF, db = db, assembly = assembly)
    saveRDS(c(gtf = GTF, genome = genome), file.path(output.dir, "outputs.rds")); c(gtf = GTF, genome = genome)
  })
  ann <- add_custom_genome("Testus testus", genome = fa, annotation = paste0(gtf, ".gz"), ref_dir = ref, star_index = FALSE)
  d <- file.path(ref, reference_folder_name("Testus testus"))
  expect_true(file.exists(seen$genome)); expect_true(file.exists(seen$GTF)); expect_match(seen$GTF, "\\.gtf$")
  expect_identical(readLines(seen$GTF), readLines(gtf))                       # gunzipped copy
  expect_equal(seen$db, "genbank")
  src_rec <- readRDS(file.path(d, "reference_source.rds"))
  expect_equal(src_rec$organism, "Testus testus"); expect_equal(src_rec$genome_md5, unname(tools::md5sum(seen$genome)))
  expect_false(file.exists(file.path(d, "poor_assembly.rds")))
  expect_error(add_custom_genome("Testus testus", genome = fa, annotation = gtf, ref_dir = ref, star_index = FALSE), "exists already")
  expect_error(add_custom_genome("Testus testus", genome = fa, ref_dir = ref), "both genome and annotation")
})

test_that("ncbi_assembly_urls validates the accession and resolves the T2T goat files", {
  expect_error(ncbi_assembly_urls("GCA_123"), "Not an NCBI assembly accession")
  skip_on_cran(); skip_if_offline("ftp.ncbi.nlm.nih.gov")
  u <- ncbi_assembly_urls("GCA_041735815.1")
  expect_match(u[["genome"]], "GCA/041/735/815/GCA_041735815.1_T2T-goat2.0/GCA_041735815.1_T2T-goat2.0_genomic.fna.gz$")
  expect_match(u[["annotation"]], "_genomic.gtf.gz$")
})

test_that("add_custom_genome defaults to the assembly pinned for the organism", {
  expect_identical(massiveNGSpipe:::pinned_accession("Capra hircus"), "GCA_041735815.1")
  expect_identical(massiveNGSpipe:::pinned_accession("Halobacterium salinarum"), "GCF_000006805.1")
  expect_null(massiveNGSpipe:::pinned_accession("Homo sapiens"))
  expect_identical(eval(formals(add_custom_genome)$accession, list(organism = "Klebsiella michiganensis"),
                        asNamespace("massiveNGSpipe")), "GCF_002090195.1")
})
