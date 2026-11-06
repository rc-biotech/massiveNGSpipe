test_that("ann_build_model places a GT..AG intron that keeps the frame", {
  # exon1: ATG + 8 codons (27 nt), intron with yeast motifs, exon2: 20 codons + TAA, flanks
  ex1 <- "ATGGCTAAAGAAGCTGCTCCAGCTAAG"
  intron <- paste0("GTATGT", strrep("A", 30), "TACTAAC", strrep("T", 12), "CAG")
  ex2 <- paste0(strrep("GCT", 20), "TAA")
  R <- Biostrings::DNAString(paste0(strrep("C", 60), ex1, intron, ex2, strrep("C", 60)))
  s1 <- 61L; e1 <- s1 + nchar(ex1) - 1L; s2 <- e1 + nchar(intron) + 1L; e2 <- s2 + 60L - 1L
  hs <- data.table::data.table(a = c(s1, s2), b = c(e1, e2), qstart = c(1L, 10L), qend = c(9L, 29L))
  m <- ann_build_model(R, hs, qlen = 30L)
  expect_false(is.null(m))
  expect_equal(start(m$exons), c(s1, s2))
  expect_equal(end(m$exons), c(e1, e2 + 3L))                     # last exon extended to the stop codon
  cds <- paste(vapply(seq_along(m$exons), function(i) as.character(Biostrings::subseq(R, start(m$exons)[i], end(m$exons)[i])), ""), collapse = "")
  expect_equal(nchar(cds) %% 3, 0)
  expect_false(grepl("*", sub("\\*$", "", as.character(Biostrings::translate(Biostrings::DNAString(cds)))), fixed = TRUE))
})

test_that("annotation_validate_cds and annotation_write_gtf round-trip a two-exon model", {
  g <- Biostrings::DNAStringSet(c(chr1 = paste0(strrep("C", 50), "ATGGCTAAA", "GTATGTAAAATACTAACTTTCAG", strrep("GCT", 10), "TAA", strrep("C", 50))))
  cds <- data.table::data.table(id = "g1", seqnames = "chr1", start = c(51L, 83L), end = c(59L, 115L), strand = "+",
                                gene_name = "G1", gene_biotype = "protein_coding", source_tag = "test")
  expect_true(annotation_validate_cds(cds, g)[["g1"]])
  f <- tempfile(fileext = ".gtf")
  annotation_write_gtf(cds, GenomicRanges::GRanges(), g, f, utr5 = 10L, utr3 = 10L)
  x <- rtracklayer::import(f)
  expect_equal(sort(unique(as.character(x$type))), sort(c("gene", "transcript", "exon", "CDS", "start_codon", "stop_codon")))
  expect_equal(range(start(x[x$type == "exon"])), c(41L, 83L))
  expect_equal(max(end(x[x$type == "CDS"])), 112L)                # CDS excludes the stop codon
  expect_equal(start(x[x$type == "stop_codon"]), 113L)
})
