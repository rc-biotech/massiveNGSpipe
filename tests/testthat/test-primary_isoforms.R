gtf_line <- function(feat, start, end, gene, tx, tags = character()) {
  tg <- paste0(vapply(tags, function(t) paste0(' tag "', t, '";'), ""), collapse = "")
  paste("1", "test", feat, start, end, ".", "+", ".",
        paste0('gene_id "', gene, '"; transcript_id "', tx, '";', tg), sep = "\t")
}
gtf_tx <- function(gene, tx, exons, cds = NULL, tags = character()) {
  st <- min(unlist(exons)); en <- max(unlist(exons))
  c(gtf_line("transcript", st, en, gene, tx, tags),
    vapply(exons, function(e) gtf_line("exon", e[1], e[2], gene, tx, tags), ""),
    if (!is.null(cds)) gtf_line("CDS", cds[1], cds[2], gene, tx, tags))
}
make_test_gtf <- function() {
  f <- tempfile(fileext = ".gtf")
  writeLines(c(
    # gene A: t1 longest CDS, t2 MANE (shorter CDS), t3 non-coding
    gtf_tx("A", "t1", list(c(100, 400)), c(110, 390)),
    gtf_tx("A", "t2", list(c(100, 300)), c(110, 250), tags = c("basic", "Ensembl_canonical", "MANE_Select")),
    gtf_tx("A", "t3", list(c(100, 800))),
    # gene B: t4 Ensembl canonical, t5 longer CDS
    gtf_tx("B", "t4", list(c(1000, 1200)), c(1010, 1150), tags = "Ensembl_canonical"),
    gtf_tx("B", "t5", list(c(1000, 1300)), c(1010, 1290)),
    # gene C: non-coding, t7 longer
    gtf_tx("C", "t6", list(c(2000, 2100))),
    gtf_tx("C", "t7", list(c(2000, 2300))),
    # gene D: same CDS, t9 longer transcript
    gtf_tx("D", "t8", list(c(3000, 3200)), c(3010, 3150)),
    gtf_tx("D", "t9", list(c(3000, 3400)), c(3010, 3150))), f)
  f
}

test_that("primary_isoforms auto: MANE > canonical > longest CDS / longest tx, one row per gene", {
  skip_if_not_installed("txdbmaker")
  f <- make_test_gtf()
  p <- suppressWarnings(primary_isoforms(f))
  expect_equal(p$gene_id, c("A", "B", "C", "D"))
  expect_equal(p$tx_name, c("t2", "t4", "t7", "t9"))
  expect_equal(p$method, c("mane", "ensembl_canonical", "longest_tx", "longest_cds"))
  expect_equal(p$coding, c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(p$n_tx, c(3L, 2L, 2L, 2L))
})

test_that("primary_isoforms single methods, fallback and expression", {
  skip_if_not_installed("txdbmaker")
  f <- make_test_gtf()
  lc <- suppressWarnings(primary_isoforms(f, "longest_cds"))
  expect_equal(lc$tx_name, c("t1", "t5", "t7", "t9"))
  m <- suppressWarnings(primary_isoforms(f, "mane"))
  expect_equal(m$method, c("mane", "longest_cds", "longest_tx", "longest_cds"))
  m0 <- suppressWarnings(primary_isoforms(f, "mane", fallback = FALSE))
  expect_equal(m0$tx_name, "t2")
  ex <- c(t8 = 50, t9 = 1, t5 = 100, t3 = 1000)
  a <- suppressWarnings(primary_isoforms(f, expression = ex))
  expect_equal(a$tx_name, c("t2", "t4", "t7", "t8"))      # tags win over expression; t3 is non-coding in a coding gene
  expect_equal(a$method[4], "expression")
  e <- suppressWarnings(primary_isoforms(f, "expression", expression = data.frame(tx_name = "t5", value = 3)))
  expect_equal(e$tx_name[2], "t5")
  cache <- tempfile(fileext = ".csv")
  p1 <- suppressWarnings(primary_isoforms(f, cache_file = cache))
  expect_true(file.exists(cache))
  expect_equal(primary_isoforms("does_not_matter", cache_file = cache)$tx_name, p1$tx_name)
})

test_that("primary_isoform_tags reads NCBI GFF3 RefSeq Select from transcript and child lines", {
  f <- tempfile(fileext = ".gff")
  writeLines(c("##gff-version 3",
    "1\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=gene-X;Name=X",
    "1\tRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=rna-NM_1.1;Parent=gene-X;Name=NM_1.1;tag=RefSeq Select",
    "1\tRefSeq\texon\t1\t500\t.\t+\t.\tID=exon-NM_2.1-1;Parent=rna-NM_2.1;tag=MANE Select",
    "1\tRefSeq\tmRNA\t1\t400\t.\t+\t.\tID=rna-NM_3.1;Parent=gene-X;Name=NM_3.1"), f)
  t <- primary_isoform_tags(f)
  expect_true(t[id == "NM_1.1", refseq_select])
  expect_true(t[id == "NM_2.1", mane])
  expect_false("NM_3.1" %in% t$id)
  expect_false("gene-X" %in% t$id)
})

test_that("subset_primary_isoforms keeps table and ranges aligned", {
  tab <- data.table::data.table(ensembl_tx_name = c("t1", "t2", "t2", "t9"), type = c("uORF", "uORF", "dORF", "NTE"))
  grl <- GenomicRanges::GRangesList(lapply(1:4, function(i) GenomicRanges::GRanges("1", IRanges::IRanges(i, i + 2), "+")))
  s <- subset_primary_isoforms(tab, data.table::data.table(tx_name = c("t2", "t9")), grl)
  expect_equal(s$table$type, c("uORF", "dORF", "NTE"))
  expect_equal(as.integer(BiocGenerics::start(unlist(s$grl))), 2:4)
  expect_error(subset_primary_isoforms(tab, "t2", grl[1:2]))
})
