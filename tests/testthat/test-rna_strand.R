strand_test_reads <- function(n, opposite=FALSE) {
 x<-GenomicAlignments::GAlignments(seqnames=Rle(factor(rep('chr1',n))),pos=rep(11L,n),cigar=rep('5M20N5M',n),strand=Rle(factor(rep(if(opposite) '-' else '+',n),levels=c('+','-','*'))));x
}
test_that('Strand inference votes once per read, not once per exon', {
 regions<-GRanges('chr1',IRanges(c(11,36),width=5),strand='+')
 a<-massiveNGSpipe:::rna_strand_classify(strand_test_reads(1200),regions);b<-massiveNGSpipe:::rna_strand_classify(strand_test_reads(1200,TRUE),regions)
 expect_equal(a$strandMode,1L);expect_equal(b$strandMode,2L);expect_equal(a$informative_reads,1200)
})
test_that('Balanced and insufficient RNA evidence are distinct', {
 regions<-GRanges('chr1',IRanges(11,40),strand='+');x<-c(strand_test_reads(600),strand_test_reads(600,TRUE))
 expect_equal(massiveNGSpipe:::rna_strand_classify(x,regions)$status,'unstranded')
 expect_equal(massiveNGSpipe:::rna_strand_classify(x[1:10],regions)$status,'unknown')
})
test_that('Reads bridging conflicting strand regions do not vote', {
 regions<-GRanges('chr1',IRanges(c(11,36),width=5),strand=c('+','-'))
 expect_equal(massiveNGSpipe:::rna_strand_classify(strand_test_reads(1200),regions)$informative_reads,0)
})
test_that('RNA collapse refreshes OFST provenance and preserves strand', {
 file<-tempfile(fileext='.ofst');sidecar<-ORFik:::rna_strand_sidecar(file);on.exit(unlink(c(file,sidecar)))
 x<-strand_test_reads(1200,TRUE);ORFik::export.ofst(x,file)
 dir.create(dirname(sidecar),recursive=TRUE,showWarnings=FALSE)
 saveRDS(list(normalized=TRUE,read_strand='opposite',decision_id='test',product_signature=massiveNGSpipe:::rna_strand_signature(file)),sidecar)
 massiveNGSpipe:::rna_strand_collapse_ofst(ORFik::fimport(file),file)
 y<-ORFik::fimport(file);expect_equal(length(y),1);expect_equal(as.character(strand(y)),'-');expect_equal(mcols(y)$score,1200)
 expect_equal(ORFik:::rna_strand_record(file)$decision_id,'test')
})
test_that('Splice gaps are not evidence for exon orientation', {
 intron<-GRanges('chr1',IRanges(20,25),strand='+')
 expect_equal(massiveNGSpipe:::rna_strand_classify(strand_test_reads(1200),intron)$informative_reads,0)
})
