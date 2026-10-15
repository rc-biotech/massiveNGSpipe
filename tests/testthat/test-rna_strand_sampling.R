test_that('Whole-BAM sampling defeats a misleading coordinate prefix without changing thresholds', {
 skip_if(Sys.which('samtools') == '', 'samtools required')
 d <- tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE))
 sam <- file.path(d,'sorted.sam');bam <- file.path(d,'sorted.bam')
 n <- 100000L
 header <- c('@HD\tVN:1.6\tSO:coordinate','@SQ\tSN:chr1\tLN:10000','@SQ\tSN:chr2\tLN:10000')
 first <- paste0('a',seq_len(n),'\t16\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNH:i:1')
 second <- paste0('b',seq_len(n),'\t0\tchr2\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNH:i:1')
 writeLines(c(header,first,second),sam)
 stopifnot(system2(Sys.which('samtools'),c('view','-b','-o',shQuote(bam),shQuote(sam)))==0L)
 regions <- GenomicRanges::GRanges(c('chr1','chr2'),IRanges::IRanges(90,120),strand='+')
 annotation <- list(regions=regions,annotation_signature='fixture')
 prefix <- rna_strand_detect(bam,annotation)
 expect_identical(prefix$strandMode,2L)
 expect_identical(prefix$scanned_reads,n)
 full <- rna_strand_detect(bam,annotation,sampling='whole_bam_random',sampling_fraction=1)
 expect_identical(full$strandMode,0L)
 expect_equal(full$scanned_reads,2*n)
 expect_equal(full$sense_fraction,.5)
 expect_identical(full$source_bam,normalizePath(bam))
 expect_identical(full$product_signature,rna_strand_signature(bam))
 sampled <- rna_strand_detect(bam,annotation,sampling='whole_bam_random',sampling_fraction=.1)
 repeated <- rna_strand_detect(bam,annotation,sampling='whole_bam_random',sampling_fraction=.1)
 expect_identical(sampled$strandMode,0L)
 expect_identical(sampled$parameters,repeated$parameters)
 expect_identical(sampled$sense_reads,repeated$sense_reads)
 expect_identical(sampled$antisense_reads,repeated$antisense_reads)
 expect_error(rna_strand_detect(bam,annotation,max_records=1000L,
   sampling='whole_bam_random',sampling_fraction=.1),'exceeds max_records')
 expect_error(rna_strand_detect(bam,annotation,sampling='whole_bam_random',sampling_fraction=0))
})

test_that('Random sampling is opt-in and its cache key records all defaults', {
 # Avoid an experiment/file conversion: exercise native ensure with a tiny
 # stand-in df and real cache files, mocking only experiment path resolution.
 d <- tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE))
 bam <- file.path(d,'x.bam');writeLines('fixture',bam)
 annotation <- list(annotation_signature='fixture')
 calls <- 0L
 testthat::local_mocked_bindings(filepath=function(...)bam,.package='ORFik')
 testthat::local_mocked_bindings(rna_strand_detect=function(bam,annotation,...){
  calls <<- calls+1L
  p <- list(...)
  list(version=1L,product_signature=rna_strand_signature(bam),
   annotation_signature=annotation$annotation_signature,parameters=p,
   status='unstranded',read_strand='unstranded',strandMode=0L)
 })
 a <- rna_strand_ensure(NULL,annotation,list(sampling='whole_bam_random'))
 b <- rna_strand_ensure(NULL,annotation,list(sampling_seed=17L,sampling='whole_bam_random'))
 expect_identical(calls,1L)
 expect_identical(a$decision_id,b$decision_id)
 expect_equal(a$parameters$sampling_fraction,.1)
 rna_strand_ensure(NULL,annotation,list(sampling='whole_bam_random',sampling_seed=31L))
 expect_identical(calls,2L)
})
