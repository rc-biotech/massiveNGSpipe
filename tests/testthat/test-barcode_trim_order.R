fastp_fixture_binary <- function() {
  p<-Sys.which('fastp')
  if(!nzchar(p))p<-path.expand('~/bin/fastp')
  testthat::skip_if_not(file.exists(p),'fastp not available for sequence regression')
  p
}

test_that('actual fastp trimming preserves the known insert and removes post-adapter UMIs', {
  tool<-fastp_fixture_binary();root<-tempfile('trim order spaces ',tmpdir=path.expand('~'));dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  body<-'ATGGAATTTCCAGCTGACCTAGACCTTACAG'
  prefix<-'TGCAGACTTCAAGGGG';tail<-'AACCTTAGAC';adapter<-'AGATCGGAAGAGCACACGTCTGAA'
  seqs<-c(rep(paste0(prefix,body,tail,adapter,strrep('C',40)),100),
          rep(paste0(prefix,substr(body,1,9),tail,adapter,strrep('C',40)),5))
  file<-file.path(root,'SRR_fixture.fastq')
  lines<-unlist(lapply(seq_along(seqs),function(i)c(paste0('@read',i),seqs[i],'+',strrep('I',nchar(seqs[i])))))
  writeLines(lines,file)
  report<-fastp_post_adapter_tail(sub(path.expand('~'),'~',file,fixed=TRUE),output.dir=sub(path.expand('~'),'~',root,fixed=TRUE),adapter=adapter,trim.front=16,trim.tail=10,threads=2,fastp_path=tool)
  result<-readLines(file.path(root,'trim/trimmed_SRR_fixture.fastq'))
  expect_identical(result[seq(2,length(result),4)],rep(body,100))
  expect_equal(nchar(result[seq(4,length(result),4)]),rep(nchar(body),100))
  metadata<-jsonlite::fromJSON(report)
  expect_equal(metadata$summary$before_filtering$total_reads,105)
  expect_equal(metadata$summary$before_filtering$read1_mean_length,floor(mean(nchar(seqs))))
  expect_equal(metadata$summary$after_filtering$total_reads,100)
  expect_equal(metadata$summary$after_filtering$read1_mean_length,nchar(body))
  expect_equal(metadata$filtering_result$too_short_reads,5)
  expect_equal(metadata$massiveNGSpipe_trim$post_adapter_tail,10)
  expect_equal(metadata$massiveNGSpipe_trim$order,'adapter_then_tail_v1')
})

test_that('an upstream fastp error cannot publish a successful canonical report', {
  root<-tempfile('trim_failure_');dir.create(root);input<-file.path(root,'bad.fastq');writeLines(c('@bad','ACGT','+','I'),input)
  expect_error(fastp_post_adapter_tail(input,output.dir=root,adapter='AGATCGGAAGAGC',trim.tail=5,fastp_path=fastp_fixture_binary()),'Ordered fastp trimming failed')
  expect_false(file.exists(file.path(root,'trim/report_bad.json')))
})


test_that('ordered streaming preserves both paired-end sequences and mate counts', {
  tool<-fastp_fixture_binary();root<-tempfile('paired_trim_');dir.create(root)
  body1<-'ATGGAATTTCCAGCTGACCTAGACCTTACAG';body2<-'CGTCTGAAACCATGCTGGACTAGTCAGGATTC'
  prefix<-'TGCAGACTTCAAGGGG';tail<-'AACCTTAGAC';adapter<-'AGATCGGAAGAGCACACGTCTGAA'
  files<-file.path(root,c('SRR_pair_1.fastq','SRR_pair_2.fastq'))
  for(j in 1:2) {
    seq<-paste0(prefix,if(j==1)body1 else body2,tail,adapter,strrep('C',40))
    writeLines(unlist(lapply(1:20,function(i)c(paste0('@pair',i,'/',j),seq,'+',strrep('I',nchar(seq))))),files[j])
  }
  report<-fastp_post_adapter_tail(files[1],files[2],output.dir=root,adapter=adapter,trim.front=16,trim.tail=10,threads=2,fastp_path=tool)
  a<-readLines(file.path(root,'trim/trimmed_SRR_pair_1.fastq'))
  b<-readLines(file.path(root,'trim/trimmed2_SRR_pair_1.fastq'))
  expect_equal(a[seq(2,length(a),4)],rep(body1,20))
  expect_equal(b[seq(2,length(b),4)],rep(body2,20))
  expect_equal(jsonlite::fromJSON(report)$summary$after_filtering$total_reads,40)
})
