test_that('single zero-length mate output is normalized before scratch copy', {
  env<-new.env(parent=globalenv());sys.source(testthat::test_path('../../R/download_sra.R'),env)
  d<-tempfile();dir.create(d);writeLines('source',file.path(d,'SRR001'))
  bin<-tempfile();dir.create(bin);file.create(file.path(bin,'fasterq-dump'))
  env$tmp_claim_root<-function()tempdir();env$fasterq_tmp_need_gb<-function(...)1
  env$tmp_space_claim<-function(...)'claim';env$tmp_space_release<-function(...)invisible(NULL)
  env$system2<-function(command,args,...){writeLines(c('@r','AC','+','II'),paste0(args[1],'_1.fastq'));file.create(paste0(args[1],'_2.fastq'));0L}
  expect_true(env$sra_to_fastq('SRR001',d,FALSE,bin,FALSE))
  expect_true(file.exists(file.path(d,'SRR001.fastq')))
  expect_identical(readLines(file.path(d,'SRR001.fastq')),c('@r','AC','+','II'))
  expect_true(file.exists(file.path(d,'SRR001')))
})

test_that('a lone _1 FASTQ (no _2 written) is normalized on the output drive', {
  env<-new.env(parent=globalenv());sys.source(testthat::test_path('../../R/download_sra.R'),env)
  d<-tempfile();dir.create(d);writeLines('source',file.path(d,'SRR001'))
  bin<-tempfile();dir.create(bin);file.create(file.path(bin,'fasterq-dump'))
  env$tmp_claim_root<-function()tempdir();env$fasterq_tmp_need_gb<-function(...)1
  env$tmp_space_claim<-function(...)NULL
  env$system2<-function(command,args,...){writeLines(c('@r','AC','+','II'),paste0(args[1],'_1.fastq'));0L}
  expect_true(env$sra_to_fastq('SRR001',d,FALSE,bin,FALSE))
  expect_true(file.exists(file.path(d,'SRR001.fastq')))
  expect_false(file.exists(file.path(d,'SRR001_1.fastq')))
})

test_that('a non-empty second read for SINGLE metadata is kept for validate_fastq_download (read 1 rule)', {
  env<-new.env(parent=globalenv());sys.source(testthat::test_path('../../R/download_sra.R'),env)
  d<-tempfile();dir.create(d);writeLines('source',file.path(d,'SRR001'))
  bin<-tempfile();dir.create(bin);file.create(file.path(bin,'fasterq-dump'))
  env$tmp_claim_root<-function()tempdir();env$fasterq_tmp_need_gb<-function(...)1
  env$tmp_space_claim<-function(...)NULL;env$system2<-function(command,args,...){for(p in paste0(args[1],c('_1.fastq','_2.fastq')))writeLines('x',p);0L}
  expect_true(env$sra_to_fastq('SRR001',d,FALSE,bin,FALSE))
  expect_true(all(file.exists(file.path(d,c('SRR001_1.fastq','SRR001_2.fastq')))))
  expect_true(file.exists(file.path(d,'SRR001')))
})
