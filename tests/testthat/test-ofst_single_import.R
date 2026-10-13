test_that('OFST worker budget caps memory and sample count independently of other stages', {
 expect_equal(massiveNGSpipe:::ofst_worker_count(rep(1e8,12),8,1e11),8L)
 expect_equal(massiveNGSpipe:::ofst_worker_count(rep(1e8,2),8,1e11),2L)
 expect_equal(massiveNGSpipe:::ofst_worker_count(rep(1e9,12),8,4e10),1L)
 expect_equal(massiveNGSpipe:::ofst_worker_count(numeric(),8,1e11),1L)
 expect_error(massiveNGSpipe:::ofst_worker_count(NA_real_,8,1e11))
})
test_that('A resumed mapper conversion skips a mode already published and flagged', {
 df<-ORFik::ORFik.template.experiment()[1,];df@listData$Run<-'SRR_FIXTURE';exp<-ORFik::name(df);run<-ORFik::runIDs(df)
 config<-list(project=tempfile('mode_flags'),all_mappers=TRUE,split_unique_mappers=TRUE)
 on.exit(unlink(config$project,recursive=TRUE));calls<-list();attempt<-0L
 testthat::local_mocked_bindings(convert_bam_to_ofst_modes=function(df,modes,after_mode){
  calls[[length(calls)+1L]]<<-modes;attempt<<-attempt+1L
  for(m in modes){if(attempt==1L && m=='unique')stop('unique failed');after_mode(m,tempfile(),data.table::data.table(read_length=10,count=1,percent=100))}
 },.package='ORFik')
 expect_error(massiveNGSpipe:::convert_mapper_modes_per_sample(df,config,BiocParallel::SerialParam()),'unique failed')
 expect_equal(massiveNGSpipe:::samples_done(config,'ofst',exp),run)
 expect_length(massiveNGSpipe:::samples_done(config,'ofst_unique',exp),0)
 massiveNGSpipe:::convert_mapper_modes_per_sample(df,config,BiocParallel::SerialParam())
 expect_equal(calls[[2]],'unique');expect_equal(massiveNGSpipe:::samples_done(config,'ofst_unique',exp),run)
})
test_that('Single mapper configurations never create flags for the other mode', {
 df<-ORFik::ORFik.template.experiment()[1,];df@listData$Run<-'SRR_FIXTURE';exp<-ORFik::name(df);seen<-character()
 testthat::local_mocked_bindings(convert_bam_to_ofst_modes=function(df,modes,after_mode){seen<<-modes;for(m in modes)after_mode(m,tempfile(),data.table::data.table(read_length=10,count=1,percent=100))},.package='ORFik')
 for(all in c(TRUE,FALSE)) {
  config<-list(project=tempfile('one_mode'),all_mappers=all,split_unique_mappers=!all)
  massiveNGSpipe:::convert_mapper_modes_per_sample(df,config,BiocParallel::SerialParam())
  expect_equal(seen,if(all)'all' else 'unique')
  expect_length(massiveNGSpipe:::samples_done(config,if(all)'ofst_unique' else 'ofst',exp),0)
  unlink(config$project,recursive=TRUE)
 }
})
