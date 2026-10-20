provenance_fixture <- function() {
  root <- tempfile('shift_binding_'); dir.create(root)
  rawdir <- file.path(root,'ofst'); out <- file.path(root,'pshifted')
  dir.create(rawdir); dir.create(out)
  ids <- c('SRR001','SRR002'); raw <- setNames(file.path(rawdir,paste0(ids,'.ofst')),ids)
  post <- setNames(file.path(out,paste0(ids,'_pshifted.ofst')),ids)
  for (p in c(raw,post)) writeLines('fixture',p)
  for (p in post) Sys.setFileTime(p,Sys.time()+10)
  table <- file.path(out,'shifting_table.rds')
  offsets <- data.frame(fraction=c(27L,28L),offsets_start=c(-12L,-13L))
  shifts <- setNames(rep(list(offsets),2),unname(raw)); saveRDS(shifts,table)
  list(df=fake_experiment_stub(run_ids=ids),raw=raw,post=post,table=table,shifts=shifts)
}
mock_provenance_paths <- function(f, env=parent.frame()) {
  testthat::local_mocked_bindings(ofst_filepath=function(df)unname(f$raw[runIDs(df)]),
    pshifted_filepath=function(df)unname(f$post[runIDs(df)]),.env=env)
}

test_that('canonical offset identity ignores row order and data.table attributes', {
  a<-data.frame(fraction=c(29,28),offsets_start=c(-13,-12));b<-data.table::as.data.table(a[2:1,]);b[,extra:='ignored']
  expect_identical(pshift_canonical_offsets(a),pshift_canonical_offsets(b))
  expect_error(pshift_canonical_offsets(data.frame(fraction=c(28,28),offsets_start=-12)),'unique')
  expect_error(pshift_canonical_offsets(data.frame(fraction=28,offsets_start=1)),'nonpositive')
})

test_that('explicit edit invalidates even when raw and P-site mtimes would pass', {
  f<-provenance_fixture();mock_provenance_paths(f)
  expect_true(pshift_sample_valid(f$df[1,]))
  changed<-f$shifts[1];changed[[1]]$offsets_start<-c(-13L,-12L)
  pshift_mark_offsets_pending(f$table,changed,'approved length-specific fix')
  expect_false(pshift_sample_valid(f$df[1,]))
  expect_true(pshift_sample_valid(f$df[2,]))
})

test_that('applied identity detects offset edits but preserves unchanged siblings', {
  f<-provenance_fixture();mock_provenance_paths(f)
  pshift_record_applied_offsets(f$df)
  expect_true(pshift_sample_valid(f$df[1,]));expect_true(pshift_sample_valid(f$df[2,]))
  changed<-f$shifts;changed[[2]]$offsets_start<-c(-13L,-12L);saveRDS(changed,f$table)
  expect_true(pshift_sample_valid(f$df[1,]));expect_false(pshift_sample_valid(f$df[2,]))
})

test_that('wrong saved offset or missing product cannot clear a pending request', {
  f<-provenance_fixture();mock_provenance_paths(f)
  desired<-f$shifts[1];desired[[1]]$offsets_start<-c(-13L,-12L)
  pshift_mark_offsets_pending(f$table,desired,'test')
  pending<-pshift_offset_record_paths(f$raw[1],f$post[1])$pending
  expect_error(pshift_record_applied_offsets(f$df[1,]),'differ from explicit request')
  expect_true(file.exists(pending));expect_false(pshift_sample_valid(f$df[1,]))
  unlink(f$post[1]);expect_error(pshift_record_applied_offsets(f$df[1,]),'missing')
  expect_true(file.exists(pending))
})

test_that('successful matching native output clears only its own pending request', {
  f<-provenance_fixture();mock_provenance_paths(f)
  pshift_mark_offsets_pending(f$table,f$shifts,'test')
  writeLines('new native output',f$post[1])
  pshift_record_applied_offsets(f$df[1,])
  expect_true(pshift_sample_valid(f$df[1,]));expect_false(pshift_sample_valid(f$df[2,]))
  writeLines('new native output',f$post[2]);pshift_record_applied_offsets(f$df[2,]);expect_true(pshift_sample_valid(f$df[2,]))
})

test_that('mapper-mode provenance lives beside its own P-sites', {
  all<-pshift_offset_record_paths('/x/ofst/SRR1.ofst','/x/pshifted/SRR1_pshifted.ofst')
  unique<-pshift_offset_record_paths('/x/unique_mappers/ofst/SRR1.ofst','/x/unique_mappers/pshifted/SRR1_pshifted.ofst')
  expect_false(identical(all$pending,unique$pending));expect_false(identical(all$applied,unique$applied))
})

test_that('changing the raw source after approval keeps the request pending', {
  f<-provenance_fixture();mock_provenance_paths(f)
  pshift_mark_offsets_pending(f$table,f$shifts[1],'test')
  writeLines('different input',f$raw[1])
  expect_error(pshift_record_applied_offsets(f$df[1,]),'source differ')
  expect_false(pshift_offset_binding_valid(f$df[1,]))
})

test_that('standard pipeline re-applies an explicit edit instead of trusting mtimes', {
  f<-provenance_fixture();mock_provenance_paths(f)
  config<-fake_config(preset='Ribo-seq',extra=list(all_mappers=TRUE,split_unique_mappers=FALSE,
    reuse_shifts_if_existing=TRUE,accepted_lengths_rpf=c(27,28),max_no_adapter_removed_pct=80))
  fake_mark_all_done(config,c('trim','collapsed','aligned','cleanbam','exp','ofst'),name(f$df))
  wanted<-f$shifts;wanted[[1]]$offsets_start<-c(-13L,-12L)
  pshift_mark_offsets_pending(f$table,wanted[1],'test');saveRDS(wanted,f$table)
  called<-character()
  testthat::local_mocked_bindings(shifts_load_safe=function(df,...)wanted,
    shiftFootprintsByExperimentSafe=function(df,shifting_table,...){called<<-runIDs(df);saveRDS(wanted,f$table);writeLines('new native output',f$post[runIDs(df)]);invisible(NULL)},
    save_pshifted_length_distributions=function(df)invisible(NULL))
  pipeline_pshift(list(f$df),config)
  expect_identical(called,'SRR001');expect_true(pshift_sample_valid(f$df[1,]));expect_true(pshift_sample_valid(f$df[2,]))
})

test_that('unchanged old product cannot be certified as a newly applied request', {
  f<-provenance_fixture();mock_provenance_paths(f)
  pshift_mark_offsets_pending(f$table,f$shifts[1],'test')
  expect_error(pshift_record_applied_offsets(f$df[1,]),'not rewritten')
  expect_false(pshift_sample_valid(f$df[1,]))
})

test_that('explicit staged mitochondrial offsets are not overridden by automatic mito detection', {
  f<-provenance_fixture();mock_provenance_paths(f);one<-f$df[1,]
  pshift_mark_offsets_pending(f$table,f$shifts[1],'explicit mito anchor review')
  auto_called<-FALSE
  testthat::local_mocked_bindings(is_mito_sample=function(df)rep(TRUE,nrow(df)),
    mito_shift_tables=function(df){auto_called<<-TRUE;stop('must not override explicit request')},
    shifts_load_safe=function(df,...)f$shifts[1],
    shiftFootprintsByExperimentSafe=function(df,...){writeLines('fresh mito native output',f$post[1]);NULL},
    save_pshifted_length_distributions=function(df)NULL)
  pshift_needed_subset(one,c(27,28),character(),BiocParallel::SerialParam(),80,TRUE)
  expect_false(auto_called);expect_true(pshift_sample_valid(one))
})

test_that('automatic mitochondrial detection remains active without an explicit request', {
  f<-provenance_fixture();mock_provenance_paths(f);one<-f$df[1,];unlink(f$post[1])
  auto_called<-FALSE
  testthat::local_mocked_bindings(is_mito_sample=function(df)rep(TRUE,nrow(df)),
    mito_shift_tables=function(df){auto_called<<-TRUE;f$shifts[1]},
    shifts_load_safe=function(df,...)f$shifts[1],
    shiftFootprintsByExperimentSafe=function(df,...){writeLines('fresh mito native output',f$post[1]);NULL},
    save_pshifted_length_distributions=function(df)NULL)
  pshift_needed_subset(one,c(27,28),character(),BiocParallel::SerialParam(),80,TRUE)
  expect_true(auto_called);expect_true(pshift_sample_valid(one))
})

test_that('extra-length production writer refreshes provenance, scratch writer does not', {
  f<-provenance_fixture();df<-f$df;out<-tempfile('extension_');dir.create(out)
  ext<-f$shifts;attr(ext,'added')<-data.frame(sample=names(ext)[1],length=23L,offset=-12L)
  records<-0L
  testthat::local_mocked_bindings(read.experiment=function(...)df,shifts_load_safe=function(...)f$shifts,
    extend_shift_table=function(...)ext,libFolder=function(...)out,name=function(...)'test',
    pshift_record_applied_offsets=function(df){records<<-records+1L},save_pshifted_length_distributions=function(...)NULL)
  testthat::local_mocked_bindings(shiftFootprintsByExperiment=function(...)NULL,.package='ORFik')
  pshift_add_extra_lengths('test',23L);expect_identical(records,1L)
  pshift_add_extra_lengths('test',23L,out_dir=tempfile('scratch_'));expect_identical(records,1L)
})

provenance_recipe <- function(lengths = c(27,28)) pshift_detection_recipe(lengths, character(), 80)

test_that('automatic reuse requires bound unchanged source and recipe, not table existence', {
  f<-provenance_fixture();mock_provenance_paths(f);one<-f$df[1,];r<-provenance_recipe()
  expect_null(pshift_reusable_offsets(one,r)) # legacy entry
  pshift_record_applied_offsets(one,recipe=r)
  expect_identical(pshift_reusable_offsets(one,r),pshift_canonical_offsets(f$shifts[[1]]))
  unlink(f$post[1]) # regeneration on unchanged source can safely reuse
  expect_identical(pshift_reusable_offsets(one,r),pshift_canonical_offsets(f$shifts[[1]]))
  expect_null(pshift_reusable_offsets(one,r,FALSE))
  writeLines('new retrimmed raw input',f$raw[1])
  expect_null(pshift_reusable_offsets(one,r))
})

test_that('mask expansion invalidates new resume and automatic reuse certificates', {
  f<-provenance_fixture();mock_provenance_paths(f);one<-f$df[1,];r<-provenance_recipe()
  pshift_record_applied_offsets(one,recipe=r)
  expect_true(pshift_recipe_valid(one,provenance_recipe(c(28,27,27))))
  expect_false(pshift_recipe_valid(one,provenance_recipe(20:40)))
  expect_null(pshift_reusable_offsets(one,provenance_recipe(20:40)))
  pshift_record_applied_offsets(one) # extra length writer keeps same-source recipe
  expect_false(pshift_recipe_valid(one,provenance_recipe(20:40)))
})

test_that('explicit requests survive no-reuse, but stale/corrupt requests stop before writing', {
  f<-provenance_fixture();mock_provenance_paths(f);one<-f$df[1,];r<-provenance_recipe()
  pshift_mark_offsets_pending(f$table,f$shifts[1],'reviewed')
  expect_identical(pshift_reusable_offsets(one,r,FALSE),pshift_canonical_offsets(f$shifts[[1]]))
  changed<-f$shifts;changed[[1]]$offsets_start[1]<- -13L;saveRDS(changed,f$table)
  expect_error(pshift_reusable_offsets(one,r,FALSE),'differ from explicit request')
  saveRDS(f$shifts,f$table);writeLines('new source',f$raw[1])
  expect_error(pshift_reusable_offsets(one,r,FALSE),'differ from explicit request')
})

test_that('retrimming enters fresh native detection and preserves a valid sibling table', {
  f<-provenance_fixture();mock_provenance_paths(f);r<-provenance_recipe()
  pshift_record_applied_offsets(f$df,recipe=r)
  writeLines('changed raw source',f$raw[1]);got<-list();fresh<-f$shifts[1]
  fresh[[1]]$offsets_start<-c(-13L,-12L)
  testthat::local_mocked_bindings(is_mito_sample=function(df)rep(FALSE,nrow(df)),
    shiftFootprintsByExperimentSafe=function(df,shifting_table,...){
      got[[length(got)+1L]]<<-list(ids=runIDs(df),shifts=shifting_table)
      saveRDS(fresh,f$table);writeLines('fresh detected native output',f$post[1]);NULL},
    save_pshifted_length_distributions=function(df)NULL)
  pshift_needed_subset(f$df,c(27,28),character(),BiocParallel::SerialParam(),80,TRUE)
  expect_identical(got[[1]]$ids,'SRR001');expect_null(got[[1]]$shifts)
  expect_identical(readRDS(f$table)[[f$raw[2]]],f$shifts[[2]])
  expect_true(pshift_sample_valid(f$df[1,]));expect_true(pshift_sample_valid(f$df[2,]))
})

test_that('no-reuse native route still applies explicit mitochondrial requests', {
  f<-provenance_fixture();mock_provenance_paths(f);one<-f$df[1,];got<-NULL
  pshift_mark_offsets_pending(f$table,f$shifts[1],'reviewed mito')
  testthat::local_mocked_bindings(is_mito_sample=function(df)TRUE,
    mito_shift_tables=function(df)stop('explicit request must override automatic mito'),
    shiftFootprintsByExperimentSafe=function(df,shifting_table,...){got<<-shifting_table;writeLines('new product',f$post[1]);NULL},
    save_pshifted_length_distributions=function(df)NULL)
  pshift_needed_subset(one,c(27,28),character(),BiocParallel::SerialParam(),80,FALSE)
  expect_identical(got[[1]],pshift_canonical_offsets(f$shifts[[1]]))
  expect_true(pshift_sample_valid(one))
})

test_that('mixed explicit and fresh detection use separate native calls without losing entries', {
  f<-provenance_fixture();mock_provenance_paths(f);got<-list()
  pshift_mark_offsets_pending(f$table,f$shifts[1],'reviewed')
  unlink(f$post[2]) # legacy missing product must detect fresh
  testthat::local_mocked_bindings(is_mito_sample=function(df)rep(FALSE,nrow(df)),
    shiftFootprintsByExperimentSafe=function(df,shifting_table,...){
      ids<-runIDs(df);got[[length(got)+1L]]<<-list(ids=ids,known=!is.null(shifting_table))
      saveRDS(f$shifts[unname(f$raw[ids])],f$table)
      writeLines('new native mixed product',f$post[ids]);NULL},
    save_pshifted_length_distributions=function(df)NULL)
  pshift_needed_subset(f$df,c(27,28),character(),BiocParallel::SerialParam(),80,FALSE)
  expect_identical(got,list(list(ids='SRR001',known=TRUE),list(ids='SRR002',known=FALSE)))
  expect_setequal(names(readRDS(f$table)),unname(f$raw))
  expect_true(pshift_sample_valid(f$df[1,]));expect_true(pshift_sample_valid(f$df[2,]))
})
