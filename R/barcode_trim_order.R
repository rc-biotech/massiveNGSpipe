# Fixed fastp tail cuts run before adapter removal. Stream adapter/front
# trimming into a second fastp process so the 3-prime cut acts on the insert.
# The pipe avoids writing/re-reading a large adapter-only FASTQ intermediate.
fastp_post_adapter_tail <- function(file, file2 = NULL, output.dir, adapter,
                                    trim.front = 0, trim.tail, min.length = 20,
                                    threads = 4L, fastp_path = install.fastp()) {
  file <- path.expand(file)
  output.dir <- path.expand(output.dir)
  if (!is.null(file2)) file2 <- path.expand(file2)
  paired <- !is.null(file2) && length(file2) && !is.na(file2)
  stopifnot(length(file)==1L,file.exists(file),trim.tail>0,trim.front>=0,
            threads>=1, min.length>=1)
  if(paired)stopifnot(file.exists(file2))
  fastp_path<-path.expand(fastp_path)
  folder<-file.path(output.dir,'trim');dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  id<-tools::file_path_sans_ext(sub('\\.gz$','',basename(file)))
  metrics<-file.path(folder,'trim_order_metrics');dir.create(metrics,showWarnings=FALSE)
  first_json<-file.path(metrics,paste0(id,'_adapter.json'))
  second_json<-file.path(metrics,paste0(id,'_tail.json'))
  final_json<-file.path(folder,paste0('report_',id,'.json'))
  first<-c('--in1',file,if(paired)c('--in2',file2),'--stdout',
           '--trim_front1',as.character(trim.front),
           if(paired)c('--trim_front2',as.character(trim.front)),
           '--length_required',as.character(min.length),'--thread',as.character(threads),
           '--disable_quality_filtering','--json',first_json,
           '--html',file.path(metrics,paste0(id,'_adapter.html')))
  first<-c(first,if(adapter=='auto')character() else if(adapter=='disable')'--disable_adapter_trimming' else c('--adapter_sequence',adapter,if(paired)c('--adapter_sequence_r2',adapter)))
  second<-c('--stdin',if(paired)'--interleaved_in',
            '--out1',file.path(folder,paste0('trimmed_',id,'.fastq')),
            if(paired)c('--out2',file.path(folder,paste0('trimmed2_',id,'.fastq'))),
            '--disable_adapter_trimming','--trim_tail1',as.character(trim.tail),
            if(paired)c('--trim_tail2',as.character(trim.tail)),
            '--length_required',as.character(min.length),'--thread',as.character(threads),
            '--disable_quality_filtering','--json',second_json,
            '--html',file.path(folder,paste0('report_',id,'.html')))
  cmd<-paste(paste(shQuote(c(fastp_path,first)),collapse=' '),'|',
             paste(shQuote(c(fastp_path,second)),collapse=' '))
  message('Adapter/front trimming then post-adapter tail removal: ',id)
  status<-suppressWarnings(system2('bash',c('-o','pipefail','-c',shQuote(cmd))))
  if(status!=0L)stop('Ordered fastp trimming failed (status ',status,') for ',file)
  stopifnot(file.exists(first_json),file.exists(second_json))
  before<-jsonlite::fromJSON(first_json,simplifyVector=FALSE)
  after<-jsonlite::fromJSON(second_json,simplifyVector=FALSE)
  # Canonical pipeline JSON describes original input -> final insert. Detailed
  # stage reports remain separate; the HTML describes the post-adapter stage.
  after$summary$before_filtering<-before$summary$before_filtering
  after$summary$sequencing<-before$summary$sequencing
  after$read1_before_filtering<-before$read1_before_filtering
  if(paired)after$read2_before_filtering<-before$read2_before_filtering
  after$adapter_cutting<-before$adapter_cutting
  for(k in setdiff(union(names(before$filtering_result),names(after$filtering_result)),
                   'passed_filter_reads')) {
    a<-before$filtering_result[[k]];b<-after$filtering_result[[k]]
    if(is.numeric(a)||is.numeric(b))after$filtering_result[[k]]<-sum(c(a,b))
  }
  after$command<-cmd
  after$massiveNGSpipe_trim<-list(order='adapter_then_tail_v1',front=trim.front,
    post_adapter_tail=trim.tail,paired=paired,adapter_stage=first_json,tail_stage=second_json)
  temp_json<-tempfile(pattern=paste0(id,'_report_'),tmpdir=folder)
  on.exit(unlink(temp_json),add=TRUE)
  jsonlite::write_json(after,temp_json,auto_unbox=TRUE,pretty=TRUE)
  if(!file.rename(temp_json,final_json))stop('Could not publish complete trimming report for ',id)
  invisible(final_json)
}
