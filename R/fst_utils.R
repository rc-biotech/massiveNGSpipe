
generateSingleMatrix <- function(grl, trx_id = names(grl), reads,libnames, output_dir, BPPARAM){
  stopifnot(class(grl) %in% c("GRangesList", "CompressedGRangesList"))
  stopifnot(length(grl) == 1)
  tx_coverage <- bplapply(reads, function(x) coveragePerTiling(grl,x),
                          BPPARAM = BPPARAM)

  names(tx_coverage) <- libnames
  tx_coverage <- rbindlist(tx_coverage, idcol = "library")[,c("library","count")][, library := as.factor(library)]

  output_path <- collection_path_from_exp(NULL, trx_id, must_exists = FALSE,
                                          collection_dir = output_dir)
  fst::write.fst(tx_coverage, output_path)
  rm(tx_coverage)
  gc()
  message(paste("done", trx_id))
  return(invisible(NULL))
}

#' General fst all samples tables
#' @param df ORFik experiment of all samples per species
#' @param trxs transcript ids, default names(grl)
#' @param output_dir RiboCrypt::collection_dir_from_exp(df)
#' @param grl loadRegion(df, "transcript")
#' @param libnames = bamVarName(df), the names of columns for fst
#' @param verbose logical, TRUE. Print name of current tx.
#' @param BPPARAM BiocParallel::bpparam(), for single thread,
#'  do: BiocParallel::SerialParam()
#' @return invisible(NULL), results saved to disc.
#' @importFrom RiboCrypt collection_dir_from_exp collection_path_from_exp
#' @export
generateAllTranscriptsMatrix <- function(df, trxs = names(grl),
                                         output_dir = RiboCrypt::collection_dir_from_exp(df),
                                         grl = loadRegion(df, "mrna"),
                                         libnames = bamVarName(df),
                                         verbose = TRUE,
                                         BPPARAM = BiocParallel::bpparam()) {
  message("Settings:")
  message("- Species: ", organism(df))
  message("- Samples: ", nrow(df))
  message("- transcripts (inserted): ", length(trxs))
  stopifnot(all(trxs %in% names(grl)))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  trxs_done <- gsub("\\.fst", "", list.files(output_dir, "\\.fst$"))
  message("- transcripts (already complete): ", length(trxs_done))
  trxs <- trxs[!(trxs %in% trxs_done)]
  message("- transcripts (to run now): ", length(trxs))

  reads <- filepath(df, "bigwig", suffix_stem = c("_pshifted",""),
                    base_folders = libFolder(df, mode = "all"))
  stopifnot(length(reads) == length(libnames))

  lapply(trxs, function(x) {
    x <- grl[x]
    generateSingleMatrix(grl=x, names(x),reads,libnames, output_dir, BPPARAM)
  })


  return(invisible(NULL))
}
