# Staging new runs for the metadata sheets (batch curation helpers, 2026-10-21).
#
# Adding many studies at once (e.g. a database search hand-off) repeats the same steps: fetch SRA metadata, keep the
# requested runs, auto-annotate, normalise organism names, drop runs we already have or that are re-deposits of a
# sample we have under another BioProject, make every sample unique on the curated columns, append to the sheet and
# read the sheet back. These helpers do that in one place; curation of the labels themselves stays a review step
# (rows are marked in CHECKED).

#' Authenticate googledrive / googlesheets4 from a cached gargle token, for non-interactive scripts
#'
#' The package helpers call \code{drive_auth(email = ...)}; in a non-interactive session that can fail even with a
#' valid cached token. This loads the token file once and makes later \code{drive_auth()} / \code{gs4_auth()} calls
#' reuse it for the rest of the session.
#' @param token_file path to a cached token (rds), default: the newest file in the gargle cache whose name contains
#'  \code{pattern}
#' @param pattern character, file name pattern in the gargle cache
#' @param cache_dir gargle cache directory
#' @return invisible(token)
#' @export
google_auth_cached <- function(token_file = NULL, pattern = "", cache_dir = "~/.cache/gargle") {
  if (is.null(token_file)) {
    f <- list.files(cache_dir, pattern, full.names = TRUE)
    f <- f[!dir.exists(f)]
    if (!length(f)) stop("No cached gargle token in ", cache_dir)
    token_file <- f[which.max(file.mtime(f))]
  }
  tok <- readRDS(token_file)
  googledrive::drive_auth(token = tok); googlesheets4::gs4_auth(token = tok)
  rd <- googledrive::drive_auth; rg <- googlesheets4::gs4_auth
  utils::assignInNamespace("drive_auth", function(email = NULL, token = NULL, ...) rd(token = tok), "googledrive")
  utils::assignInNamespace("gs4_auth", function(email = NULL, token = NULL, ...) rg(token = tok), "googlesheets4")
  utils::assignInNamespace("token_fetch", function(...) tok, "gargle")
  invisible(tok)
}

#' Curated columns that make a sample unique within a study and organism
#' @noRd
metadata_unique_columns <- function() c("LIBRARYTYPE", "CELL_LINE", "TISSUE", "CONDITION", "GENE", "INHIBITOR",
                                        "TIMEPOINT", "FRACTION", "BATCH", "REPLICATE")

#' Make samples unique on the curated columns within study x organism
#'
#' Duplicated label sets first get REPLICATE 1..n (only where REPLICATE is empty); what is still duplicated gets an
#' \code{auto_<n>} FRACTION placeholder (appended to an existing FRACTION), the same convention as the pipeline's
#' validation. Placeholders mark rows to curate, not biology.
#' @param dt data.table with study_accession, ScientificName and the curated columns (missing columns are added empty)
#' @return dt (modified copy) with character curated columns
#' @export
metadata_make_unique <- function(dt) {
  dt <- copy(as.data.table(dt)); k <- metadata_unique_columns()
  for (cl in k) {
    if (!cl %in% names(dt)) dt[, (cl) := ""]
    set(dt, j = cl, value = { v <- as.character(dt[[cl]]); v[is.na(v)] <- ""; v })
  }
  key <- function(d) do.call(paste, c(lapply(c("study_accession", "ScientificName", k),
                                             function(cl) tolower(as.character(d[[cl]]))), sep = "|"))
  dt[, .key := key(dt)][, .n := .N, by = .key]
  dt[.n > 1 & REPLICATE == "", REPLICATE := as.character(seq_len(.N)), by = .key]
  dt[, .key := key(dt)][, .n := .N, by = .key]
  dt[.n > 1, FRACTION := ifelse(FRACTION == "", paste0("auto_", seq_len(.N)), paste0(FRACTION, "_auto_", seq_len(.N))),
     by = .key]
  dt[, c(".key", ".n") := NULL]
  stopifnot(!anyDuplicated(key(dt)))
  dt[]
}

#' Runs whose BioSample already belongs to a run of another BioProject (re-deposits / duplicate uploads)
#' @param dt data.table with Run, BioProject, BioSample (new rows)
#' @param existing data.table with BioProject, BioSample (rows we have, e.g. all sheets)
#' @return data.table(Run, BioProject, BioSample, other_bioprojects)
#' @export
metadata_biosample_duplicates <- function(dt, existing) {
  empty <- data.table(Run = character(), BioProject = character(), BioSample = character(), other_bioprojects = character())
  ex <- unique(as.data.table(existing)[!is.na(BioSample) & BioSample != "", .(BioSample, BioProject)])
  d <- as.data.table(dt)[!is.na(BioSample) & BioSample != "" & BioSample %in% ex$BioSample, .(Run, BioProject, BioSample)]
  if (!nrow(d)) return(empty)
  d[, other_bioprojects := vapply(seq_len(.N), function(i)
    paste(setdiff(ex[BioSample == d$BioSample[i]]$BioProject, d$BioProject[i]), collapse = ";"), character(1))]
  d[other_bioprojects != ""]
}

#' Assay hints from sample titles that contradict a footprint (RFP) library
#'
#' Search hand-offs mark whole studies as Ribo-seq; titles often show disome, small-subunit, RNA/input or
#' immunoprecipitation libraries inside them (e.g. "Disome 60mers, WT, rep1", "MC4100 ..., IP").
#' @param titles character, sample titles
#' @return character, "" or one of "disome", "SSU", "RNA", "IP"
#' @export
metadata_title_assay <- function(titles) {
  t <- tolower(ifelse(is.na(titles), "", titles))
  out <- rep("", length(t))
  out[grepl("\\b(di|tri)-?somes?\\b|60 ?mers?|disome", t)] <- "disome"
  out[out == "" & grepl("\\b40s\\b|\\bssu\\b|small subunit|tcp-?seq|scanning", t)] <- "SSU"
  ribo <- grepl("translatome|ribosome profil|ribo-?seq|footprint|\\brpfs?\\b", t)
  out[out == "" & !ribo & grepl("rna-?seq|\\bmrna\\b|\\binput\\b|total rna|poly-?a\\b", t)] <- "RNA"
  out[out == "" & grepl("(^|[^a-z-])(ip|rip|chip)([^a-z-]|$)|immunoprecip|pull-?down|\\btrap\\b", t)] <- "IP"
  out
}

#' Fill CONDITION and REPLICATE from sample titles where the automatic annotation left them empty
#'
#' Within each study x organism: a trailing replicate number ("rep1", "_2", "-3", "Rep 1", "r2") is moved to
#' REPLICATE (if empty), the words shared by all titles of the study are removed, and the rest becomes CONDITION
#' (if CONDITION is empty for the whole study and the rest varies between samples). Example: "w-ct4-1" /
#' "bmal2-ct4-1" -> CONDITION "w" / "bmal2", REPLICATE 1. Rows changed this way still need review.
#' @param dt data.table with study_accession, ScientificName, sample_title, CONDITION, REPLICATE
#' @return dt (copy)
#' @export
metadata_title_labels <- function(dt) {
  dt <- copy(as.data.table(dt))
  for (cl in c("CONDITION", "REPLICATE", "sample_title")) {
    if (!cl %in% names(dt)) dt[, (cl) := ""]
    set(dt, j = cl, value = { v <- as.character(dt[[cl]]); v[is.na(v)] <- ""; v })
  }
  rep_re <- "(?i)[ _.-]*(?:rep(?:licate)?|r|bio|br)?[ _.-]*(\\d{1,2})$"
  dt[, .grp := paste(study_accession, ScientificName)]
  for (g in unique(dt$.grp)) {
    i <- which(dt$.grp == g); t <- trimws(dt$sample_title[i])
    if (length(i) < 2 || any(t == "")) next
    rep <- ifelse(grepl(rep_re, t, perl = TRUE), sub(paste0(".*?", rep_re), "\\1", t, perl = TRUE), "")
    base <- trimws(sub(rep_re, "", t, perl = TRUE))
    tok <- strsplit(tolower(base), "[^a-z0-9()+]+")
    shared <- Reduce(intersect, tok)
    # replicate-like and bare numeric tokens never describe the condition
    drop_tok <- function(x) x[!x %in% shared & nzchar(x) & !grepl("^(rep(licate)?|r|br|bio)?[0-9]+$", x)]
    rest <- vapply(tok, function(x) paste(drop_tok(x), collapse = "_"), character(1))
    cond <- dt$CONDITION[i]
    # use the title part when CONDITION is empty, or when one CONDITION value covers samples whose titles differ
    # (lossy annotation, e.g. "mock" and "WT" both annotated WT)
    lossy <- any(tapply(rest, cond, function(x) length(unique(x))) > 1)
    if (length(unique(rest)) > 1 && all(nzchar(rest)) && (all(cond == "") || lossy)) set(dt, i, "CONDITION", rest)
    fill <- dt$REPLICATE[i] == "" & nzchar(rep)
    if (any(fill)) set(dt, i[fill], "REPLICATE", rep[fill])
  }
  dt[, .grp := NULL][]
}

#' Stage new runs for a sheet: fetch, annotate, normalise, filter and make unique
#'
#' @param runs data.table with at least Run and BioProject (the runs to add)
#' @param config the mNGSp config (online mode)
#' @param libtype library type to assign (e.g. "RFP")
#' @param sheets data.table of rows we already have (all sheets), with Run, BioProject, BioSample; runs in it are
#'  dropped and its BioSamples are used for the re-deposit check
#' @param organism_map named character, extra ScientificName mapping applied after organism_name_cleanup()
#' @param keep_types library types from pipeline_metadata_annotate() accepted as \code{libtype}: an annotation that
#'  says e.g. "disome" or "RNA" is dropped for review instead of being overwritten
#' @param checked value for CHECKED, marks auto-curated rows
#' @param comment text appended to COMMENT
#' @param mito_runs runs that get FRACTION "mito" when FRACTION is empty (mitoribosome profiling)
#' @param biosample_labels logical, fill empty CELL_LINE / TISSUE from BioSample attributes (author-declared)
#' @param biosample_cache folder for cached BioSample records
#' @param vocab optional vocabulary (vocabulary_load()) to canonicalise filled values
#' @return list(rows = data.table, dropped = data.table(reason, Run, ...))
#' @export
stage_new_runs <- function(runs, config, libtype = "RFP", sheets = data.table(Run = character(), BioProject = character(),
                           BioSample = character()), organism_map = character(), keep_types = c("", libtype),
                           checked = "auto", comment = "", mito_runs = character(), biosample_labels = TRUE,
                           biosample_cache = file.path(config$project, "biosample_cache"), vocab = NULL) {
  runs <- as.data.table(runs)
  # one BioProject at a time: an SRA error for one study (empty error / throttling) must not stop the others
  fetch <- lapply(unique(runs$BioProject), function(b) tryCatch(
    pipeline_metadata(b, config, BPPARAM = BiocParallel::SerialParam()),
    error = function(e) { message("- metadata failed for ", b, ": ", conditionMessage(e)); b }))
  failed <- unlist(Filter(is.character, fetch))
  md <- rbindlist(Filter(is.data.frame, fetch), fill = TRUE)
  if (!nrow(md)) stop("SRA metadata failed for all ", length(fetch), " BioProjects")
  dropped <- list(metadata_failed = runs[BioProject %in% failed, .(Run, BioProject)],
                  not_in_sra = data.table(Run = setdiff(runs[!(BioProject %in% failed)]$Run, md$Run)))
  md <- md[Run %in% runs$Run][!duplicated(Run)]
  md <- pipeline_metadata_annotate(md)[!duplicated(Run)]
  bad_type <- !(as.character(md$LIBRARYTYPE) %in% keep_types | is.na(md$LIBRARYTYPE))
  dropped$annotated_other_type <- md[bad_type, .(Run, BioProject, annotated = LIBRARYTYPE, sample_title)]
  md <- md[!bad_type]
  ta <- metadata_title_assay(md$sample_title)
  other <- ta != "" & ta != libtype
  dropped$title_says_other_assay <- md[other, .(Run, BioProject, sample_title)][, annotated := ta[other]]
  md <- md[!other]
  md[, ScientificName := organism_name_cleanup(ScientificName)]
  if (!"COMMENT" %in% names(md)) md[, COMMENT := ""]
  md[, COMMENT := as.character(COMMENT)][is.na(COMMENT), COMMENT := ""]
  hit <- md$ScientificName %in% names(organism_map)
  md[hit, COMMENT := trimws(paste0(COMMENT, " SRA taxon ", ScientificName, "; mapped to species reference"))]
  md[hit, ScientificName := organism_map[ScientificName]]
  md[, study_accession := BioProject][, LIBRARYTYPE := libtype][, KEEP := TRUE][, CHECKED := checked]
  if (!"FRACTION" %in% names(md)) md[, FRACTION := ""]
  md[Run %in% mito_runs & (is.na(FRACTION) | FRACTION == ""), FRACTION := "mito"]
  if (nzchar(comment)) md[, COMMENT := trimws(paste(COMMENT, comment))]
  sheets <- as.data.table(sheets)
  dropped$in_sheet <- md[Run %in% sheets$Run, .(Run, BioProject)]
  md <- md[!(Run %in% sheets$Run)]
  redep <- metadata_biosample_duplicates(md, sheets)
  dropped$biosample_in_other_bioproject <- redep
  md <- md[!(Run %in% redep$Run)]
  if (biosample_labels && "BioSample" %in% names(md))
    md <- metadata_biosample_labels(md, metadata_biosample_attributes(md$BioSample, biosample_cache), vocab)
  md <- metadata_title_labels(md)
  md <- metadata_make_unique(md)
  list(rows = md, dropped = unique(rbindlist(dropped, idcol = "reason", fill = TRUE)))
}

#' Append staged rows to a project's sheet and read them back
#' @param rows data.table from stage_new_runs()$rows
#' @param project project folder name, e.g. "NGS_pipeline"
#' @param base_dir Bio_data folder
#' @return invisible(list(appended, present, duplicated_runs))
#' @export
stage_append <- function(rows, project = "NGS_pipeline", base_dir = "~/livemount/Bio_data") {
  sheet <- modality_read(project, base_dir, TRUE)
  rows <- rows[, intersect(names(sheet), names(rows)), with = FALSE]
  stopifnot(!any(rows$Run %in% sheet$Run))
  modality_append(project, rows, base_dir, TRUE)
  s <- modality_read(project, base_dir, TRUE)
  res <- list(appended = nrow(rows), present = sum(rows$Run %in% s$Run), duplicated_runs = sum(duplicated(s$Run)))
  message("- ", project, " read-back: ", res$present, "/", res$appended, " rows, duplicated runs in sheet: ", res$duplicated_runs)
  if (res$present != res$appended || res$duplicated_runs > 0) stop("sheet read-back failed")
  invisible(res)
}

#' Refresh the FINAL rows of some studies from the live sheet, all other rows unchanged
#'
#' Validates only \code{studies} (pipeline_validate_metadata(only_curated = TRUE)) into a staging file and splices
#' them into config$complete_metadata. Asserts that every row of other studies is byte-identical, keeps a gzipped
#' backup of the old FINAL in \code{backup_dir}, and writes atomically.
#' @param config the mNGSp config
#' @param studies character, study accessions (BioProject column) to refresh
#' @param backup_dir folder for the FINAL backup
#' @param expected_runs optional character, the exact set of runs the refreshed studies must have
#' @return invisible(data.table of the refreshed rows)
#' @export
final_refresh_studies <- function(config, studies, backup_dir, expected_runs = NULL) {
  txt <- function(v) { v <- as.character(v); v[is.na(v)] <- ""; v }
  oldc <- fread(config$complete_metadata, colClasses = "character")
  dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(oldc, file.path(backup_dir, sprintf("FINAL_before_%s.csv.gz", format(Sys.time(), "%Y%m%d_%H%M%S"))))
  s <- read_sheet_safe(config$google_url); staging <- tempfile(fileext = ".csv")
  ok <- pipeline_validate_metadata(s, config, accessions = studies, only_curated = TRUE, output_file = staging)
  if (!isTRUE(ok)) stop("validation failed for the refreshed studies (see messages above)")
  st <- fread(staging, colClasses = "character")[BioProject %in% studies]
  suppressWarnings({ strip_double_quotes(st); strip_double_quotes(oldc) })
  if (!is.null(expected_runs) && !setequal(st$Run, expected_runs))
    stop("refreshed runs differ from expected: missing ", length(setdiff(expected_runs, st$Run)),
         ", extra ", length(setdiff(st$Run, expected_runs)))
  keep_old <- oldc[!(BioProject %in% studies)]; cand <- rbind(keep_old, st, fill = TRUE)
  for (cl in names(keep_old)) stopifnot(identical(txt(cand[[cl]][seq_len(nrow(keep_old))]), txt(keep_old[[cl]])))
  tmp <- paste0(config$complete_metadata, ".tmp"); fwrite(cand, tmp); file.rename(tmp, config$complete_metadata)
  message("- FINAL: ", nrow(oldc), " -> ", nrow(cand), " rows (", length(studies), " studies refreshed, ", nrow(st), " rows)")
  invisible(st)
}

#' NCBI BioSample attributes for a set of BioSample accessions (read-only, cached per accession)
#' @param ids character, BioSample accessions (SAM...)
#' @param cache_dir folder for per-accession rds files
#' @param batch integer, accessions per efetch request
#' @return data.table(BioSample, attribute, value); attribute is the harmonized name when NCBI gives one
#' @export
metadata_biosample_attributes <- function(ids, cache_dir = file.path(tempdir(), "biosample_cache"), batch = 100) {
  ids <- unique(ids[!is.na(ids) & grepl("^SAM", ids)])
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cpath <- function(i) file.path(cache_dir, paste0(i, ".rds"))
  todo <- ids[!file.exists(cpath(ids))]
  for (k in seq(1, length(todo), by = batch)[seq_len(ceiling(length(todo) / batch))]) {
    b <- todo[k:min(k + batch - 1, length(todo))]
    # efetch id= takes numeric UIDs: "SAMEA7425700" would be read as UID 7425700 (another sample). Resolve the
    # accessions with esearch [accn] first, then fetch the UIDs; records are matched back by their accession.
    get <- function(url) { for (a in 1:5) { x <- tryCatch(xml2::read_xml(url), error = function(e) NULL); if (!is.null(x)) return(x); Sys.sleep(3 * a) }; NULL }
    term <- utils::URLencode(paste(paste0(b, "[accn]"), collapse = " OR "), reserved = TRUE)
    es <- get(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&retmax=", 2 * batch, "&term=", term))
    uids <- if (is.null(es)) character() else xml2::xml_text(xml2::xml_find_all(es, "//IdList/Id"))
    if (!length(uids)) { if (is.null(es)) warning("BioSample esearch failed for ", length(b), " accessions"); next }
    doc <- get(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id=", paste(uids, collapse = ",")))
    if (is.null(doc)) { warning("BioSample efetch failed for ", length(b), " accessions"); next }
    for (bs in xml2::xml_find_all(doc, "//BioSample")) {
      if (!xml2::xml_attr(bs, "accession") %in% b) next
      at <- xml2::xml_find_all(bs, "./Attributes/Attribute")
      nm <- xml2::xml_attr(at, "harmonized_name"); raw <- xml2::xml_attr(at, "attribute_name")
      nm[is.na(nm)] <- raw[is.na(nm)]
      saveRDS(data.table(BioSample = xml2::xml_attr(bs, "accession"), attribute = tolower(nm),
                         value = trimws(xml2::xml_text(at))), cpath(xml2::xml_attr(bs, "accession")))
    }
    Sys.sleep(0.4)
  }
  have <- ids[file.exists(cpath(ids))]
  if (!length(have)) return(data.table(BioSample = character(), attribute = character(), value = character()))
  rbindlist(lapply(cpath(have), readRDS))
}

#' Fill empty CELL_LINE / TISSUE from BioSample attributes (author-declared), optionally canonicalised
#'
#' Only empty fields are filled ("cell_line"/"cell line" -> CELL_LINE, "tissue"/"organism_part" -> TISSUE); values
#' like "not applicable" are ignored. With \code{vocab} (\code{vocabulary_load()}), a value that is a literal alias of
#' one vocabulary term is replaced by that term. Filled rows get a COMMENT note: author-declared identity, not
#' verified (misidentified lines such as HeLa derivatives need data checks).
#' @param dt data.table with BioSample, CELL_LINE, TISSUE, COMMENT
#' @param attrs output of metadata_biosample_attributes()
#' @param vocab optional vocabulary table
#' @return dt (copy)
#' @export
metadata_biosample_labels <- function(dt, attrs, vocab = NULL) {
  dt <- copy(as.data.table(dt))
  for (cl in c("CELL_LINE", "TISSUE", "COMMENT")) {
    if (!cl %in% names(dt)) dt[, (cl) := ""]
    set(dt, j = cl, value = { v <- as.character(dt[[cl]]); v[is.na(v)] <- ""; v })
  }
  if (!nrow(attrs)) return(dt)
  junk <- "^(n/?a|na|none|not applicable|not collected|not provided|missing|unknown|-|\\.)$"
  pick <- function(keys) {
    a <- attrs[attribute %in% keys & !grepl(junk, value, ignore.case = TRUE) & nzchar(value)]
    a[, .(value = value[1]), by = BioSample]
  }
  canon <- function(x, column) {
    if (is.null(vocab)) return(x)
    v <- vocab[Column == column & Category != "convertToTissue", .(n = uniqueN(Main), Main = Main[1]), by = .(key = tolower(alias))]
    h <- v[match(tolower(x), key)]
    ifelse(!is.na(h$n) & h$n == 1, h$Main, x)
  }
  note <- character(nrow(dt))
  for (m in list(c("CELL_LINE", "cell_line", "cell line", "cellline"), c("TISSUE", "tissue", "organism_part", "tissue_type"))) {
    p <- pick(m[-1]); i <- which(dt[[m[1]]] == "" & dt$BioSample %in% p$BioSample)
    if (length(i)) {
      set(dt, i, m[1], canon(p$value[match(dt$BioSample[i], p$BioSample)], m[1]))
      note[i] <- paste0(note[i], m[1], " ")
    }
  }
  i <- which(nzchar(note))
  if (length(i)) set(dt, i, "COMMENT", trimws(paste0(dt$COMMENT[i], " ", trimws(note[i]), " from BioSample (author-declared, not verified)")))
  dt[]
}
