#' Controlled vocabulary for sample metadata (RiboSeqOrg vocabularies + local extension)
#'
#' Loads the RiboSeqOrg vocabulary (Content tab: Category, Column, Main Name, All Names) from a local
#' snapshot directory (one csv per tab) or the Google sheet, plus an optional local extension csv with
#' the same columns (terms we use that RiboSeqOrg lacks; add a \code{reason} column).\cr
#' "All Names" entries serve two purposes: literal aliases of a value (normalisation) and regex patterns
#' to find the value in free text such as sample titles (derivation, e.g. \code{^HEK$}, \code{_harr_}).
#' @param dir directory with \code{Content.csv} (and optional \code{Content_hardcoded.csv}), or NULL to
#'  download from \code{url}
#' @param url the RiboSeqOrg vocabulary Google sheet
#' @param extension optional csv path with columns Category, Column, Main Name, All Names (+ reason)
#' @return data.table with Category, Column, Main, alias (one row per alias, Main included)
#' @export
vocabulary_load <- function(dir = NULL,
                            url = "https://docs.google.com/spreadsheets/d/1P9otyhFHY-BTkP8F7mUQS4ukx0f2juIeAGXDrFSPJxY",
                            extension = NULL) {
  read_tab <- function(name) {
    if (!is.null(dir)) {
      f <- file.path(dir, paste0(name, ".csv"))
      return(if (file.exists(f)) fread(f, colClasses = "character") else NULL)
    }
    as.data.table(googlesheets4::range_read(url, sheet = name, col_types = "c", .name_repair = "minimal"))
  }
  tabs <- list(read_tab("Content"), read_tab("Content_hardcoded"))
  if (!is.null(extension) && file.exists(extension)) tabs <- c(tabs, list(fread(extension, colClasses = "character")))
  tabs <- lapply(Filter(Negate(is.null), tabs), function(x) {
    x <- x[, 1:4]; setnames(x, c("Category", "Column", "Main", "All")); x
  })
  v <- rbindlist(tabs)
  v <- v[!is.na(Column) & Column != "" & !is.na(Main) & Main != ""]
  v[, Main := trimws(Main)]
  aliases <- v[, .(alias = vocabulary_split_aliases(All)), by = .(Category, Column, Main)]
  out <- unique(rbind(aliases, unique(v[, .(Category, Column, Main, alias = Main)])))
  out[alias != ""][]
}

#' @noRd
vocabulary_split_aliases <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]            # a Main name can have several vocabulary rows
  if (!length(x)) return(character())
  a <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  a <- gsub('^"|"$', "", a)
  unique(a[nzchar(a)])
}

#' Check sample metadata against the controlled vocabulary
#'
#' Per sample and column: \code{empty}, \code{canonical} (equals a Main name), \code{alias} (a literal alias:
#' proposed = Main), \code{other_column} (the value is a term of a related column, e.g. a tissue written in
#' CELL_LINE), or \code{unknown}. Empty values get a proposal derived from free-text columns when exactly
#' one Main name's patterns match (\code{derived}); TISSUE is also proposed from a canonical CELL_LINE via
#' the vocabulary's convertToTissue category (\code{from_cell_line}).
#' @param dt metadata table (e.g. the RFP sheet)
#' @param vocab output of \code{vocabulary_load()}
#' @param columns metadata columns to check
#' @param text_cols free-text columns used to derive missing values
#' @param related named list: extra vocabulary columns searched for a sheet column
#' @return data.table report: Run, BioProject, column, value, status, proposed, proposed_column, evidence
#' @export
vocabulary_check <- function(dt, vocab,
                             columns = c("CELL_LINE", "TISSUE", "INHIBITOR", "TIMEPOINT", "CONDITION", "REPLICATE",
                                         "FRACTION", "BATCH", "LIBRARYTYPE"),
                             text_cols = c("sample_title", "LibraryName", "SampleName", "sample_source"),
                             related = list(CELL_LINE = "TISSUE", TISSUE = "CELL_LINE")) {
  dt <- as.data.table(dt)
  columns <- intersect(columns, names(dt))
  text_cols <- intersect(text_cols, names(dt))
  vocab <- copy(vocab)[, key := tolower(alias)]
  fields <- lapply(text_cols, function(cl) { x <- trimws(as.character(dt[[cl]])); x[is.na(x)] <- ""; x })
  text <- if (length(fields)) do.call(paste, c(fields, sep = " | ")) else rep("", nrow(dt))
  # candidates per row: each field as is (anchored patterns like ^HEK$ describe a whole field) and "_"-padded
  # with spaces as "_" (underscore-delimited patterns like _harr_ hit at token boundaries)
  candidates <- lapply(seq_len(nrow(dt)), function(k) {
    f <- vapply(fields, `[`, character(1), k); f <- f[nzchar(f)]
    unique(c(f, paste0("_", gsub("[[:space:]]+", "_", f), "_")))
  })
  norm <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }

  rbindlist(lapply(columns, function(cl) {
    v <- norm(dt[[cl]])
    voc <- vocab[Column == cl & Category != "convertToTissue"]
    mains <- unique(voc$Main)
    status <- ifelse(v == "", "empty", ifelse(v %in% mains, "canonical", NA_character_))
    proposed <- rep(NA_character_, length(v)); proposed_col <- rep(NA_character_, length(v)); evidence <- rep(NA_character_, length(v))
    # literal alias (case-insensitive, unique Main)
    amb <- voc[, .(n = uniqueN(Main), Main = Main[1]), by = key]
    i <- which(is.na(status))
    hit <- amb[match(tolower(v[i]), key)]
    ok <- !is.na(hit$n) & hit$n == 1
    status[i[ok]] <- "alias"; proposed[i[ok]] <- hit$Main[ok]; proposed_col[i[ok]] <- cl; evidence[i[ok]] <- "literal alias"
    # value belongs to a related column (e.g. tissue written in CELL_LINE)
    for (rc in related[[cl]]) {
      j <- which(is.na(status))
      if (!length(j)) break
      rvoc <- vocab[Column == rc & Category != "convertToTissue", .(n = uniqueN(Main), Main = Main[1]), by = key]
      h <- rvoc[match(tolower(v[j]), key)]
      ok <- !is.na(h$n) & h$n == 1
      status[j[ok]] <- "other_column"; proposed[j[ok]] <- h$Main[ok]; proposed_col[j[ok]] <- rc
      evidence[j[ok]] <- paste("term of", rc)
    }
    status[is.na(status)] <- "unknown"
    # derive empty values from free text: exactly one Main whose patterns match
    e <- which(status == "empty" & nzchar(text))
    if (length(e) && nrow(voc)) {
      pats <- unique(voc[, .(Main, alias)])
      found <- vapply(e, function(k) {
        m <- unique(pats$Main[vapply(pats$alias, function(p) any(vocabulary_match(p, candidates[[k]])), logical(1))])
        if (length(m) == 1) m else if (length(m) > 1) paste0("AMBIGUOUS:", paste(m, collapse = "|")) else NA_character_
      }, character(1))
      ok <- !is.na(found) & !startsWith(found, "AMBIGUOUS:")
      status[e[ok]] <- "derived"; proposed[e[ok]] <- found[ok]; proposed_col[e[ok]] <- cl
      evidence[e[ok]] <- paste("free text:", substr(text[e[ok]], 1, 120))
      amb_i <- !is.na(found) & startsWith(found, "AMBIGUOUS:")
      evidence[e[amb_i]] <- found[amb_i]
    }
    data.table(Run = as.character(dt$Run), BioProject = as.character(dt$BioProject), column = cl, value = v,
               status = status, proposed = proposed, proposed_column = proposed_col, evidence = evidence)
  }))[, rbind(.SD, vocabulary_tissue_from_cell_line(dt, vocab))[!duplicated(paste(Run, column, status))]]
}

#' Vocabulary pattern match: regex (case-sensitive, as written in the vocabulary), literal on regex error
#' @noRd
vocabulary_match <- function(pattern, text) {
  if (!length(text)) return(logical())
  r <- tryCatch(grepl(pattern, text, perl = TRUE), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(r)) grepl(pattern, text, fixed = TRUE) else r
}

#' TISSUE proposals from canonical CELL_LINE values (vocabulary category convertToTissue)
#' @noRd
vocabulary_tissue_from_cell_line <- function(dt, vocab) {
  if (!all(c("CELL_LINE", "TISSUE") %in% names(dt))) return(NULL)
  conv <- vocab[Category == "convertToTissue"]
  if (!nrow(conv)) return(NULL)
  cl <- trimws(as.character(dt$CELL_LINE)); cl[is.na(cl)] <- ""
  ti <- trimws(as.character(dt$TISSUE)); ti[is.na(ti)] <- ""
  conv <- conv[, .(n = uniqueN(Main), Main = Main[1]), by = .(key = tolower(alias))]
  h <- conv[match(tolower(cl), key)]
  i <- which(ti == "" & cl != "" & !is.na(h$n) & h$n == 1)
  if (!length(i)) return(NULL)
  data.table(Run = as.character(dt$Run[i]), BioProject = as.character(dt$BioProject[i]), column = "TISSUE", value = "",
             status = "from_cell_line", proposed = h$Main[i], proposed_column = "TISSUE",
             evidence = paste("cell line", cl[i]))
}

#' Per-study summary of a vocabulary_check() report
#' @param report output of \code{vocabulary_check()}
#' @return data.table, one row per study and column with counts per status
#' @export
vocabulary_study_summary <- function(report) {
  dcast(report[, .N, by = .(BioProject, column, status)], BioProject + column ~ status, value.var = "N", fill = 0)
}

#' Species whose samples are single cells: CELL_LINE holds the strain, TISSUE is not required
#' @return character vector of scientific names
#' @export
unicellular_species <- function() {
  c("Saccharomyces cerevisiae", "Schizosaccharomyces pombe", "Candida albicans", "Neurospora crassa",
    "Escherichia coli", "Bacillus subtilis", "Mycobacterium tuberculosis", "Mycobacterium smegmatis",
    "Caulobacter crescentus", "Vibrio cholerae", "Staphylococcus aureus", "Pseudomonas aeruginosa",
    "Salmonella enterica", "Streptococcus pneumoniae", "Haloferax volcanii", "Sulfolobus acidocaldarius",
    "Chlamydomonas reinhardtii", "Plasmodium falciparum", "Trypanosoma brucei", "Tetrahymena thermophila",
    "Dictyostelium discoideum", "Toxoplasma gondii", "Zymomonas mobilis", "Synechocystis sp. PCC 6803")
}

#' CELL_LINE / TISSUE convention proposals (see project decision 2026-10-13_cell_line_conventions)
#'
#' Multicellular: primary tissue samples have CELL_LINE "primary" (placeholders like primary_fibroblast,
#' embryonic, fibroblast become "primary" + TISSUE); empty CELL_LINE with a TISSUE -> "primary" (review).
#' Unicellular: CELL_LINE is the strain; missing TISSUE is not an issue.
#' @param dt metadata table with ScientificName, CELL_LINE, TISSUE, Run, BioProject
#' @param unicellular species treated as unicellular
#' @return data.table of proposals in the vocabulary_check() report format
#' @export
metadata_conventions <- function(dt, unicellular = unicellular_species()) {
  dt <- as.data.table(dt)
  norm <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; x }
  cl <- norm(dt$CELL_LINE); ti <- norm(dt$TISSUE); multi <- !dt$ScientificName %in% unicellular
  out <- list()
  add <- function(i, column, value, proposed, status, evidence) {
    if (length(i)) out[[length(out) + 1L]] <<- data.table(Run = as.character(dt$Run[i]), BioProject = as.character(dt$BioProject[i]),
      column = column, value = value[i], status = status, proposed = proposed, proposed_column = column, evidence = evidence)
  }
  placeholder_tissue <- c(primary_fibroblast = "fibroblast", fibroblast = "fibroblast", embryonic = "embryo")
  i <- which(multi & tolower(cl) %in% names(placeholder_tissue))
  add(i, "CELL_LINE", cl, "primary", "convention", "placeholder cell line -> primary")
  j <- i[ti[i] == ""]
  if (length(j)) add(j, "TISSUE", ti, unname(placeholder_tissue[tolower(cl[j])]), "convention", paste("from placeholder", cl[j]))
  add(which(multi & cl == "Primary"), "CELL_LINE", cl, "primary", "convention", "case")
  add(which(multi & cl == "" & ti != ""), "CELL_LINE", cl, "primary", "convention_review",
      "multicellular sample with TISSUE but no CELL_LINE: primary tissue? (confirm not an unrecorded cell line)")
  if (!length(out)) return(data.table())
  rbindlist(out)
}

#' Add our cell line sheet to the vocabulary (our names are canonical)
#'
#' The cell line sheet ("CELL_LINE | TISSUE" pairs with cancer status etc., used by
#' \code{extend_done_samples_metadata()}) defines the cell lines of our data. Each cell line becomes a CELL_LINE
#' Main name and its tissue a convertToTissue mapping. Where a RiboSeqOrg Main denotes the same line under another
#' spelling (our name equals one of its aliases, case-insensitive), our name becomes the Main and the RiboSeqOrg
#' Main and aliases become aliases of ours, so data already using our names is not renamed.
#' @param vocab output of \code{vocabulary_load()}
#' @param cell_lines the cell line sheet (column \code{Tissue...Cell.line} = "CELL_LINE | TISSUE")
#' @return vocab with cell line terms added/re-keyed
#' @export
vocabulary_add_cell_lines <- function(vocab, cell_lines) {
  pair <- as.character(cell_lines[["Tissue...Cell.line"]])
  pair <- pair[!is.na(pair) & grepl("\\|", pair)]
  ours <- data.table(cell = trimws(sub("\\|.*", "", pair)), tissue = trimws(sub(".*\\|", "", pair)))
  ours <- unique(ours[cell != ""])
  v <- copy(vocab)
  for (cl in unique(ours$cell)) {
    same <- unique(v[Column == "CELL_LINE" & Category != "convertToTissue" & tolower(alias) == tolower(cl)]$Main)
    if (length(same) == 1 && same != cl) v[Column == "CELL_LINE" & Main == same, Main := cl]
  }
  tissues <- unique(ours$tissue[ours$tissue != ""])
  add <- rbind(data.table(Category = "Cell Lines (our sheet)", Column = "CELL_LINE", Main = ours$cell, alias = ours$cell),
               data.table(Category = "Tissue (our sheet)", Column = "TISSUE", Main = tissues, alias = tissues),
               data.table(Category = "convertToTissue", Column = "TISSUE", Main = ours$tissue[ours$tissue != ""],
                          alias = ours$cell[ours$tissue != ""]))
  unique(rbind(v, add))
}
