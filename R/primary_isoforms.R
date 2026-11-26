# Primary (canonical) isoform per gene, for subsetting isoform-level outputs (translon tables, tracks) to one transcript
# per gene. Evidence order (method = "auto"): MANE Select > Ensembl canonical > RefSeq Select > expression (optional)
# > longest CDS (coding genes; ties: longest transcript) / longest transcript (non-coding genes).
# Tags are read from the annotation text (GTF or GFF3) next to the TxDb: Ensembl GTFs carry tag "MANE_Select" (human)
# and tag "Ensembl_canonical" (Ensembl >= 104, all genes); NCBI RefSeq GFFs carry tag=MANE Select (human) and
# tag=RefSeq Select. Species without such tags get the longest-CDS rule, which is what Ensembl canonical mostly reduces
# to for coding genes without other evidence.

#' Fast scan of primary-isoform tags in a GTF / GFF3 file
#'
#' Reads only the annotation lines that carry a MANE Select, Ensembl canonical or RefSeq Select tag (grep prefilter,
#' plain R fallback) and returns the transcript identifiers named on those lines.
#' @param annotation path to a GTF or GFF3 file (optionally .gz)
#' @return data.table with columns id (transcript identifier as written: transcript_id, ID, Name or Parent, with the
#' NCBI "rna-" / Ensembl "transcript:" prefixes also stripped), mane, ensembl_canonical, refseq_select (logical).
#' @export
primary_isoform_tags <- function(annotation) {
  stopifnot(is.character(annotation), length(annotation) == 1, file.exists(annotation))
  pattern <- "MANE_Select|MANE Select|Ensembl_canonical|RefSeq Select"
  gz <- grepl("\\.gz$", annotation)
  lines <- if (nzchar(Sys.which(if (gz) "zgrep" else "grep"))) {
    suppressWarnings(system2(if (gz) "zgrep" else "grep", c("-E", shQuote(pattern), shQuote(annotation)),
                             stdout = TRUE, stderr = FALSE))
  } else {
    l <- readLines(annotation); l[grepl(pattern, l)]
  }
  empty <- data.table(id = character(), mane = logical(), ensembl_canonical = logical(), refseq_select = logical())
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) return(empty)
  attr <- sub("^([^\t]*\t){8}", "", lines)
  feat <- sub("^[^\t]*\t[^\t]*\t([^\t]*)\t.*", "\\1", lines)
  gff <- !grepl('transcript_id "', attr[1], fixed = TRUE)
  grab <- function(key) {
    re <- if (gff) paste0("(^|;)", key, "=([^;]+)") else paste0(key, ' "([^"]+)"')
    m <- regmatches(attr, regexec(re, attr))
    vapply(m, function(x) if (length(x)) x[length(x)] else NA_character_, "")
  }
  flags <- data.table(mane = grepl("MANE_Select|MANE Select", attr), ensembl_canonical = grepl("Ensembl_canonical", attr),
                      refseq_select = grepl("RefSeq Select", attr))
  ids <- if (gff) {
    # transcript features name themselves (ID / Name); exon / CDS lines point to their transcript (Parent)
    is_child <- feat %in% c("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon", "UTR")
    list(ifelse(is_child, NA, grab("ID")), ifelse(is_child, NA, grab("Name")), ifelse(is_child, grab("Parent"), NA))
  } else list(grab("transcript_id"))
  res <- rbindlist(lapply(ids, function(id) cbind(data.table(id = id), flags)))
  res <- res[!is.na(id)]
  res <- res[, .(id = unlist(strsplit(id, ",", fixed = TRUE))), by = .(mane, ensembl_canonical, refseq_select)]
  stripped <- copy(res)[, id := sub("^(rna-|transcript:|transcript-)", "", id)]
  res <- rbind(res, stripped)
  res <- res[, .(mane = any(mane), ensembl_canonical = any(ensembl_canonical), refseq_select = any(refseq_select)), by = id]
  setcolorder(res, c("id", "mane", "ensembl_canonical", "refseq_select"))
  res[]
}

#' Transcript table (gene, lengths) from a TxDb, GTF/GFF path or ORFik experiment
#' @param annotation a TxDb, a path to a TxDb (.db) or GTF/GFF file (a ".db" next to it is used when present), or
#' an ORFik experiment
#' @return list(tx = data.table(tx_name, gene_id, tx_len, cds_len), annotation_file = path of the text annotation or NA)
#' @keywords internal
primary_isoform_tx_table <- function(annotation) {
  ann_file <- NA_character_
  if (is(annotation, "experiment")) annotation <- annotation@txdb
  if (is.character(annotation)) {
    stopifnot(length(annotation) == 1, file.exists(annotation))
    if (grepl("\\.db$", annotation)) {
      ann_file <- sub("\\.db$", "", annotation)
      txdb <- AnnotationDbi::loadDb(annotation)
    } else {
      ann_file <- annotation
      txdb <- if (file.exists(paste0(annotation, ".db"))) AnnotationDbi::loadDb(paste0(annotation, ".db")) else {
        if (!requireNamespace("txdbmaker", quietly = TRUE)) stop("Install txdbmaker to build a TxDb from ", annotation)
        suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(annotation)))
      }
    }
    if (!file.exists(ann_file)) ann_file <- NA_character_
  } else if (is(annotation, "TxDb")) {
    txdb <- annotation
  } else stop("annotation must be a TxDb, a .db / GTF / GFF path or an ORFik experiment")
  tl <- setDT(GenomicFeatures::transcriptLengths(txdb, with.cds_len = TRUE))
  tl <- tl[, .(tx_name, gene_id = as.character(gene_id), tx_len, cds_len)]
  tl <- tl[!is.na(tx_name) & tx_name != ""]                        # unnamed transcripts cannot be referenced downstream
  tl[is.na(gene_id) | gene_id == "", gene_id := tx_name]            # unlinked transcripts are their own gene
  list(tx = tl, annotation_file = ann_file)
}

#' Choose one primary (canonical) isoform per gene
#'
#' One transcript per gene, for subsetting isoform-level outputs (translon tables, coverage tracks) to a
#' non-redundant set. Coding genes (any transcript with a CDS) only consider CDS-bearing transcripts; non-coding genes
#' consider all their transcripts. With method "auto" the first available evidence wins:
#' MANE Select > Ensembl canonical > RefSeq Select > expression (only if \code{expression} is given) > longest CDS
#' (coding; ties: longest transcript, then name) / longest transcript (non-coding). A single method ("mane",
#' "ensembl_canonical", "refseq_select", "expression") uses only that evidence and falls back to the longest rule for
#' genes without it (\code{fallback = TRUE}) or drops them (\code{fallback = FALSE}). The rule used for every gene is
#' recorded in column \code{method}.
#'
#' Tags come from the GTF / GFF3 text: Ensembl GTFs (release >= 104) tag every gene's canonical transcript
#' ("Ensembl_canonical") and human MANE Select transcripts ("MANE_Select"); NCBI RefSeq GFFs tag "RefSeq Select" and
#' "MANE Select". Species from other sources (GenBank, custom builds, prokaryotes) use the longest-CDS rule.
#' @param annotation a TxDb, a path to a TxDb (.db) or a GTF / GFF3 (the text annotation is the .db path without
#' ".db"), or an ORFik experiment
#' @param method one of "auto", "mane", "ensembl_canonical", "refseq_select", "longest_cds", "expression"
#' @param expression optional per-transcript expression (e.g. merged RNA or Ribo-seq counts): a named numeric vector
#' (names = transcript names) or a data.frame with columns tx_name and value. Genes whose candidates all have zero /
#' missing expression fall through to the longest rule.
#' @param annotation_file path to the GTF / GFF3 holding the tags, default: derived from \code{annotation}
#' @param fallback logical, for single methods: use the longest rule for genes without evidence (default TRUE)
#' @param cache_file optional csv path: read when it exists (and \code{overwrite = FALSE}), else written
#' @param overwrite logical, recompute even if \code{cache_file} exists
#' @return data.table, one row per gene: gene_id, tx_name, method (evidence used), coding, tx_len, cds_len,
#' n_tx (transcripts of the gene)
#' @export
#' @examples
#' \dontrun{
#' df <- ORFik::read.experiment("all_merged-Homo_sapiens")
#' prim <- primary_isoforms(df)
#' prim[, .N, by = method]
#' }
primary_isoforms <- function(annotation, method = c("auto", "mane", "ensembl_canonical", "refseq_select",
                                                    "longest_cds", "expression"),
                             expression = NULL, annotation_file = NULL, fallback = TRUE,
                             cache_file = NULL, overwrite = FALSE) {
  method <- match.arg(method)
  if (!is.null(cache_file) && file.exists(cache_file) && !overwrite) return(fread(cache_file))
  if (method == "expression" && is.null(expression)) stop("method = 'expression' needs the expression argument")
  tt <- primary_isoform_tx_table(annotation)
  tx <- tt$tx
  if (is.null(annotation_file)) annotation_file <- tt$annotation_file
  tag_methods <- c("mane", "ensembl_canonical", "refseq_select")
  use_tags <- if (method == "auto") tag_methods else intersect(method, tag_methods)
  tx[, (tag_methods) := FALSE]
  if (length(use_tags) > 0) {
    if (is.na(annotation_file) || !file.exists(annotation_file)) {
      warning("No GTF/GFF text annotation found for tags, using fallback rules only")
    } else {
      tags <- primary_isoform_tags(annotation_file)
      for (m in use_tags) tx[, (m) := tx_name %in% tags$id[tags[[m]]]]
    }
  }
  tx[, expr := NA_real_]
  if (!is.null(expression) && method %in% c("auto", "expression")) {
    if (is.numeric(expression)) {
      stopifnot(!is.null(names(expression)))
      tx[, expr := as.numeric(expression[tx_name])]
    } else {
      e <- as.data.table(expression); stopifnot(all(c("tx_name", "value") %in% names(e)))
      tx[, expr := e$value[match(tx_name, e$tx_name)]]
    }
  }
  tx[, n_tx := .N, by = gene_id]
  tx[, coding := any(cds_len > 0), by = gene_id]
  cand <- tx[!coding | cds_len > 0]                                # coding genes: CDS-bearing transcripts only
  # deterministic order inside a gene: longest CDS, longest transcript, name
  setorder(cand, gene_id, -cds_len, -tx_len, tx_name)
  order_rules <- switch(method, auto = c(tag_methods, "expression"), longest_cds = character(), method)
  cand[, `:=`(chosen = NA_character_, selected = FALSE, row = seq_len(.N))]
  # row ids, not transcript names: some annotations reuse a transcript name in several genes
  for (rule in order_rules) {
    pick <- if (rule == "expression") {
      e <- cand[is.na(chosen) & is.finite(expr) & expr > 0]
      e[order(gene_id, -expr)][!duplicated(gene_id)]$row
    } else cand[is.na(chosen) & get(rule)][!duplicated(gene_id)]$row   # cand is sorted: first = longest
    genes <- unique(cand$gene_id[pick])
    cand[gene_id %in% genes & is.na(chosen), chosen := rule]
    cand[pick, selected := TRUE]
  }
  if (method == "longest_cds" || fallback || method == "auto") {
    fb <- cand[is.na(chosen)][!duplicated(gene_id)]$row
    cand[fb, `:=`(chosen = ifelse(coding, "longest_cds", "longest_tx"), selected = TRUE)]
  }
  res <- cand[selected %in% TRUE, .(gene_id, tx_name, method = chosen, coding, tx_len, cds_len, n_tx)]
  stopifnot(!anyDuplicated(res$gene_id))
  setorder(res, gene_id)
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    fwrite(res, cache_file)
  }
  res[]
}

#' Subset a translon table (and its ranges) to primary isoforms
#'
#' Keeps rows whose transcript is a primary isoform (see \code{\link{primary_isoforms}}), as the published human
#' translon table does (transcripts of the canonical set only).
#' @param table data.table of translons with a transcript column
#' @param primary character vector of primary transcript names, or the output of \code{primary_isoforms}
#' @param grl optional GRangesList parallel to \code{table} (same order and length)
#' @param tx_col transcript column of \code{table}, default "ensembl_tx_name"
#' @return list(table, grl) subset
#' @export
subset_primary_isoforms <- function(table, primary, grl = NULL, tx_col = "ensembl_tx_name") {
  if (is.data.frame(primary)) primary <- primary$tx_name
  stopifnot(tx_col %in% names(table))
  if (!is.null(grl)) stopifnot(length(grl) == nrow(table))
  keep <- table[[tx_col]] %in% primary
  list(table = table[keep], grl = if (is.null(grl)) NULL else grl[keep])
}
