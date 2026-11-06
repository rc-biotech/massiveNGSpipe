#' Build a gene annotation for a poorly annotated genome
#'
#' A standard three-layer annotation for species whose public annotation is missing or unusable (e.g. Saccharomyces jurei,
#' GCA_900290405.1: genes without CDS/exon features), in the spirit of Ensembl's projection pipelines:
#' \enumerate{
#'  \item homology: every protein of a reference species (e.g. S. cerevisiae) is aligned to the target genome (tblastn); HSPs are chained
#'   into exons, introns are placed at GT..AG sites that keep the frame (scored by donor, acceptor and branch-point motifs), the start is
#'   snapped to the aligned or nearest in-frame ATG and the end to the first in-frame stop. A short first exon that tblastn cannot see
#'   (e.g. ribosomal protein genes) is searched upstream. One model per locus (best bitscore); orthologue = the projected reference gene.
#'  \item sequence ORFs fill loci without a homology model (\code{annotation_orf_fill()}), and ncRNA genes are placed by blastn.
#'  \item data: RNA-seq junctions and coverage refine introns and UTRs (\code{annotation_rnaseq_refine()}); Ribo-seq P-sites refine the
#'   start codon (\code{annotation_riboseq_refine()}).
#' }
#' \code{annotation_write_gtf()} writes an Ensembl-style GTF that \code{ORFik::makeTxdbFromGenome()} and STAR accept.
#' BLAST+ (tblastn, blastn, blastp, makeblastdb) must be available in \code{blast_bin}.
#' @name annotation_build
NULL

ann_blast <- function(blast_bin, prog, args) {
  st <- suppressWarnings(system2(file.path(blast_bin, prog), args, stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(st, "status"))) stop(prog, " failed: ", paste(head(st, 5), collapse = "\n"))
  invisible(st)
}
ann_fmt6 <- c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "bitscore")
ann_read6 <- function(f) {
  if (!file.exists(f) || file.size(f) == 0) return(setnames(data.table(matrix(nrow = 0, ncol = 12)), ann_fmt6))
  fread(f, header = FALSE, col.names = ann_fmt6)
}

#' Reference proteins, one per gene (longest CDS), from a GTF and genome
#' @param ref_gtf reference GTF path
#' @param ref_genome reference genome fasta path
#' @param exclude_seqnames seqnames to skip (default the mitochondrion; it uses another genetic code)
#' @return list(protein = AAStringSet named by transcript, info = data.table transcript_id, gene_id, gene_name, exons, aa)
#' @export
annotation_reference_proteins <- function(ref_gtf, ref_genome, exclude_seqnames = c("Mito", "MT", "chrM", "mitochondrion")) {
  g <- rtracklayer::import(ref_gtf)
  cds <- g[g$type == "CDS" & !(as.character(seqnames(g)) %in% exclude_seqnames)]
  fa <- readDNAStringSet(ref_genome); names(fa) <- sub(" .*", "", names(fa))
  cl <- split(cds, cds$transcript_id)
  cl <- endoapply(cl, function(x) x[order(start(x), decreasing = as.character(strand(x))[1] == "-")])
  seqs <- GenomicFeatures::extractTranscriptSeqs(fa, cl)
  # Ensembl GTF CDS excludes the stop codon; translate what is there and drop a trailing stop if present
  prot <- as.character(translate(seqs[width(seqs) %% 3 == 0], if.fuzzy.codon = "X"))
  prot <- sub("\\*$", "", prot)
  info <- unique(as.data.table(mcols(cds))[, .(transcript_id, gene_id, gene_name)])[match(names(prot), transcript_id)]
  info[, `:=`(exons = lengths(cl)[transcript_id], aa = nchar(prot))]
  info[is.na(gene_name) | gene_name == "", gene_name := gene_id]
  ok <- !grepl("\\*", prot) & nchar(prot) >= 30
  info <- info[ok]; prot <- prot[ok]
  keep <- info[order(gene_id, -aa)][!duplicated(gene_id)]$transcript_id
  list(protein = AAStringSet(prot[keep]), info = info[match(keep, transcript_id)])
}

# Chain tblastn HSPs of one query on one contig/strand: best collinear chain by summed bitscore (dynamic programming)
ann_chain <- function(h, max_intron) {
  h <- h[order(qstart)]
  n <- nrow(h); if (n == 1) return(h)
  plus <- h$send[1] >= h$sstart[1]
  gs <- if (plus) h$sstart else -h$sstart; ge <- if (plus) h$send else -h$send      # genome position increasing with transcript
  best <- h$bitscore; prev <- rep(NA_integer_, n)
  for (j in seq_len(n)[-1]) for (i in seq_len(j - 1)) {
    gap <- gs[j] - ge[i]
    if (h$qstart[j] >= h$qend[i] - 10 && gap >= -30 && gap <= max_intron && best[i] + h$bitscore[j] > best[j]) {
      best[j] <- best[i] + h$bitscore[j]; prev[j] <- i }
  }
  k <- which.max(best); idx <- k
  while (!is.na(prev[k])) { k <- prev[k]; idx <- c(k, idx) }
  h[idx]
}

ann_motif_score <- function(s, motif) sum(strsplit(s, "")[[1]] == strsplit(motif, "")[[1]])

# Build one CDS model in local (+) coordinates of region sequence R from chained HSPs (local nt coords, query aa coords).
# Returns list(exons = IRanges incl. stop codon, flags) or NULL.
ann_build_model <- function(R, hs, qlen, ref_exon1_aa = NA_integer_, max_intron = 2000L, max_shift = 30L) {
  Rs <- as.character(R); nR <- nchar(Rs)
  sub3 <- function(p) substr(Rs, p, p + 2L)
  stops <- c("TAA", "TAG", "TGA")
  ex_s <- hs$a; ex_e <- hs$b; flags <- character()
  # merge HSPs separated by a small in-frame gap, drop the model on a frame shift between them
  if (nrow(hs) > 1) {
    keep <- rep(TRUE, nrow(hs))
    for (i in seq_len(nrow(hs) - 1)) {
      gap <- ex_s[i + 1] - ex_e[i] - 1
      if (gap < 40) {
        if ((ex_s[i + 1] - ex_s[i]) %% 3 != 0) return(NULL)                          # frame shift: pseudogene / assembly error
        ex_s[i + 1] <- ex_s[i]; keep[i] <- FALSE
      }
    }
    ex_s <- ex_s[keep]; ex_e <- ex_e[keep]; hs <- hs[keep]
  }
  # introns between chained exons: donor GT near the end of exon i, acceptor AG near the start of exon i+1, same frame shift mod 3
  if (length(ex_s) > 1) for (i in seq_len(length(ex_s) - 1)) {
    best <- NULL; bs <- -Inf
    for (d in max(ex_s[i] + 3L, ex_e[i] - max_shift):min(ex_s[i + 1] - 40L, ex_e[i] + max_shift + 1L)) {
      if (substr(Rs, d, d + 1L) != "GT") next
      for (a in max(d + 40L, ex_s[i + 1] - max_shift - 1L):min(ex_e[i + 1] - 3L, ex_s[i + 1] + max_shift)) {
        if (substr(Rs, a - 1L, a) != "AG") next
        if (((d - 1L - ex_e[i]) - (a + 1L - ex_s[i + 1])) %% 3L != 0L) next
        intr <- substr(Rs, d, a)
        sc <- ann_motif_score(substr(Rs, d, d + 5L), "GTATGT") + 2 * grepl("TACTAAC", intr, fixed = TRUE) +
          ann_motif_score(substr(Rs, a - 2L, a), "YAG" |> sub(pattern = "Y", replacement = substr(Rs, a - 2L, a - 2L))) -
          (abs(d - 1L - ex_e[i]) + abs(a + 1L - ex_s[i + 1])) / 15
        if (sc > bs) { bs <- sc; best <- c(d, a) }
      }
    }
    if (is.null(best)) return(NULL)
    ex_e[i] <- best[1] - 1L; ex_s[i + 1] <- best[2] + 1L
  }
  # start codon: aligned Met, else nearest in-frame ATG upstream (before an in-frame stop) or within the first 15 codons
  s1 <- ex_s[1]; q1 <- hs$qstart[1]
  upstream_atg <- function(pos, max_codons) { p <- pos; found <- NA_integer_
    for (k in seq_len(max_codons)) { p <- p - 3L; if (p < 1L) break; cd <- sub3(p); if (cd %in% stops) break; if (cd == "ATG") found <- p }
    found }
  if (sub3(s1) == "ATG" && q1 <= 2) {
    start <- s1
  } else {
    up <- upstream_atg(s1, q1 + 20L)
    exon1 <- NULL
    if (!is.na(ref_exon1_aa) && q1 > 1 && q1 <= 60) {
      # short first exon not seen by tblastn: ATG + ~(q1 - 1) codons + GT .. AG just before the first aligned exon
      L1 <- (q1 - 1L) * 3L; bs <- -Inf
      for (a in max(1L, s1 - max_shift - 1L):(s1 + max_shift)) {
        if (a < 2 || substr(Rs, a - 1L, a) != "AG") next
        for (d in max(1L, a - max_intron):(a - 40L)) {
          if (substr(Rs, d, d + 1L) != "GT") next
          for (l1 in seq(max(3L, L1 - 9L), L1 + 9L)) {
            st <- d - l1; if (st < 1L || sub3(st) != "ATG") next
            if (((d - st) + (ex_e[1] - a)) %% 3L != 0L) next
            cds_try <- paste0(substr(Rs, st, d - 1L), substr(Rs, a + 1L, ex_e[1]))
            aa <- as.character(translate(DNAString(substr(cds_try, 1, 3 * (nchar(cds_try) %/% 3))), if.fuzzy.codon = "X"))
            if (grepl("*", aa, fixed = TRUE)) next
            sc <- ann_motif_score(substr(Rs, d, d + 5L), "GTATGT") + 2 * grepl("TACTAAC", substr(Rs, d, a), fixed = TRUE) - abs(l1 - L1) / 3
            if (sc > bs) { bs <- sc; exon1 <- c(st, d - 1L, a + 1L) }
          }
        }
      }
    }
    if (!is.null(exon1)) {
      ex_s <- c(exon1[1], exon1[3], ex_s[-1]); ex_e <- c(exon1[2], ex_e)
      start <- exon1[1]; flags <- c(flags, "short_exon1")
    } else if (!is.na(up)) {
      start <- up
    } else {
      dn <- NA_integer_; p <- s1
      for (k in 0:15) { pp <- p + 3L * k; if (pp + 2L > ex_e[1]) break; if (sub3(pp) == "ATG") { dn <- pp; break } }
      if (is.na(dn)) return(NULL)
      start <- dn; flags <- c(flags, "start_downstream")
    }
    ex_s[1] <- start
  }
  # stop codon: first in-frame stop from the start through the spliced model, extending the last exon
  last <- length(ex_s)
  spliced_pos <- unlist(lapply(seq_len(last), function(i) ex_s[i]:ex_e[i]))
  ph <- (length(spliced_pos)) %% 3L
  p <- ex_e[last] - ph + 1L                                     # first nt of the next codon in frame after the last full codon
  if (ph == 0L) p <- ex_e[last] + 1L
  cod <- seq(1L, length(spliced_pos) - 2L, by = 3L)
  internal <- which(vapply(cod, function(k) paste(substr(Rs, spliced_pos[k], spliced_pos[k]), substr(Rs, spliced_pos[k + 1], spliced_pos[k + 1]),
                                                   substr(Rs, spliced_pos[k + 2], spliced_pos[k + 2]), sep = "") %in% stops, logical(1)))
  if (length(internal)) {
    # stop inside the aligned part: end there if the protein is still >= 50% of the query, else reject
    k <- cod[internal[1]]; endpos <- spliced_pos[k + 2]
    if ((k - 1) / 3 < 0.5 * qlen) return(NULL)
    keepx <- ex_s <= endpos; ex_s <- ex_s[keepx]; ex_e <- pmin(ex_e[keepx], endpos)
    flags <- c(flags, "early_stop")
  } else {
    repeat {
      if (p + 2L > nR) return(NULL)
      if (sub3(p) %in% stops) { ex_e[last] <- p + 2L; break }
      p <- p + 3L
      if (p - ex_e[last] > 3L * (qlen + 100L)) return(NULL)
    }
  }
  len <- sum(ex_e - ex_s + 1L)
  if (len %% 3L != 0L) return(NULL)
  list(exons = IRanges(ex_s, ex_e), flags = paste(flags, collapse = ";"))
}

#' Project reference proteins onto a target genome (tblastn, chaining, splice-site refinement)
#' @param target_genome target genome fasta
#' @param ref list from \code{annotation_reference_proteins()}
#' @param work_dir working directory for BLAST files
#' @param blast_bin directory with BLAST+ binaries
#' @param threads BLAST threads
#' @param max_intron maximum intron length (nt)
#' @param min_coverage minimum fraction of the reference protein covered by the chained HSPs
#' @return data.table, one row per exon (model_id, seqnames, start, end, strand, exon_rank) with model columns
#'  (ref_gene, ref_tx, ref_name, bitscore, qcov, pident, flags)
#' @export
annotation_homology_project <- function(target_genome, ref, work_dir, blast_bin, threads = 8, max_intron = 2000L,
                                        min_coverage = 0.5) {
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  fa <- readDNAStringSet(target_genome); names(fa) <- sub(" .*", "", names(fa))
  gfa <- file.path(work_dir, "target_genome.fa"); writeXStringSet(fa, gfa)
  pfa <- file.path(work_dir, "ref_proteins.fa"); writeXStringSet(ref$protein, pfa)
  db <- file.path(work_dir, "target_genome_db")
  ann_blast(blast_bin, "makeblastdb", c("-in", gfa, "-dbtype", "nucl", "-out", db))
  out <- file.path(work_dir, "tblastn_ref_vs_target.tsv")
  if (!file.exists(out) || file.size(out) == 0)
    ann_blast(blast_bin, "tblastn", c("-query", pfa, "-db", db, "-evalue", "1e-5", "-outfmt", shQuote(paste("6", paste(ann_fmt6, collapse = " "))),
                                      "-num_threads", threads, "-max_target_seqs", 5, "-seg", "no", "-out", out))
  h <- ann_read6(out)
  h[, strand := fifelse(send >= sstart, "+", "-")]
  ch <- h[, ann_chain(.SD, max_intron), by = .(qseqid, sseqid, strand)]
  sc <- ch[, .(bits = sum(bitscore), qcov = sum(qend - qstart + 1) / qlen[1]), by = .(qseqid, sseqid, strand)]
  sc <- sc[qcov >= min_coverage][order(qseqid, -bits)][!duplicated(qseqid)]          # best locus per reference protein
  ch <- ch[sc, on = .(qseqid, sseqid, strand)]
  info <- ref$info
  models <- list()
  for (q in unique(ch$qseqid)) {
    x <- ch[qseqid == q]; chr <- x$sseqid[1]; st <- x$strand[1]; L <- width(fa[chr])
    lo <- max(1L, min(x$sstart, x$send) - max_intron - 300L); hi <- min(L, max(x$sstart, x$send) + 3000L)
    R <- subseq(fa[[chr]], lo, hi); if (st == "-") R <- reverseComplement(R)
    loc <- function(p) if (st == "+") p - lo + 1L else hi - p + 1L
    hs <- x[, .(a = pmin(loc(sstart), loc(send)), b = pmax(loc(sstart), loc(send)), qstart, qend)][order(a)]
    hs[, b := b - ((b - a + 1L) %% 3L)]
    ex1 <- NA_integer_
    if (info[transcript_id == q, exons] > 1 && hs$qstart[1] > 1) ex1 <- hs$qstart[1] - 1L
    m <- tryCatch(ann_build_model(R, hs, qlen = x$qlen[1], ref_exon1_aa = ex1, max_intron = max_intron), error = function(e) NULL)
    if (is.null(m)) next
    e <- m$exons
    gs <- if (st == "+") start(e) + lo - 1L else hi - end(e) + 1L
    ge <- if (st == "+") end(e) + lo - 1L else hi - start(e) + 1L
    models[[q]] <- data.table(ref_tx = q, seqnames = chr, start = gs, end = ge, strand = st, exon_rank = seq_along(gs),
                              bitscore = sc[qseqid == q, bits], qcov = sc[qseqid == q, qcov], pident = x[, weighted.mean(pident, length)],
                              flags = m$flags)
  }
  mm <- rbindlist(models)
  if (!nrow(mm)) return(mm)
  mm <- merge(mm, info[, .(ref_tx = transcript_id, ref_gene = gene_id, ref_name = gene_name, ref_exons = exons)], by = "ref_tx")
  # one model per locus: models overlapping a better-scoring model on the same strand by > 30% of their length are dropped
  span <- mm[, .(seqnames = seqnames[1], strand = strand[1], start = min(start), end = max(end), bitscore = bitscore[1]), by = ref_tx][order(-bitscore)]
  gr <- GRanges(span$seqnames, IRanges(span$start, span$end), span$strand); keep <- logical(length(gr)); acc <- GRanges()
  for (i in seq_along(gr)) {
    ov <- findOverlaps(gr[i], acc)
    frac <- if (length(ov)) max(width(pintersect(rep(gr[i], length(ov)), acc[subjectHits(ov)]))) / width(gr[i]) else 0
    if (frac <= 0.3) { keep[i] <- TRUE; acc <- c(acc, gr[i]) }
  }
  mm <- mm[ref_tx %in% span$ref_tx[keep]]
  mm[]
}

#' Validate CDS models: ATG start, stop end, no internal stop, length multiple of 3
#' @param models data.table from \code{annotation_homology_project()} (one row per exon, with model id column \code{id})
#' @param genome DNAStringSet of the target genome
#' @return logical per model id (named)
#' @export
annotation_validate_cds <- function(models, genome) {
  gl <- split(GRanges(models$seqnames, IRanges(models$start, models$end), models$strand), models$id)
  gl <- endoapply(gl, function(x) x[order(start(x), decreasing = as.character(strand(x))[1] == "-")])
  s <- as.character(GenomicFeatures::extractTranscriptSeqs(genome, gl))
  ok <- nchar(s) %% 3 == 0 & substr(s, 1, 3) == "ATG" & substr(s, nchar(s) - 2, nchar(s)) %in% c("TAA", "TAG", "TGA")
  aa <- as.character(translate(DNAStringSet(s[ok]), if.fuzzy.codon = "X"))
  ok[ok] <- !grepl("*", sub("\\*$", "", aa), fixed = TRUE)
  setNames(ok, names(gl))
}

#' Sequence ORFs (longest ATG..stop per stop) in loci without a model
#' @param genome DNAStringSet
#' @param occupied GRanges of existing genes (models and ncRNA)
#' @param min_codons minimum length in codons
#' @param max_overlap maximum overlapped fraction with occupied or longer ORFs (either strand)
#' @return GRanges of ORFs (start codon to stop codon inclusive)
#' @export
annotation_orf_fill <- function(genome, occupied, min_codons = 100L, max_overlap = 0.4) {
  one <- function(strand) {
    s <- if (strand == "+") genome else reverseComplement(genome)
    o <- findORFs(s, startCodon = "ATG", longestORF = TRUE, minimumLength = min_codons - 1L)
    unlist(GRangesList(lapply(seq_along(genome), function(i) { ch <- names(genome)[i]                # findORFs names results by index
      r <- o[[as.character(i)]]; if (is.null(r) || !length(r)) return(GRanges())
      if (strand == "-") { L <- width(genome[ch]); r <- IRanges(L - end(r) + 1L, L - start(r) + 1L) }
      GRanges(ch, r, strand) })))
  }
  orf <- c(one("+"), one("-")); orf <- orf[order(-width(orf))]
  acc <- granges(occupied); mcols(acc) <- NULL; keep <- logical(length(orf))
  for (i in seq_along(orf)) {
    ov <- findOverlaps(orf[i], acc, ignore.strand = TRUE)
    frac <- if (length(ov)) sum(width(reduce(pintersect(rep(orf[i], length(ov)), acc[subjectHits(ov)], ignore.strand = TRUE)))) / width(orf[i]) else 0
    if (frac <= max_overlap) { keep[i] <- TRUE; acc <- c(acc, granges(orf[i])) }
  }
  sort(orf[keep])
}

#' Place ncRNA genes of a reference species by blastn
#' @param target_genome target genome fasta
#' @param ref_gtf reference GTF
#' @param ref_genome reference genome fasta
#' @param work_dir working directory
#' @param blast_bin BLAST+ directory
#' @param biotypes gene biotypes to place
#' @param min_cov,min_pident hit filters (fraction of query, percent identity)
#' @param threads BLAST threads
#' @return GRanges with ref_gene, gene_name, gene_biotype
#' @export
annotation_ncrna_blastn <- function(target_genome, ref_gtf, ref_genome, work_dir, blast_bin,
                                    biotypes = c("rRNA", "tRNA", "snoRNA", "snRNA", "ncRNA"), min_cov = 0.8, min_pident = 75, threads = 8) {
  g <- rtracklayer::import(ref_gtf); genes <- g[g$type == "gene" & g$gene_biotype %in% biotypes & !(as.character(seqnames(g)) %in% c("Mito", "MT", "chrM"))]
  rf <- readDNAStringSet(ref_genome); names(rf) <- sub(" .*", "", names(rf))
  q <- getSeq(rf, genes); names(q) <- genes$gene_id; qf <- file.path(work_dir, "ref_ncrna.fa"); writeXStringSet(q, qf)
  db <- file.path(work_dir, "target_genome_db")
  if (!file.exists(paste0(db, ".nsq")) && !file.exists(paste0(db, ".ndb"))) {
    fa <- readDNAStringSet(target_genome); names(fa) <- sub(" .*", "", names(fa)); gfa <- file.path(work_dir, "target_genome.fa"); writeXStringSet(fa, gfa)
    ann_blast(blast_bin, "makeblastdb", c("-in", gfa, "-dbtype", "nucl", "-out", db))
  }
  out <- file.path(work_dir, "blastn_ncrna_vs_target.tsv")
  ann_blast(blast_bin, "blastn", c("-task", "blastn", "-query", qf, "-db", db, "-evalue", "1e-10", "-outfmt", shQuote(paste("6", paste(ann_fmt6, collapse = " "))),
                                   "-num_threads", threads, "-max_target_seqs", 50, "-out", out))
  h <- ann_read6(out)[length / qlen >= min_cov & pident >= min_pident][order(-bitscore)]
  if (!nrow(h)) return(GRanges())
  gr <- GRanges(h$sseqid, IRanges(pmin(h$sstart, h$send), pmax(h$sstart, h$send)), fifelse(h$send >= h$sstart, "+", "-"), ref_gene = h$qseqid)
  keep <- logical(length(gr)); acc <- GRanges()
  for (i in seq_along(gr)) if (!overlapsAny(gr[i], acc)) { keep[i] <- TRUE; acc <- c(acc, gr[i]) }
  gr <- sort(gr[keep]); m <- match(gr$ref_gene, genes$gene_id)
  gr$gene_biotype <- genes$gene_biotype[m]; gr$gene_name <- ifelse(is.na(genes$gene_name[m]) | genes$gene_name[m] == "", genes$gene_id[m], genes$gene_name[m])
  gr
}

#' Write an Ensembl-style GTF from CDS models, ncRNA genes and UTRs
#' @param cds data.table one row per CDS exon incl. stop codon: id, seqnames, start, end, strand, gene_name, gene_biotype, source_tag
#' @param ncrna GRanges with id, gene_name, gene_biotype
#' @param genome DNAStringSet (contig lengths)
#' @param out_gtf output path
#' @param utr5,utr3 named integer vectors (by id) of UTR lengths, or a single default (pseudo-UTRs)
#' @param header character lines for the header
#' @return out_gtf, invisibly
#' @export
annotation_write_gtf <- function(cds, ncrna, genome, out_gtf, utr5 = 100L, utr3 = 100L, header = character()) {
  L <- setNames(width(genome), names(genome))
  ids <- unique(cds$id)
  u5 <- if (length(utr5) == 1) setNames(rep(as.integer(utr5), length(ids)), ids) else utr5[ids]
  u3 <- if (length(utr3) == 1) setNames(rep(as.integer(utr3), length(ids)), ids) else utr3[ids]
  u5[is.na(u5)] <- 0L; u3[is.na(u3)] <- 0L
  rows <- list()
  for (gid in ids) {
    x <- cds[id == gid][order(start)]; st <- x$strand[1]; ch <- x$seqnames[1]
    attr <- sprintf('gene_id "%s"; transcript_id "%s_t1"; gene_name "%s"; gene_biotype "%s"; transcript_biotype "%s"; annotation_source "%s";',
                    gid, gid, x$gene_name[1], x$gene_biotype[1], x$gene_biotype[1], x$source_tag[1])
    ex <- x[, .(start, end)]
    if (st == "+") { ex$start[1] <- max(1L, ex$start[1] - u5[gid]); ex$end[nrow(ex)] <- min(L[ch], ex$end[nrow(ex)] + u3[gid]) }
    else { ex$start[1] <- max(1L, ex$start[1] - u3[gid]); ex$end[nrow(ex)] <- min(L[ch], ex$end[nrow(ex)] + u5[gid]) }
    # CDS without the stop codon; start/stop codon features (stop codon may not be split in yeast-like genomes; assumed contiguous)
    cs <- copy(x[, .(start, end)])
    if (st == "+") { stop_c <- c(cs$end[nrow(cs)] - 2L, cs$end[nrow(cs)]); cs$end[nrow(cs)] <- cs$end[nrow(cs)] - 3L; start_c <- c(cs$start[1], cs$start[1] + 2L) }
    else { stop_c <- c(cs$start[1], cs$start[1] + 2L); cs$start[1] <- cs$start[1] + 3L; start_c <- c(cs$end[nrow(cs)] - 2L, cs$end[nrow(cs)]) }
    cs <- cs[end >= start]
    # CDS frame (phase) along the transcript
    ord <- if (st == "+") seq_len(nrow(cs)) else rev(seq_len(nrow(cs)))
    ph <- integer(nrow(cs)); cum <- 0L
    for (k in ord) { ph[k] <- (3L - cum %% 3L) %% 3L; cum <- cum + cs$end[k] - cs$start[k] + 1L }
    mk <- function(type, s, e, fr = ".") data.table(ch, "massiveNGSpipe_annotation", type, s, e, ".", st, fr, attr)
    rows[[gid]] <- rbind(mk("gene", min(ex$start), max(ex$end)), mk("transcript", min(ex$start), max(ex$end)),
                         mk("exon", ex$start, ex$end), mk("CDS", cs$start, cs$end, as.character(ph)),
                         mk("start_codon", start_c[1], start_c[2], "0"), mk("stop_codon", stop_c[1], stop_c[2], "0"))
  }
  for (i in seq_along(ncrna)) {
    attr <- sprintf('gene_id "%s"; transcript_id "%s_t1"; gene_name "%s"; gene_biotype "%s"; transcript_biotype "%s"; annotation_source "blastn";',
                    ncrna$id[i], ncrna$id[i], ncrna$gene_name[i], ncrna$gene_biotype[i], ncrna$gene_biotype[i])
    mk <- function(type) data.table(as.character(seqnames(ncrna[i])), "massiveNGSpipe_annotation", type, start(ncrna[i]), end(ncrna[i]), ".",
                                    as.character(strand(ncrna[i])), ".", attr)
    rows[[ncrna$id[i]]] <- rbind(mk("gene"), mk("transcript"), mk("exon"))
  }
  g <- rbindlist(rows, use.names = FALSE); setnames(g, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
  g[, gid := sub('gene_id "([^"]+)".*', "\\1", attribute)]
  fo <- c(gene = 1, transcript = 2, exon = 3, CDS = 4, start_codon = 5, stop_codon = 6)
  g[, gstart := min(start), by = gid]
  g <- g[order(match(seqname, names(genome)), gstart, gid, fo[feature], start)][, c("gid", "gstart") := NULL]
  writeLines(header, out_gtf)
  fwrite(g, out_gtf, sep = "\t", quote = FALSE, col.names = FALSE, append = TRUE)
  invisible(out_gtf)
}

#' Valid CDS spans from an existing (partial) annotation
#'
#' Some submitted annotations give coding genes only as 'gene' features (e.g. S. jurei GCA_900290405.1: 5,170 of the 5,672 loci of
#' >= 300 nt are exact ATG..stop spans, the rest are intron-gene spans, transposons and microsatellite 'genes'). This keeps CDS features
#' when present, else gene spans, and returns those that are complete single-exon CDSs (ATG, stop, frame, no internal stop).
#' @param gff existing GFF/GTF path
#' @param genome DNAStringSet of the target genome
#' @param min_codons minimum CDS length in codons
#' @return GRanges with id and valid (only valid spans are returned; \code{attr(, "n_input")} gives the number tested)
#' @export
annotation_existing_cds <- function(gff, genome, min_codons = 50L) {
  g <- rtracklayer::import(gff)
  x <- if (any(g$type == "CDS")) g[g$type == "CDS"] else g[g$type == "gene"]
  if (any(g$type == "CDS")) {                                   # multi-exon CDS: keep single-segment ones only (introns come from homology/RNA)
    pid <- if (!is.null(x$Parent)) vapply(x$Parent, function(p) if (length(p)) p[1] else NA_character_, "") else x$ID
    x <- x[!(pid %in% pid[duplicated(pid)])]
  }
  s <- as.character(getSeq(genome, x)); w <- width(x)
  ok <- w %% 3 == 0 & w >= 3 * (min_codons + 1) & substr(s, 1, 3) == "ATG" & substr(s, nchar(s) - 2, nchar(s)) %in% c("TAA", "TAG", "TGA")
  ok[ok] <- !grepl("*", sub("\\*$", "", as.character(translate(DNAStringSet(s[ok]), if.fuzzy.codon = "X"))), fixed = TRUE)
  out <- granges(x[ok]); out$id <- if (!is.null(x$ID)) x$ID[ok] else paste0("existing_", which(ok))
  attr(out, "n_input") <- length(x)
  out
}

#' Merge homology models with an existing CDS annotation
#'
#' A homology model is 'confirmed' when an existing CDS ends at the same stop codon on the same strand; a differing existing start is
#' kept as \code{alt_start} for the Ribo-seq start-codon step. Existing CDSs that do not overlap a homology model (> 30% of their length,
#' either strand) are added as single-exon models (source 'existing'); those overlapping with a different stop are reported, not used.
#' @param hm data.table of homology CDS exons (id, seqnames, start, end, strand, ...) incl. the stop codon
#' @param existing GRanges from \code{annotation_existing_cds()}
#' @return list(models = data.table of exons with evidence and alt_start, report = data.table per existing CDS)
#' @export
annotation_merge_evidence <- function(hm, existing) {
  span <- hm[, .(seqnames = seqnames[1], strand = strand[1], start = min(start), end = max(end)), by = id]
  span[, `:=`(stop_end = fifelse(strand == "+", end, start), cds_start = fifelse(strand == "+", start, end))]
  ex <- data.table(eid = existing$id, seqnames = as.character(seqnames(existing)), strand = as.character(strand(existing)),
                   start = start(existing), end = end(existing))
  ex[, `:=`(stop_end = fifelse(strand == "+", end, start), ex_start = fifelse(strand == "+", start, end))]
  same_stop <- merge(span[, .(id, seqnames, strand, stop_end, cds_start)], ex[, .(eid, seqnames, strand, stop_end, ex_start)],
                     by = c("seqnames", "strand", "stop_end"))
  hm <- copy(hm)
  hm[, evidence := fifelse(id %in% same_stop$id, "homology+existing", "homology")]
  alt <- same_stop[ex_start != cds_start, .(id, alt_start = ex_start)][!duplicated(id)]
  hm <- merge(hm, alt, by = "id", all.x = TRUE)
  sg <- GRanges(span$seqnames, IRanges(span$start, span$end), span$strand)
  eg <- GRanges(ex$seqnames, IRanges(ex$start, ex$end), ex$strand)
  ov <- findOverlaps(eg, sg, ignore.strand = TRUE)
  frac <- numeric(length(eg))
  if (length(ov)) {
    w <- width(pintersect(eg[queryHits(ov)], sg[subjectHits(ov)], ignore.strand = TRUE))
    frac <- pmax(frac, tapply(w, factor(queryHits(ov), levels = seq_along(eg)), max, default = 0) / width(eg))
  }
  ex[, class := fifelse(eid %in% same_stop$eid, "confirms_homology", fifelse(frac > 0.3, "conflicts_homology", "existing_only"))]
  add <- ex[class == "existing_only"]
  add_models <- add[, .(id = eid, seqnames, start, end, strand, exon_rank = 1L, evidence = "existing", alt_start = NA_integer_)]
  list(models = rbind(hm, add_models, fill = TRUE), report = ex)
}
