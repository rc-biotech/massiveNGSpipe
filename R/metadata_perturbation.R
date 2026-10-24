#' Detect gene perturbations (KO, KD, OE, mutant, degron) and their controls in sample titles
#'
#' Rules on the title (keywords case-insensitive except the morpholino "MO"), applied in this order; the first rule
#' with a valid gene wins:
#' \itemize{
#'  \item control: siCtrl/siNT/siNC/siScr, shCtrl/shScr/shNT/shGFP/shLuc, sgCtrl/sgNT/sgGal4/sgAAVS1, scrambled,
#'   non-targeting, NTC, mock-transfected -> CONDITION "control"; empty vector / EV -> "EV"
#'  \item mutant: <GENE> <point mutation> (RPL12 S38A, asc1-M1X) -> "mutant", GENE "<GENE>_<mutation>"
#'  \item degron: <GENE>-AID / -mAID / -dTAG / -FKBP12F36V / -degron -> "degron"
#'  \item DKO: <G1>/<G2> DKO, double knockout -> "DKO", GENE "<G1>_<G2>"
#'  \item KD: si<GENE>, sh<GENE>, <GENE> siRNA/shRNA/RNAi/knockdown/KD, <GENE> MO (morpholino), siRNA/KD <GENE>,
#'   sg<GENE> with CRISPRi/dCas9-KRAB
#'  \item KO: <GENE>-/-, <GENE>\u0394, \u0394<GENE>, delta<GENE>, <GENE> KO/knockout/null/deletion, KO <GENE>, sg<GENE>
#'  \item OE: <GENE> OE/overexpression/OX, overexpressing <GENE>
#' }
#' A gene token is accepted when it is a gene symbol of the sample's organism (case-insensitive, an optional clone
#' suffix such as "-1" removed). Without a symbol set only tightly attached forms (si/sh/sg prefix, \u0394, -/-,
#' attached KO/KD) are accepted.
#' @param titles character vector of sample titles
#' @param organisms character vector (same length) of scientific names, used to pick symbol sets
#' @param symbols named list: scientific name -> character vector of gene symbols (lower case)
#' @return data.table(CONDITION, GENE, rule); empty strings where nothing was detected
#' @export
metadata_perturbation_labels <- function(titles, organisms = rep("", length(titles)), symbols = list()) {
  t <- ifelse(is.na(titles), "", as.character(titles))
  # SRA "Platform sequencing: GSM123: <title>; Organism; strategy" decoration
  t <- sub("^.*?sequencing: *GSM[0-9]+: *", "", t, perl = TRUE)
  for (o in unique(organisms[nzchar(organisms)])) t[organisms == o] <- gsub(o, " ", t[organisms == o], fixed = TRUE)
  t <- gsub("(?i)\\b(RNA-?Seq|Ribo-?Seq|OTHER|ribosome profiling)\\s*$", " ", t, perl = TRUE)
  n <- length(t)
  out <- data.table::data.table(CONDITION = rep("", n), GENE = rep("", n), rule = rep("", n))
  not_gene <- c("control", "ctrl", "con", "nc", "nt", "neg", "scr", "scramble", "scrambled", "luc", "gfp", "lacz", "empty",
                "mock", "rep", "replicate", "sample", "cell", "cells", "wt", "total", "input", "rpf", "chx", "and", "with",
                "the", "of", "no", "day", "min", "cas9", "vector", "ev", "double", "cre", "flox", "ip", "ko", "kd", "oe",
                "rna", "mrna", "seq", "ribo", "set", "max", "per", "ice", "cat", "can", "act", "gal4", "aavs1", "renilla",
                "mcherry", "yfp", "rfp", "egfp", "puro", "dox")
  norm <- function(g) sub("[-_.][0-9]{1,2}$", "", g)
  valid <- function(g, i, strict) {
    ok <- !is.na(g) & nzchar(g) & grepl("^[A-Za-z][A-Za-z0-9.-]{0,14}$", g) & !tolower(g) %in% not_gene
    has_set <- organisms[i] %in% names(symbols)
    in_set <- vapply(seq_along(i), function(k) has_set[k] && ok[k] &&
                       (tolower(g[k]) %in% symbols[[organisms[i[k]]]] || tolower(norm(g[k])) %in% symbols[[organisms[i[k]]]]),
                     logical(1))
    ifelse(has_set, in_set, ok & strict & nchar(g) >= 3)
  }
  # all matches of capture group 1 per title
  cands <- function(re) lapply(regmatches(t, gregexpr(re, t, perl = TRUE)), function(m)
    if (length(m)) sub(re, "\\1", m, perl = TRUE) else character())
  apply_rule <- function(re, cond, rule, strict, mk = identity, only = rep(TRUE, n)) {
    todo <- which(out$CONDITION == "" & only)
    if (!length(todo)) return(invisible())
    cl <- cands(re)[todo]
    for (k in seq_along(todo)) {
      g <- cl[[k]]; if (!length(g)) next
      v <- valid(g, rep(todo[k], length(g)), strict)
      if (any(v)) { gg <- g[which(v)[1]]
        data.table::set(out, todo[k], c("CONDITION", "GENE", "rule"), list(cond, mk(norm(gg)), rule)) }
    }
  }
  # 1. controls
  ctrl <- grepl("(?i)(?<![a-z])(si|sh|sg)[-_ ]?(ctrl|control|con|nc|nt|neg|scr|scramble[d]?|luc|gfp|lacz|non[-_ ]?target|gal4|aavs1|renilla|empty)", t, perl = TRUE) |
    grepl("(?i)scrambled?|non[-_ ]?targeting|(?<![a-z])ntc(?![a-z])|mock[-_ ]?transfect", t, perl = TRUE)
  # negated perturbation ("untreated cells (no HTS1 knock-down)", "without siRNA") is the control series
  ctrl <- ctrl | grepl("(?i)(?<![a-z])(no|without|w/o|non)[-_ ]+(?:[A-Za-z0-9.]+[-_ ])?(knock-?down|knock-?out|KD|KO|siRNA|shRNA|overexpression|OE)(?![a-z])", t, perl = TRUE)
  data.table::set(out, which(ctrl), c("CONDITION", "rule"), list("control", "control"))
  ev <- out$CONDITION == "" & grepl("(?i)empty[-_ ]?vector|(?<![a-z0-9])ev(?![a-z0-9])", t, perl = TRUE)
  data.table::set(out, which(ev), c("CONDITION", "rule"), list("EV", "empty_vector"))
  # 2. point mutants (gene + mutation captured together, validated on the gene part)
  pm <- "(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9]{1,9}[-_ ][ACDEFGHIKLMNPQRSTVWY][0-9]{1,4}(?:[ACDEFGHIKLMNPQRSTVWY]|X|del|\\*))(?![A-Za-z0-9])"
  todo <- which(out$CONDITION == ""); pm_all <- cands(pm)
  for (k in todo) { m <- pm_all[[k]]; if (!length(m)) next
    gene <- sub("[-_ ].*$", "", m); mut <- sub("^[^-_ ]*[-_ ]", "", m)
    # 2A self-cleaving peptides are tags, not mutations; without a symbol set require a >= 2-digit position
    v <- valid(gene, rep(k, length(m)), TRUE) & !grepl("^[PTEF]2A$", mut) &
      (organisms[k] %in% names(symbols) | grepl("^[A-Z][0-9]{2,}", mut))
    if (any(v)) { j <- which(v)[1]
      data.table::set(out, k, c("CONDITION", "GENE", "rule"), list("mutant", sub("[-_ ]", "_", m[j]), "point_mutant")) } }
  # 3. degron
  apply_rule("(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9]{1,14})[-_ ]?(?i:mAID|AID|dTAG|FKBP12F36V|degron)(?![A-Za-z])", "degron", "degron", TRUE)
  # 4. double knockout
  todo <- which(out$CONDITION == "")
  dko <- "(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9]{1,9})[/_+-]([A-Za-z][A-Za-z0-9]{1,9})[-_ ]?(?i:DKO|double[-_ ]?(?:KO|knock-?out|mutant))(?![A-Za-z])"
  for (k in todo) { m <- regmatches(t[k], regexec(dko, t[k], perl = TRUE))[[1]]; if (length(m) < 3) next
    if (all(valid(m[2:3], c(k, k), TRUE)))
      data.table::set(out, k, c("CONDITION", "GENE", "rule"), list("DKO", paste(m[2], m[3], sep = "_"), "double_ko")) }
  # 5. knockdown
  crispri <- grepl("(?i)crispri|dcas9[-_ ]?krab", t, perl = TRUE)
  # si/sh prefix: without a symbol set the target must look like a symbol (capital first letter or a digit), not a word
  si_words <- grepl("(?i)(?<![a-z])(sh(ock|ort|ift|am|oot|eet|ell|ape|ear|aker)|si(ngle|mple|ze|lent|gnal|gma|te|des|lico|milar|ster|x))", t, perl = TRUE)
  apply_rule("(?<![A-Za-z0-9])(?:si|sh)[-_]?([A-Z][A-Za-z0-9.-]{1,14}|[a-z][a-z]*[0-9][A-Za-z0-9.-]*)(?![A-Za-z0-9])", "KD", "si_sh_prefix", TRUE, only = !si_words)
  apply_rule("(?<![A-Za-z0-9])(?:si|sh)[-_]?([a-z][A-Za-z0-9.-]{1,14})(?![A-Za-z0-9])", "KD", "si_sh_prefix", FALSE, only = !si_words)
  apply_rule("(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9.-]{1,14})[-_ ]?(?:(?i:siRNA|shRNA|RNAi|knock-?down|KD)(?![A-Za-z])|MO(?![A-Za-z]))", "KD", "gene_kd_suffix", FALSE)
  apply_rule("(?:(?i:siRNA|shRNA|RNAi|dsRNA|knock-?down|KD)|MO)[-_ ](?i:against[-_ ]|of[-_ ]|targeting[-_ ])?([A-Za-z][A-Za-z0-9.-]{1,14})(?![A-Za-z0-9])", "KD", "kd_prefix", FALSE)
  apply_rule("(?<![A-Za-z0-9])sg[-_]?([A-Za-z][A-Za-z0-9.-]{1,14})(?![A-Za-z0-9])", "KD", "crispri_sg", TRUE, only = crispri)
  # 6. knockout
  ko_null <- "(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9._]{1,14})(?:-/-|\u0394|\u2206)"
  ko_delta <- "(?:\u0394|\u2206|(?<![A-Za-z])(?i:delta)[-_ ]?)([A-Za-z][A-Za-z0-9._]{1,14})(?![A-Za-z0-9])"
  # several deleted genes in one title (cth1-delta cth2-delta) -> DKO
  todo <- which(out$CONDITION == ""); kn <- cands(ko_null); kd <- cands(ko_delta)
  for (k in todo) { g <- unique(norm(c(kn[[k]], kd[[k]]))); if (length(g) < 2) next
    v <- valid(g, rep(k, length(g)), TRUE); if (sum(v) >= 2)
      data.table::set(out, k, c("CONDITION", "GENE", "rule"), list("DKO", paste(g[v][1:2], collapse = "_"), "multi_deletion")) }
  apply_rule(ko_null, "KO", "gene_null", TRUE)
  apply_rule(ko_delta, "KO", "delta_prefix", TRUE)
  apply_rule("(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9.]{1,14})[-_]?(?i:KO|knock-?out|null)(?![A-Za-z])", "KO", "gene_ko_attached", TRUE)
  apply_rule("(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9.]{1,14}) (?i:KO|knock-?out|null|deletion)(?![A-Za-z])", "KO", "gene_ko_suffix", FALSE)
  apply_rule("(?<![A-Za-z])(?i:KO|knock-?out)[-_ ](?i:of[-_ ])?([A-Za-z][A-Za-z0-9.-]{1,14})(?![A-Za-z0-9])", "KO", "ko_prefix", FALSE)
  apply_rule("(?<![A-Za-z0-9])sg[-_]?([A-Za-z][A-Za-z0-9.-]{1,14})(?![A-Za-z0-9])", "KO", "crispr_sg", TRUE, only = !crispri)
  # 7. overexpression
  apply_rule("(?<![A-Za-z0-9])([A-Za-z][A-Za-z0-9.-]{1,14})[-_ ]?(?i:OE|overexpress(?:ion|ed|ing)?|OX)(?![A-Za-z])", "OE", "gene_oe_suffix", FALSE)
  apply_rule("(?i:overexpress(?:ion|ed|ing)?|OE)[-_ ](?i:of[-_ ])?([A-Za-z][A-Za-z0-9.-]{1,14})(?![A-Za-z0-9])", "OE", "oe_prefix", FALSE)
  out[]
}

#' Gene symbol sets (lower case) per scientific name from the reference folders (gene_symbol_tx_table.fst)
#' @param organisms scientific names
#' @param ref_dir references folder
#' @return named list; organisms without a symbol table are left out
#' @export
metadata_symbol_sets <- function(organisms, ref_dir = "~/livemount/Bio_data/references") {
  organisms <- unique(organisms[!is.na(organisms) & nzchar(organisms)])
  res <- list()
  for (o in organisms) {
    f <- file.path(ref_dir, gsub(" ", "_", tolower(o)), "gene_symbol_tx_table.fst")
    if (!file.exists(f)) next
    x <- fst::read_fst(f, columns = intersect(c("external_gene_name", "ensembl_gene_id"), fst::metadata_fst(f)$columnNames))
    s <- unique(tolower(c(x$external_gene_name, x$ensembl_gene_id)))
    res[[o]] <- s[!is.na(s) & nzchar(s)]
  }
  res
}

#' Fill empty CONDITION / GENE from detected gene perturbations
#'
#' Only empty fields are filled; existing labels are never changed.
#' @param dt data.table with sample_title, ScientificName, CONDITION, GENE
#' @param symbols see \code{\link{metadata_symbol_sets}}
#' @return list(dt = updated copy, changes = data.table(row, column, old, new, rule))
#' @export
metadata_fill_perturbations <- function(dt, symbols = metadata_symbol_sets(dt$ScientificName)) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  for (cl in c("CONDITION", "GENE")) data.table::set(dt, j = cl, value = { v <- as.character(dt[[cl]]); v[is.na(v)] <- ""; v })
  p <- metadata_perturbation_labels(dt$sample_title, as.character(dt$ScientificName), symbols)
  vocab <- c("WT", "KO", "KD", "OE", "DKO", "MUTANT", "KI", "EV", "CONTROL", "DEGRON")
  cond_u <- toupper(dt$CONDITION)
  has_gene_in_cond <- mapply(function(cnd, g) nzchar(g) && grepl(tolower(g), tolower(cnd), fixed = TRUE), dt$CONDITION, p$GENE)
  fill_c <- which(dt$CONDITION == "" & p$CONDITION != "")
  gene_ok <- p$GENE != "" & (dt$CONDITION == "" | cond_u == toupper(p$CONDITION) | (!cond_u %in% vocab & has_gene_in_cond))
  fill_g <- which(dt$GENE == "" & gene_ok)
  # GENE left by the old KO rule (title minus "KO", e.g. "MV4-11_RNAseq_3_RPL22"): contains the detected gene plus title text
  old_rule <- dt$GENE != "" & gene_ok & toupper(dt$GENE) != toupper(p$GENE) &
    mapply(function(old, g) nzchar(g) && grepl(tolower(g), tolower(old), fixed = TRUE) && grepl("[^A-Za-z0-9]", old) &&
             nchar(old) > nchar(g) + 1, dt$GENE, p$GENE)
  # keep compound values naming a second gene (ski2_dom34KO, HSP70_PELO): they describe double perturbations
  second_gene <- mapply(function(old, g, o) {
    tok <- tolower(unlist(strsplit(gsub("(?i)ko|kd|oe", " ", old, perl = TRUE), "[^A-Za-z0-9.]+")))
    tok <- tok[nzchar(tok) & tok != tolower(g)]
    any(tok %in% symbols[[o]]) || (!o %in% names(symbols) && length(tok) > 0)
  }, dt$GENE, p$GENE, as.character(dt$ScientificName))
  fix_g <- which(old_rule & !second_gene)
  ch <- data.table::rbindlist(list(
    data.table::data.table(row = fill_c, column = rep("CONDITION", length(fill_c)), old = "", new = p$CONDITION[fill_c], rule = p$rule[fill_c]),
    data.table::data.table(row = fill_g, column = rep("GENE", length(fill_g)), old = "", new = p$GENE[fill_g], rule = p$rule[fill_g]),
    data.table::data.table(row = fix_g, column = rep("GENE", length(fix_g)), old = dt$GENE[fix_g], new = p$GENE[fix_g],
                           rule = paste0(p$rule[fix_g], "+old_KO_rule_gene"))))
  data.table::set(dt, fill_c, "CONDITION", p$CONDITION[fill_c]); data.table::set(dt, c(fill_g, fix_g), "GENE", p$GENE[c(fill_g, fix_g)])
  list(dt = dt[], changes = ch)
}
