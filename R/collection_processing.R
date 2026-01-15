UMAP_by_gene_counts <- function(all_exp = list.experiments(validate = FALSE, pattern = "all_samples", libtypeExclusive = lib_type),
                                lib_type = "RFP",
                                organisms = unique(all_exp[!(organism %in% "all_species")]$name),
                                meta_rc = fread("~/livemount/Bio_data/NGS_pipeline/FINAL_LIST.csv"),
                                min_samples = 10) {

  if (!requireNamespace("matrixStats", quietly = TRUE) ||
      packageVersion("matrixStats") >= "0.63.0") {
    message("Installing matrixStats 0.62.0 for legacy scater compatibility...")
    remotes::install_version("matrixStats", version = "0.62.0", upgrade = "never")
  }
  library(matrixStats)
  library(SummarizedExperiment)
  library(scater)     # For logNormCounts and runUMAP
  library(uwot)       # Or use scater::runUMAP wrapper
  library(SingleCellExperiment)

  org <- organisms[1]
  for (org in unique(organisms)) {
    df_all <- read.experiment(org, validate = FALSE)
    message("- ", organism(df_all))
    if (nrow(df_all) < min_samples) {
      message("-- Too few samples to make UMAP, must be at least ", min_samples, "!")
      next
    }
    canonical_txs <- canonical_isoforms(df_all)
    count_table <- countTable(df_all, type = "summarized")
    count_table_canonical <- count_table[rownames(count_table) %in% canonical_txs]
    sce <- as(count_table_canonical, "SingleCellExperiment")
    # Normalize counts (or use your own transformation)
    sce <- logNormCounts(sce)
    sce <- runPCA(sce, ncomponents = min(50, nrow(df_all)-1))
    sce <- runUMAP(sce, dimred = "PCA", n_neighbors = min(30, nrow(df_all)-1), min_dist = 2)

    # Fix metadata
    srrs <- ORFik:::runIDs(df_all)
    rc_m <- meta_rc[match(srrs, Run),]
    stopifnot(nrow(rc_m) > 0)
    m <- rc_m

    m[AUTHOR %in% c("", "Makar"), AUTHOR := "Unknown"]
    m[CELL_LINE %in% c(""), CELL_LINE := "Unknown"]
    m[TISSUE %in% c(""), TISSUE := "Unknown"]
    m[INHIBITOR %in% c(""), INHIBITOR := "Unknown"]
    cell_lines <- m$CELL_LINE
    tissues <- m$TISSUE
    inhibitors <- m$INHIBITOR
    tissues_cell_lines <- paste0(tissues, " | ", cell_lines)
    author <- m$AUTHOR
    studies <- m$BioProject

    dt_umap <- as.data.table(reducedDim(sce, "UMAP"))
    colnames(dt_umap) <- c("UMAP 1", "UMAP 2")
    dt_umap[, `:=`(sample = m$Run,
                   cell_line = tolower(cell_lines),
                   tissue = tolower(tissues),
                   tissues_cell_lines = tolower(tissues_cell_lines),
                   inhibitors = tolower(inhibitors),
                   BioProject = studies,
                   author = m$AUTHOR,
                   studies = m$BioProject)]
    dir <- file.path(refFolder(df_all), "UMAP")
    dir.create(dir, FALSE, TRUE)
    type <- ifelse(lib_type == "RFP", "", paste0("_", lib_type))
    path <- file.path(dir, paste0("UMAP_by_gene_counts", type, ".fst"))
    fst::write_fst(dt_umap, path)
    message("Done , saved UMAP fst to location: ", path)
  }
  return(invisible(NULL))
}

save_coverage_tiles <- function(filepath, chromosomes, chromosome_lengths,
                                output_dir, tile_size, colnames, verbose = TRUE,
                                remake = FALSE, BPPARAM) {
  message("-- Loading file: ", basename(filepath))
  covrle_list_by_strand <- qs::qread(filepath, nthreads = 5)
  message("-- Processing all chromosomes: ", basename(filepath))
  index_dt <- rbindlist(lapply(chromosomes, function(chr) {
    fst::threads_fst(5, reset_after_fork = FALSE)

    total <- chromosome_lengths[chr]
    n_tiles <- ceiling(total / tile_size)
    indices <- seq_len(n_tiles)
    if (verbose) message("Processing chromosome: ", chr, " (Tile pages: ", n_tiles, ")")
    part_files <- file.path(output_dir, paste0("coverage_", chr, "_part", indices))
    if (attr(covrle_list_by_strand, "strand") == "+") {
      file_paths <- paste0(part_files, "_forward",".fst")
    } else {
      file_paths <- paste0(part_files, "_reverse",".fst")
    }

    irl <- List(lapply(indices, function(i) {
      idx_from <- ((i - 1) * tile_size) + 1
      idx_to <- min(i * tile_size, total)
      IRanges(idx_from, idx_to)
    }))
    chr_subset <- List(lapply(covrle_list_by_strand, function(covrle) covrle[[chr]]))
    if (verbose) cat("Tile Page:")
    lapply(indices, function(i) {
      if (verbose) cat(" ", i, ", ", sep = "")
      if (!file.exists(file_paths[i]) | remake) {
        # TODO: resolutions go here ->
        write_fst(setnames(setDT(lapply(chr_subset, function(x) as.integer(x[irl[[i]]]))), colnames),
                  file_paths[i])
      }
      return(invisible(NULL))
    })

    data.table(
      chr = chr,
      part = indices,
      start = start(unlist(irl)),
      end = end(unlist(irl)),
      file = file_paths
    )}))
}

# TODO: add multiple resolutions 1,3,9,27 (powers of 3 to 1 million  ish)
fst_megaformat_make <- function(all_exp = list.experiments(validate = FALSE, pattern = "all_samples", libtypeExclusive = lib_type),
                                lib_type = "RFP",
                                organisms = unique(all_exp[!(organism %in% "all_species")][samples > 10][order(samples, decreasing = TRUE)]$name),
                                min_samples = 10,
                                tile_size = 2e6,
                                force_remake_megacovrle = TRUE,
                                fst_threads = 5,
                                bp2 = 5, BPPARAM = BiocParallel::MulticoreParam(2, exportglobals = FALSE, log = FALSE)) {
  fst::threads_fst(fst_threads, reset_after_fork = FALSE)
  for (organism in organisms) {
    # ---- PARAMETERS ----
    # df <- read.experiment("all_samples-Homo_sapiens_1_9_2025", validate = FALSE)
    df <- read.experiment(organism, validate = FALSE)
    # df <- read.experiment("all_samples-Saccharomyces_cerevisiae", validate = FALSE)
    message("- ", organism(df))
    # df <- df[sample(nrow(df), 20, replace = FALSE),]
    output_dir <- file.path(resFolder(df), "collection_tables_indexed")
    index_path <- file.path(output_dir, "coverage_index.fst")
    if (file.exists(index_path)) next
    dir.create(output_dir, showWarnings = FALSE)

    file_paths <- file.path(resFolder(df), "megalist_covrle", paste0("all_covrles_", sub(" ", "_", organism(df)), c("_forward", "_reverse"), ".qs"))
    if (!all(file.exists(file_paths))) {
      message("- Creating megalist covrle")
      covrle_collection_path <- file.path(resFolder(df), "megalist_covrle", "all_covrles.qs")
      if (!file.exists(covrle_collection_path) | force_remake_megacovrle) {
        message("-- Loading all covrles..")
        cov_paths <- try(filepath(df, "cov", suffix_stem = c("", "_pshifted"), base_folders = libFolder(df, mode = "all")), silent = TRUE)
        if (is(cov_paths, "try-error")) {
          message("-- Organism does not have cov RLEs for all libraries, skipping!")
          next
        }
        all_covrles <- lapply(cov_paths, function(path) read_RDSQS(path, nthread = 5))
      } else {
        message("-- Loading premade file of all covrles..")
        all_covrles <- read_RDSQS(covrle_collection_path)
      }

      #
      covs_split <- list(lapply(all_covrles, function(x) f(x)),
                         lapply(all_covrles, function(x) r(x)))
      rm(all_covrles)
      gc()
      attr(covs_split[[1]], "strand") <- "+"
      attr(covs_split[[2]], "strand") <- "-"
      dir.create(dirname(file_paths[1]), showWarnings = FALSE, recursive = TRUE)
      qs2::qs_save(covs_split[[1]], file_paths[1], nthreads = 5)
      qs2::qs_save(covs_split[[2]], file_paths[2], nthreads = 5)
      rm(covs_split)
      gc()
    }

    stopifnot(length(file_paths) == 2 & all(file.exists(file_paths)))
    verbose <- TRUE
    colnames <- runIDs(df)
    seqinfo <- seqinfo(df)
    seqinfo <- seqinfo[seqnames(seqinfo)[order(seqlengths(seqinfo), decreasing = TRUE)]]
    seqinfo <- seqinfo[seqnames(seqinfo)[order(rep(seq(bp2), length.out = length(seqlengths(seqinfo))))]]
    chromosomes <- seqnames(seqinfo)
    chromosome_lengths <- seqlengths(seqinfo)

    message("-- Creating chromosome fst pages")
    message("--- Total chromosomes: ", length(chromosomes))
    message("--- Total pages: ", sum(ceiling(seqlengths(seqinfo) / tile_size)))
    index_dt <- lapply(file_paths, function(filepath, chromosomes,
                                            chromosome_lengths,
                                            output_dir, tile_size, colnames,
                                            save_coverage_tiles, bp, verbose) {
      save_coverage_tiles(filepath, chromosomes, chromosome_lengths,
                          output_dir, tile_size, colnames, verbose,
                          BPPARAM = BiocParallel::MulticoreParam(bp, exportglobals = FALSE, log = FALSE))
    },
    chromosomes = chromosomes,
    chromosome_lengths = chromosome_lengths,
    output_dir = output_dir, tile_size = tile_size, verbose = verbose,
    colnames = colnames,
    save_coverage_tiles = save_coverage_tiles,
    bp = bp2)

    colnames(index_dt[[1]])[colnames(index_dt[[1]]) == "file"] <- "file_forward"
    colnames(index_dt[[2]])[colnames(index_dt[[2]]) == "file"] <- "file_reverse"
    index_dt <- cbind(index_dt[[1]], file_reverse = index_dt[[2]]$file_reverse)

    # Save index file
    write_fst(index_dt, index_path)
  }
  message("Done")
}

meta_meta_motif_motif <- function(all_exp = list.experiments(validate = FALSE, pattern = "all_samples", libtypeExclusive = lib_type),
                                  lib_type = "RFP",
                                  organisms = unique(all_exp[!(organism %in% "all_species")][samples > 10][order(samples, decreasing = TRUE)]$name),
                                  force_remake_megacovrle = TRUE,
                                  letters_vec = c("P", "R")) {
  remake <- force_remake_megacovrle
  for (org in organisms) {
    df <- read.experiment(org, validate = FALSE)
    filepaths <- try(filepath(df, "cov", base_folders = libFolder(df, mode = "all"), suffix_stem = c("_pshifted", "")))
    if (is(filepaths, "try-error")) {
      warning("Missing covrle files for organism, going to next")
      next
    }
    lib_names <- bamVarName(df)

    org_short <- gsub(" ", "_", tolower(organism(df)))
    message("- ", org_short)
    org_dir <- file.path(dirname(QCfolder(df)), "meta_collection_tables")
    TIS_file <- file.path(org_dir, "meta_collection_table_TIS.fst")
    TTS_file <- file.path(org_dir, "meta_collection_table_TTS.fst")
    cds_anchor_files <- c(TIS = TIS_file, TTS = TTS_file)

    all_AA_letters <- Biostrings::AA_STANDARD
    permute <- function(x) {
      if (length(x) == 1) return(x)
      res <- c()
      for (i in seq_along(x)) {
        sub_perms <- permute(x[-i])
        res <- c(res, paste0(x[i], sub_perms))
      }
      return(res)
    }

    # Get permutations as strings

    # Generate all 2-character combinations with repetition
    pairs <- as.vector(outer(letters_vec, letters_vec, paste0))
    aa_pair_files <- file.path(org_dir, paste0("meta_collection_table_",  pairs,".fst"))
    names(aa_pair_files) <- pairs
    files <- c(cds_anchor_files, aa_pair_files)
    if (all(file.exists(files)) & !remake) next
    if (!dir.exists(org_dir)) dir.create(org_dir, recursive = TRUE)

    txdb <- loadTxdb(df)
    subset <- integer()
    if (org_short == "homo_sapiens") {
      mane <- canonical_isoforms(df)
      subset <- mane
    } else try(subset <- filterTranscripts(txdb, NULL, 30, NULL, longestPerGene = TRUE), silent = TRUE)
    all_cds <- cds <- loadRegion(txdb, "cds")
    mrna_original <- loadRegion(txdb, "mrna")
    cds <- all_cds[names(all_cds) %in% subset]
    mrna <- mrna_original[names(mrna_original) %in% subset]
    stopifnot(length(cds) == length(mrna)); length(cds)
    extension <- 300
    mrna <- extendTrailersUntil(extendLeadersUntil(mrna, cds, extension), cds, extension)

    anchor_site <- "TIS"
    for (anchor_site in names(files)) {
      message(anchor_site)
      if (file.exists(files[anchor_site]) & !remake) next
      anchor <- if (anchor_site == "TIS") {
        startSites(cds, TRUE, TRUE, TRUE)
      } else if (anchor_site == "TTS") {
        stopSites(cds, TRUE, TRUE, TRUE)
      } else if (anchor_site %in% pairs) {
        seqs <- txSeqsFromFa(cds, df, is.sorted = TRUE, keep.names = TRUE)
        codon_seqs <- translate(seqs, if.fuzzy.codon = "solve")
        pattern_match_irl <- IRangesList(Biostrings::vmatchPattern(anchor_site, codon_seqs))
        names(pattern_match_irl) <- seq_along(pattern_match_irl)
        pattern_match_irl <- pattern_match_irl[lengths(pattern_match_irl) > 0]
        width(pattern_match_irl) <- 1
        end(pattern_match_irl) <- IntegerList((end(pattern_match_irl)*3)-2L)
        start(pattern_match_irl) <- IntegerList((start(pattern_match_irl)*3)-2L)
        anchor <- unlistGrl(pmapFromTranscriptF(pattern_match_irl, cds, removeEmpty = T))

      } else stop("Invite anchor site")
      # cds_stops <- split(cds_stops, names(cds_stops))
      length_original <- length(anchor)
      windows_all <- windowPerGroup(anchor, tx = mrna, upstream = extension, downstream = extension - 1)
      summary(widthPerGroup(windows_all))
      windows_all <- windows_all[widthPerGroup(windows_all) == extension*2]
      length_filtered <- length(windows_all)
      length_diff <- length_original - length_filtered
      message("Lost ", length_diff, " anchor sites (", round(100 - (100*(length_filtered/length_original)), 2), "%)")
      max_ranges <- 30000
      if (length_filtered > max_ranges) {
        message("Original anchor sites: ", length_original)
        windows_all <- windows_all[sample(length_filtered, max_ranges, replace = FALSE)]
        subsampling_length <- length(windows_all)
        length_diff <- length_original - subsampling_length
        message("Lost ", length_diff, " anchor sites (", round(100 - (100*(subsampling_length/length_original)), 2), "%)")
      }


      print(Sys.time())
      dt <- rbindlist(
        bplapply(seq_along(filepaths), function(i, windows_all, filepaths, lib_names) {
          file <- filepaths[i]
          lib_name <- lib_names[i]
          d <-
            data.table(library = factor(lib_name),
                       count = coverageScorings(coveragePerTiling(windows_all, fimport(file), as.data.table = TRUE, is.sorted = TRUE), "transcriptNormalized")$score)
        }, windows_all = windows_all, filepaths = filepaths, lib_names = lib_names))
      print(Sys.time())
      fst::write_fst(dt, files[anchor_site])
    }
  }
}

