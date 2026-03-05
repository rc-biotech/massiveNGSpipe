#' For docker usage copy
#'
#' If you have download docker and release docker.
#' This is useful to copy files from download to release
#' @param config a config object
#' @param pipelines pipeline objects, default:
#' \code{pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE)}
#' @param new_exp_dir character, default "~/livemount/Bio_data/ORFik_experiments/"
#' @param docker_conversion the sed conversion string: 's/livemount\///g'
#' @param individual_studies logical, default TRUE. Should raw ORFik experiments
#' be exported per study.
#' @param merged_by_organism logical, default TRUE.
#' Create a new study, merge the superset of all per organism into single file.
#' This requires that you have run: pipeline_merge_org
#' @param all_by_organism logical, default TRUE.
#' Create a new study, the superset of all studies of all per organism.
#' This requires that you have run: pipeline_merge_org
#' @param all_modalities_by_org logical, default TRUE.
#' Should the all modalities tracks be exported per organism.
#' @return logical, TRUE if sucessful for all.
#' @export
docker_copy_done_experiments <- function(config,
                                         pipelines = pipeline_init_all(config, gene_symbols = FALSE, only_complete_genomes = TRUE),
                                         new_exp_dir = "~/livemount/Bio_data/ORFik_experiments/",
                                         docker_conversion = "'s/livemount\\///g'",
                                         individual_studies = TRUE,
                                         merged_by_organism = TRUE,
                                         all_by_organism = TRUE,
                                         all_modalities_by_org = TRUE) {
  if (length(pipelines) == 0) stop("Empty pipelines object given!")
  message("- Individual")
  done_stats <- progress_report(pipelines, config, return_progress_vector = TRUE)
  max_step <- length(unlist(config$flag_steps))
  done_all_steps <- done_stats == max_step
  res <- c()
  csv_names <- paste0(pipelines_names(pipelines), ".csv")
  csv_names <- csv_names[done_all_steps]
  old_exp_dir <- config$config["exp"]

  if (length(csv_names) == 0) {
    warning("No experiments are done, returning directly!")
    return(FALSE)
  }
  if (individual_studies) {
    res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
  }

  if (merged_by_organism | all_by_organism | all_modalities_by_org) {
    names(pipelines) <- NULL
    exp <- lapply(pipelines, function(x) lapply(x$organisms, function(o) o$conf["exp"]))
    exp <- unlist(exp, recursive = FALSE)
    # Merge all per organism
    done_exp_list <- exp[done_all_steps]
    done_organisms <- unique(names(done_exp_list))

    if (merged_by_organism) {
      message("- Merged by organism")
      csv_names <- paste0(organism_merged_exp_name(done_organisms), ".csv")
      existing_exps <- list.files(config$config["exp"])
      merged_is_made <- csv_names %in% existing_exps
      csv_names <- csv_names[merged_is_made]
      res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
    }

    if (all_by_organism) {
      message("- All by organism")
      csv_names <- paste0(organism_collection_exp_name(done_organisms), ".csv")
      existing_exps <- list.files(config$config["exp"])
      merged_is_made <- csv_names %in% existing_exps
      csv_names <- csv_names[merged_is_made]

      res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
    }

    if (all_modalities_by_org) {
      message("- All modalities by organism")
      csv_names <- paste0(organism_merged_modalities_exp_name(done_organisms), ".csv")
      existing_exps <- list.files(config$config["exp"])
      merged_is_made <- csv_names %in% existing_exps
      csv_names <- csv_names[merged_is_made]

      res <- c(res, copy_experiments_to(csv_names, old_exp_dir, new_exp_dir, docker_conversion))
    }
  }
  return(all(res == 0))
}

#' Copy done experiment csvs to new directory
#' @examples
#' csvs <- list.files(ORFik::config()["exp"], pattern = "modalities\\.csv$")
#' copy_experiments_to(csvs)
copy_experiments_to <- function(csv_names, old_exp_dir = ORFik::config()["exp"],
                                new_exp_dir = sub("_local$|_local/$", "", old_exp_dir),
                                docker_conversion = "'s/livemount\\///g'") {
  stopifnot(old_exp_dir != new_exp_dir)
  stopifnot(dir.exists(old_exp_dir)); stopifnot(dir.exists(new_exp_dir))
  done_exp_old_path <- file.path(old_exp_dir, csv_names)
  done_exp_new_path <- file.path(new_exp_dir, csv_names)
  res <- c()
  for (i in seq_along(csv_names)) {
    file.copy(done_exp_old_path[i], done_exp_new_path[i], overwrite = TRUE)
    cmd <- system(paste0("sed -i ", docker_conversion, " ",done_exp_new_path[i]))
    res <- c(res, cmd)
  }
  return(res == 0)
}

is_config <- function(config) return(is(config, "list") && !is.null(config$preset))

#' Get name of function
#' @importFrom rlang ensyms
#' @importFrom purrr map
#' @param ... functions
#' @return character, name of function
#' @export
name_of_function <- function(...) unlist(purrr:::map(rlang::ensyms(...) , as.character), use.names = FALSE)

#' @inherit rstudioapi::filesPaneNavigate
#' @importFrom rstudioapi filesPaneNavigate
#' @return invisible(NULL)
#' @export
navigate <- rstudioapi::filesPaneNavigate

#' Navigate to a experiment Directory in the Files Pane
#'
#' Navigate to a directory in the Files pane.
#' The contents of that directory will be listed and shown in the Files pane.
#' @param exp ORFik experiment or character name of experiment,
#' a ORFik experiment name of folder to go to, e.g. "GSE91068-saccharomyces_cerevisiae".
#' @param rel_dir character, relative directory, "aligned", "trim", c("", "processed") or "raw"
#' @return invisible(NULL)
#' @export
navigateExp <- function(exp, rel_dir = "aligned") {
  navigate(dirExp(exp, rel_dir))
}

#' @inherit navigateExp
#' @title Get path to dir type of exp
#' @return a character path
#' @export
dirExp <- function(exp, rel_dir = "aligned") {
  if (!is.character(exp)) exp <- name(exp)
  sub_dir <- if(grepl("_SSU$|_RNA-seq$|_disome$", exp)) {
    sub <- sub(".*_", "", exp)
    exp <- sub(paste0("_", sub), "", exp)
    sub
  } else ""

  main_dir <- file.path(ORFik::config()["bam"], sub_dir)
  path <-
    if (rel_dir == "aligned") {
      file.path(main_dir, exp, "aligned")
    } else if (rel_dir == "trim") {
      file.path(main_dir, exp, "trim")
    } else if (rel_dir %in% c("", "processed")) {
      file.path(main_dir, exp)
    } else if (rel_dir %in% c("collapsed")) {
      file.path(main_dir, exp, "trim", "SINGLE")
    } else {
      file.path(ORFik::config()["fastq"], sub_dir, exp)
    }
  return(path)
}

#' Navigate to a reference Directory in the Files Pane
#'
#' Navigate to a directory in the Files pane.
#' The contents of that directory will be listed and shown in the Files pane.
#' @inheritParams navigateExp
#' @param rel_dir character, relative directory, default: ""
#' @return invisible(NULL)
#' @export
navigateRef <- function(exp, rel_dir = "") {
  if (!is.character(exp)) exp <- organism(exp)
  organism <- sub(" |-", "_" , tolower(exp))
  path <- file.path(ORFik::config()["ref"], organism)

  navigate(path)
}

#' Get all experiment names to search with auto completion
#' @inheritParams ORFik::list.experiments
#' @return a named list of characters,
#'  so you can search experiment names with auto completion
#' @export
list.experiments.list <- function(pattern = "*", libtypeExclusive = NULL,
                                  dir = ORFik::config()["exp"]) {
  all_exp <- list.experiments(dir = dir, validate = FALSE, pattern = pattern,
                              libtypeExclusive = libtypeExclusive)
  return(split(all_exp$name, all_exp$name))
}

#' Get massiveNGSpipe related file counts and total drive usage
#' @param config config = ORFik::config()
#' @param include_references = TRUE
#' @return a data.table
#' @export
#' @examples
#' config <- ORFik::config()
#' config_RNA <- file.path(config, "RNA-seq")
#' config_SSU <- file.path(config, "SSU")
#' config_disome <- file.path(config, "disome")
#' names(config_RNA) <- names(config_SSU) <- names(config_disome) <- names(config)
#' info_dt <- rbindlist(
#' list(RFP =    file_statistics_all_types(config),
#'      RNA =    file_statistics_all_types(config_RNA, FALSE),
#'      SSU =    file_statistics_all_types(config_SSU, FALSE),
#'      disome = file_statistics_all_types(config_disome, FALSE)), idcol = TRUE)
#'  info_dt[, sum(total_size_GB), by = .id]
file_statistics_all_types <- function(config = ORFik::config(), include_references = TRUE) {
  stopifnot(!is.null(names(config)))
  rbindlist(list(file_statistics_processed(config),
                 file_statistics_raw(config),
                 if (include_references) {file_statistics_references(config)}))
}

file_statistics_processed <- function(config = ORFik::config(), processed_dir = config["bam"]) {
  l_main <- list.files(processed_dir, full.names = T)
  types <- c("aligned", "trim")
  res <- lapply(types, function(type) {
    l <- file.path(l_main, type)
    l <- l[dir.exists(l)]

    if (type == "aligned") {
      dir_types <- c("", "ofst", "pshifted", "cov_RLE", "cov_RLE_List","bigwig")
      formats <- c("bam", "ofst", "ofst", "covqs", "covqs","bigWig")
    } else {
      dir_types <- c("", "SINGLE", "PAIRED")
      formats <- c("fastq", "fasta", "fasta")
    }
    file_statistics_internal(dir_types, formats, type, l)
  })
  return(rbindlist(res))
}

file_statistics_references <- function(config = ORFik::config(),
                                       reference_dir = config["ref"]) {
  l <- list.files(reference_dir, full.names = T)
  type <- "reference"
  dir_types <- c("", "", "", "STAR_index/genomeDir/")
  formats <- c("(fasta|fa|fna)", "(gtf|gff|gff2|gff3)", "(db)", "")

  file_statistics_internal(dir_types, formats, type, l)
}

file_statistics_raw <- function(config = ORFik::config(),
                                raw_dir = config["fastq"]) {
  l <- list.files(raw_dir, full.names = T)
  type <- "raw"
  dir_types <- c("", "multiplexed")
  formats <- c("(fasta|fastq)(\\.gz)?", "(fasta|fastq)(\\.gz)?")

  file_statistics_internal(dir_types, formats, type, l)
}

file_statistics_internal <- function(dir_types, formats, type, l) {
  res <- mapply(function(dir_type, format) {
    d <- file.path(l, dir_type)
    d <- d[dir.exists(d)]
    if (dir_type == "") dir_type <- format
    format <- paste0(ifelse(format != "", "\\.", ""), format, "$")
    all_files <- list.files(d, format, full.names = TRUE)
    info <- as.data.table(file.info(all_files))
    data.table(type, dir_type, n_dirs = length(d), n_files = length(all_files), total_size_GB = round(sum(info$size) / 1e9))
  }, dir_types, formats, SIMPLIFY = FALSE)
  rbindlist(res)
}

conda_init_env_string <- function(env = "blast", conda = "~/miniconda3/bin/conda") {
  conda_works <- system(paste(conda, "info"), intern = TRUE)
  conda_works <- is.null(attr(conda_works, "status"))
  if (!conda_works) stop("Conda is was not found in path specified!")
  conda_init_env <- paste(conda, 'run -n', env)
  return(conda_init_env)
}
