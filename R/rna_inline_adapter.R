# Cache the adapter decision for the actual Run/mate being aligned. A study
# can mix long RNA-seq libraries with fragmented small-RNA-kit inputs.
rna_inline_adapter <- function(file, run, output_dir) {
  stopifnot(length(file) == 1L, file.exists(file),
            length(run) == 1L, grepl("^[A-Za-z0-9_-]+$", run))
  # A manual row always wins (as in detect_adapter()): the FastQC call can pick a downstream adapter segment and leave
  # the upstream part on the reads (GSE119903 medaka 2026-11-27: cut at TCGTATGCCGTC, AGATCGGAAGAG left, 0% aligned).
  manual <- rna_manual_adapter(run, file.path(output_dir, "trim"))
  if (!is.na(manual)) return(manual)
  signature <- list(path = normalizePath(file), size = file.info(file)$size,
                    mtime = as.numeric(file.info(file)$mtime), version = 1L)
  cache_dir <- file.path(output_dir, "aligned", ".inline_trim")
  cache <- file.path(cache_dir, paste0(run, "_adapter.rds"))
  old <- if (file.exists(cache)) readRDS(cache) else NULL
  if (!is.null(old) && identical(old$signature, signature)) return(old$adapter)
  adapter <- short_insert_adapter(file)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = cache_dir)
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(list(signature = signature, adapter = adapter), tmp)
  if (!file.rename(tmp, cache)) stop("Cannot save per-Run inline adapter decision")
  adapter
}

#' Adapter of this Run in <trim>/adapters_manual.csv (written by apply_trim_fix_to_sample()), else NA
#' @noRd
rna_manual_adapter <- function(run, trimmed_dir) {
  f <- file.path(trimmed_dir, "adapters_manual.csv")
  if (!file.exists(f)) return(NA_character_)
  m <- data.table::fread(f, header = TRUE, colClasses = "character")
  if (!all(c("Run", "adapter") %in% names(m))) return(NA_character_)
  a <- m[["adapter"]][m[["Run"]] == run]
  if (length(a) != 1L || is.na(a) || !nzchar(a)) return(NA_character_)
  stats::setNames(a, "manual")
}
