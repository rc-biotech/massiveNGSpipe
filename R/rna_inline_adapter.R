# Cache the adapter decision for the actual Run/mate being aligned. A study
# can mix long RNA-seq libraries with fragmented small-RNA-kit inputs.
rna_inline_adapter <- function(file, run, output_dir) {
  stopifnot(length(file) == 1L, file.exists(file),
            length(run) == 1L, grepl("^[A-Za-z0-9_-]+$", run))
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
