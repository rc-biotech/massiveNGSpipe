# Per-Run, per-mapper offset binding. Legacy outputs without a record retain
# their historical resume rule; only newly applied offsets are certified here.
# Explicit requests bridge table edits to the regular pipeline without changing
# raw OFST mtimes or deleting P-site files.
pshift_canonical_offsets <- function(x) {
  if (!is.data.frame(x) || !all(c("fraction", "offsets_start") %in% names(x)) || !nrow(x))
    stop("P-site offsets require nonempty fraction/offsets_start columns")
  f <- as.numeric(x$fraction); o <- as.numeric(x$offsets_start)
  if (anyNA(f) || anyNA(o) || any(!is.finite(f) | !is.finite(o)) ||
      any(f != as.integer(f) | o != as.integer(o) | f <= 0 | o > 0 | -o >= f) || anyDuplicated(f))
    stop("P-site offsets require unique positive integer lengths and valid nonpositive integer offsets")
  i <- order(f)
  data.frame(fraction = as.integer(f[i]), offsets_start = as.integer(o[i]))
}

pshift_file_identity <- function(path) {
  if (length(path) != 1L || !file.exists(path)) return(NULL)
  z <- file.info(path)
  list(path = normalizePath(path, mustWork = TRUE), size = unname(z$size),
       mtime = as.numeric(z$mtime), ctime = as.numeric(z$ctime))
}

pshift_offset_record_paths <- function(raw, shifted) {
  id <- sub("\\.ofst$", "", basename(raw))
  list(pending = file.path(dirname(shifted), ".offsets_pending", paste0(id, ".rds")),
       applied = file.path(dirname(shifted), ".offset_provenance", paste0(id, ".rds")))
}

pshift_atomic_record <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path)); on.exit(unlink(tmp))
  saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish P-site provenance: ", path)
  invisible(NULL)
}

pshift_table_entry <- function(raw, table_path) {
  if (!file.exists(table_path)) return(NULL)
  tab <- readRDS(table_path)
  keys <- normalizePath(names(tab), mustWork = FALSE)
  hit <- which(keys == normalizePath(raw, mustWork = FALSE))
  if (length(hit) != 1L) return(NULL)
  pshift_canonical_offsets(tab[[hit]])
}

# Called before an offset-table edit. A crash between marker and table write is
# safe: the normal pipeline cannot mark that request applied against old offsets.
pshift_mark_offsets_pending <- function(table_path, shifts, reason = "") {
  if (!is.list(shifts) || !length(shifts) || is.null(names(shifts)) ||
      anyNA(names(shifts)) || anyDuplicated(names(shifts))) stop("Offsets must be named by exact raw OFST paths")
  entries <- lapply(shifts, pshift_canonical_offsets) # validate all before writing
  for (i in seq_along(entries)) {
    raw <- names(entries)[i]
    if (!file.exists(raw)) stop("Raw OFST missing for requested offsets: ", raw)
    shifted <- file.path(dirname(table_path), paste0(sub("\\.ofst$", "", basename(raw)), "_pshifted.ofst"))
    paths <- pshift_offset_record_paths(raw, shifted)
    pshift_atomic_record(list(schema = 1L, raw = pshift_file_identity(raw),
      offsets = entries[[i]], previous_shifted = pshift_file_identity(shifted),
      reason = reason, requested = Sys.time()), paths$pending)
  }
  invisible(NULL)
}

# Internal staging API for a standard pipeline rerun. It does not change flags,
# run shifting, or certify the existing P-site file.
pshift_stage_offsets <- function(df, shifts, reason = "") {
  raw <- ofst_filepath(df)
  if (!all(names(shifts) %in% raw)) stop("Requested offset keys do not match this mapper branch")
  table_path <- file.path(dirname(pshifted_filepath(df[1, ])), "shifting_table.rds")
  pshift_mark_offsets_pending(table_path, shifts, reason)
  ORFik::shifts_save(shifts, dirname(table_path))
  invisible(NULL)
}

# Cheap recipe/source guards: never infer biological suitability from a saved
# table's existence. Missing legacy provenance permits historical resume only,
# not reuse of offsets for a product that needs rebuilding.
pshift_detection_recipe <- function(accepted_lengths, allowed_hard12_species,
                                    max_no_adapter_removed_pct) {
  list(version = "raw-bound-v1", accepted_lengths = sort(unique(as.integer(accepted_lengths))),
       allowed_hard12_species = sort(unique(as.character(allowed_hard12_species))),
       max_no_adapter_removed_pct = as.numeric(max_no_adapter_removed_pct))
}

pshift_recipe_valid <- function(df_one_row, recipe) {
  paths <- pshift_offset_record_paths(ofst_filepath(df_one_row), pshifted_filepath(df_one_row))
  if (!file.exists(paths$applied)) return(TRUE) # legacy resume, not certification
  tryCatch({
    record <- readRDS(paths$applied)
    is.null(record$recipe) || identical(record$recipe, recipe)
  }, error = function(e) FALSE)
}

# Explicit requests override the reuse switch, but must match the current raw
# source and exact saved entry. Corruption/stale requests stop before any writer.
pshift_reusable_offsets <- function(df_one_row, recipe, reuse = TRUE) {
  raw <- ofst_filepath(df_one_row); shifted <- pshifted_filepath(df_one_row)
  paths <- pshift_offset_record_paths(raw, shifted)
  table_path <- file.path(dirname(shifted), "shifting_table.rds")
  if (file.exists(paths$pending)) {
    pending <- readRDS(paths$pending)
    entry <- pshift_table_entry(raw, table_path)
    if (!identical(pending$schema, 1L) || is.null(entry) ||
        !identical(pending$raw, pshift_file_identity(raw)) ||
        !identical(pending$offsets, entry))
      stop("Pending P-site offsets/source differ from explicit request: ", raw)
    return(entry)
  }
  if (!isTRUE(reuse) || !file.exists(paths$applied)) return(NULL)
  tryCatch({
    record <- readRDS(paths$applied)
    entry <- pshift_table_entry(raw, table_path)
    if (identical(record$schema, 1L) && !is.null(entry) &&
        identical(record$raw, pshift_file_identity(raw)) &&
        identical(record$recipe, recipe) && identical(record$offsets, entry)) entry else NULL
  }, error = function(e) NULL)
}

pshift_offset_binding_valid <- function(df_one_row) {
  raw <- ofst_filepath(df_one_row); shifted <- pshifted_filepath(df_one_row)
  paths <- pshift_offset_record_paths(raw, shifted)
  if (file.exists(paths$pending)) return(FALSE)
  if (!file.exists(paths$applied)) return(TRUE) # legacy: not retroactively certified
  tryCatch({
    record <- readRDS(paths$applied)
    current <- pshift_table_entry(raw, file.path(dirname(shifted), "shifting_table.rds"))
    identical(record$schema, 1L) && !is.null(current) &&
      identical(record$raw, pshift_file_identity(raw)) &&
      identical(record$shifted, pshift_file_identity(shifted)) &&
      identical(record$offsets, current)
  }, error = function(e) FALSE)
}

# Called only after the native shifting writer succeeds. Pending requests are
# cleared only for existing products whose saved offsets match the request.
pshift_record_applied_offsets <- function(df, recipe = NULL) {
  for (i in seq_len(nrow(df))) {
    one <- df[i, ]; raw <- ofst_filepath(one); shifted <- pshifted_filepath(one)
    paths <- pshift_offset_record_paths(raw, shifted)
    pending <- if (file.exists(paths$pending)) readRDS(paths$pending) else NULL
    table_path <- file.path(dirname(shifted), "shifting_table.rds")
    current <- pshift_table_entry(raw, table_path)
    if (!file.exists(raw) || !file.exists(shifted) || is.null(current)) {
      if (!is.null(pending)) stop("Requested P-site output or saved table missing: ", raw)
      next
    }
    if (!is.null(pending) && (!identical(pending$raw, pshift_file_identity(raw)) ||
                             !identical(pending$offsets, current)))
      stop("Applied P-site offsets/source differ from explicit request: ", raw)
    if (!is.null(pending) && identical(pending$previous_shifted, pshift_file_identity(shifted)))
      stop("P-site product was not rewritten after explicit offset request: ", raw)
    # Extra-length writers retain a same-source recipe; legacy records do not
    # acquire a fresh-detection certificate merely by being read again.
    bound_recipe <- recipe
    if (is.null(bound_recipe) && file.exists(paths$applied)) {
      old <- readRDS(paths$applied)
      if (identical(old$raw, pshift_file_identity(raw))) bound_recipe <- old$recipe
    }
    pshift_atomic_record(list(schema = 1L, raw = pshift_file_identity(raw),
      shifted = pshift_file_identity(shifted), offsets = current, recipe = bound_recipe,
      offset_source = if (!is.null(pending)) "explicit_request" else "native_pipeline",
      applied = Sys.time()), paths$applied)
    if (!is.null(pending)) {
      if (!identical(pending, readRDS(paths$pending))) stop("P-site request changed during apply: ", raw)
      unlink(paths$pending)
    }
  }
  invisible(NULL)
}
