# base64url (no padding) encoding for raw vectors
.base64url <- function(raw_bytes) {
  s <- openssl::base64_encode(raw_bytes)
  s <- sub("=+$", "", s)                 # strip '=' padding
  chartr("+/", "-_", s)                  # URL-safe alphabet
}

# GA4GH-style sha512t24u: SHA-512, first 24 bytes, base64url
sha512t24u <- function(x) {
  dt <- data.table(str = x)
  dt[, enc := .base64url(digest::digest(str, algo = "sha512", serialize = FALSE, raw = TRUE)[seq_len(24L)]), by = str]$enc
}

# Length bucket Lk with decade-centered bins (10^(k-0.5) .. 10^(k+0.5))
length_bucket <- function(total_len_nt) {
  if (any(!is.finite(total_len_nt) | total_len_nt <= 0)) {
    stop("Non-positive total length for translon.")
  }
  k <- round(log10(total_len_nt))
  paste0("L", k)
}

#' Convert GRangesList of translons to generic strings
#'
#' This is input for make_translon_id
#' @param grl a GRangesList
#' @param assembly character, the assembly name (e.g. GRCh38)
#' @return a character vector
#' @examples
#' grl <- makeGRangesListFromCharacter("chr22:101231230-101231500:+;chr22:101240000-101240962:+")
#' genome(grl) <- "GRCh38"
#' translon_string <- as_translon_info_list(grl)
#' translon_string
as_translon_info_list <- function(grl, assembly = unique(genome(grl))) {
  stopifnot(is(grl, "GRangesList"))
  stopifnot(!is.na(assembly) && !is.null(assembly) && is.character(assembly) &&
              length(assembly) == 1)

  chr <- seqnamesPerGroup(grl, keep.names = FALSE)
  strand <- strandPerGroup(grl, keep.names = FALSE)
  width <- widthPerGroup(grl, keep.names = FALSE)
  region_string <- as.character(ranges(grl))

  return(list(asm = assembly, ctg = chr,
              strand = strand, total_len = width,
              region_string = region_string))
}

setMethod("as.character", "IRangesList", function(x, ...) {
  if (length(x) == 0) return(character())
  u <- unlist(x, use.names = FALSE)
  per <- paste0(start(u), "-", end(u))
  cl <- relist(per, x)
  res <- unstrsplit(cl, sep = ";")
  # names(res) <- names(x)
  return(res)
})

#' Make translon id
#'
#' Make the readable, stable ID:
#' \code{TRN:asm:ctg:Lk:strand:short-hash}
#' where 'short-hash' = first 'short_len' chars of
#' \code{sha512t24u(canonical_string$region_string)}, i.e string of
#' string of genomic coordinates without chromosome or strand.
#' @importFrom digest digest
#' @importFrom openssl base64_encode
#' @param translon_info_list a list translon information needed of length 5.
#' @param short_len integer, default 10
#' @return character vector of length equal to \code{lengths(translon_info_list)[2]}
#' @export
#' @examples
#' grl <- makeGRangesListFromCharacter("chr22:101231230-101231500:+;chr22:101240000-101240962:+")
#' genome(grl) <- "GRCh38"
#' as_translon_info_list <- as_translon_info_list(grl)
#' as_translon_info_list
#' make_translon_id(as_translon_info_list)
#' # output: "TRN:GRCh38:chr22:L3:+:AB3FZK2HQY"   # (your hash will differ)
#' make_translon_machine_id(canon)
#' # output: "TRN:V1:sha512t24u:Oqk...<longer url-safe digest>"
make_translon_id <- function(translon_info_list, short_len = 10) {
  meta <- translon_info_list
  Lk   <- length_bucket(meta$total_len)

  # computed identifier over the EXACT canonical string (assembly|ctg|strand|blocks)
  full_digest <- sha512t24u(meta$region_string)

  # For the readable suffix, take the first N chars (tweak N if you want more entropy)
  short_hash <- substr(full_digest, 1L, short_len)

  paste("TRN", meta$asm, meta$ctg, Lk, meta$strand, short_hash, sep = ":")
}

# (Optional) also expose the authoritative machine ID if you want to store it:
# TRN:V1:sha512t24u:<digest>
make_translon_machine_id <- function(canonical_string) {
  paste0("TRN:V1:sha512t24u:", sha512t24u(canonical_string))
}

