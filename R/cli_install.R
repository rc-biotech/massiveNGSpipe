#' Put the mngsp command on the PATH
#'
#' Symlinks the package's \code{bin/mngsp} wrapper into \code{dest}. The link points into the installed package, so
#' it follows package updates. The wrapper unsets \code{LC_ALL}, raises the open-file limit (STAR) up to the hard
#' limit and the C stack to 64 MB, then calls \code{\link{mngsp_cli}}.
#' @param dest folder on the PATH, default "~/.local/bin"
#' @param name command name, default "mngsp"
#' @param overwrite logical, replace an existing file or link
#' @return the link path, invisibly
#' @export
#' @examples
#' \dontrun{
#' install_cli()                         # ~/.local/bin/mngsp
#' install_cli("/usr/local/bin")         # system-wide (needs write access)
#' }
install_cli <- function(dest = "~/.local/bin", name = "mngsp", overwrite = FALSE) {
  src <- system.file("bin", "mngsp", package = "massiveNGSpipe")
  if (!nzchar(src)) stop("bin/mngsp not found in the massiveNGSpipe installation; reinstall the package")
  if (file.access(src, 1) != 0 && !isTRUE(Sys.chmod(src, "0755", use_umask = FALSE)))
    stop(src, " is not executable and cannot be made executable (chmod +x it)")
  dest <- path.expand(dest)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  link <- file.path(dest, name)
  if (file.exists(link) || nzchar(Sys.readlink(link))) {
    if (identical(normalizePath(link, mustWork = FALSE), normalizePath(src))) {
      message(link, " already points to ", src); return(invisible(link))
    }
    if (!overwrite) stop(link, " exists; use overwrite = TRUE to replace it")
    unlink(link)
  }
  if (!file.symlink(src, link)) stop("could not create ", link)
  on_path <- normalizePath(dest) %in% normalizePath(strsplit(Sys.getenv("PATH"), ":")[[1]], mustWork = FALSE)
  message("Installed ", link, " -> ", src,
          if (!on_path) paste0("\n", dest, " is not on your PATH; add to ~/.bashrc:  export PATH=\"", dest, ":$PATH\""),
          "\nTry: mngsp --help")
  invisible(link)
}
