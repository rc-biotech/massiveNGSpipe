#' Search / define a list of accessions
#'
#' There are the primary ways:\cr
#' - Define the accessions yourself
#' - massive search over NCBI
#' - Get from rpfdb
#' @param user_specified character, default character(), user specified accessions
#' @param organisms character vector, default "all" (Use all organisms found).
#'  Else binomial latin name with capital letter for genus: "Homo sapiens" etc.
#' @param massive_search logical, default TRUE. Massive accession search over
#' NCBI using search terms defined over organisms and libtypes.
#' @param libtypes character, default "Ribosome profiling",
#' which library types to use: Other examples could be
#' "Guide seq", "mRNA seq", "CHIP seq", etc.
#' @param rpfdb logical, default "Ribosome profiling" %in% libtypes.
#' If so will get all rpfdb accessions for given organisms too.
#' @importFrom jsonlite fromJSON
#' @export
accessions_to_use <- function(user_specified = character(), organisms = "Homo sapiens",
                              massive_search = TRUE, libtypes = "Ribosome profiling",
                              rpfdb = "Ribosome profiling" %in% libtypes) {

  rpfdb_accessions <- search_accessions <- character()
  if (rpfdb) {
    rpfdb_studies <- as.data.table(jsonlite::fromJSON(
      system.file("extdata", "rpfdb-studies.json", package = "massiveNGSpipe")))
    organisms_rpfdb <- rpfdb_org_name(organisms)
    rpfdb_accessions <- trimws(rpfdb_studies[Species %in% organisms_rpfdb, ]$Study) # SacCer
  }

  if (massive_search) {
    message("- Massive NCBI search for:")
    for (org in organisms) {
      for (libtype in libtypes) {
        term <- paste(org, libtype)
        message("-- ", term)
        search_accessions <- c(search_accessions,
                               ORFik:::get_bioproject_candidates(term))
      }
    }
  }

  accessions <- c(user_specified, unique(c(search_accessions, rpfdb_accessions)))
  if (any(duplicated(accessions)))
    warning("Duplicated accessions found, you added user_specified accessions
            already found by either massive_search or rpfdb!")
  return(unique(c(accessions)))
}

rpfdb_org_name <- function(organisms) {
  non_standard_names <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                          "Drosophila melanogaster", "Arabidopsis thaliana",
                          "Danio rerio", "Cricetulus griseus")
  needs_renaming <- organisms %in% non_standard_names
  if (any(needs_renaming)) {
    rename_to_use <- non_standard_names %in% organisms
    organisms[needs_renaming] <- c("Human", "Mouse", "Rat", "Drosophila", "Arabidopsis",
                                   "Zebrafish", "Hamster")[rename_to_use]
  }
  name <- organisms[!needs_renaming]
  organisms[!needs_renaming] <- gsub(".* ", paste0(substr(name,1,1), "."), name)
  return(organisms)
}
