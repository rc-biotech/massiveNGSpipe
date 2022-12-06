pipeline_merge_per_study <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "merged_lib", name(df))) next
    all_merged <- outputLibs(df, type = "pshifted", output.mode = "list")
    names(all_merged) <- NULL
    libs <- bamVarName(df)
    message("- Merging..")
    final_merged <- collapseDuplicatedReads(do.call("c", all_merged), addSizeColumn = TRUE)
    merged_folder <- file.path(libFolder(df), "pshifted_merged")
    ofst.path <- paste0(merged_folder, "/", "RFP", "_merged.ofst")
    dir.create(merged_folder, showWarnings = FALSE)
    export.ofst(final_merged, ofst.path)
    remove.experiments(df)
    set_flag(config, "merged_lib", name(df)) 
  }
}