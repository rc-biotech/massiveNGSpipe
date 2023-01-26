pipeline_merge_study <- function(df_list, config) {
  for (df in df_list) {
    if (!step_is_next_not_done(config, "merged_lib", name(df))) next
    ORFik::mergeLibs(df, file.path(libFolder(df), "pshifted_merged"),
                     "lib", "pshifted", FALSE)
    set_flag(config, "merged_lib", name(df))
  }
}
