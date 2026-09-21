#' Map a stage-group name to its sample-marker step id, if it has one
#'
#' Only stages with a plain in-process per-sample loop we control directly
#' get sample-level detail (see R/pipeline_sample_flags.R); everything else
#' (stages that hand off to ORFik's own internal BiocParallel dispatch, or
#' process a whole folder in one call) only ever gets study-level counts.
#' @param stage_name character, e.g. "pipe_trim_collapse"
#' @return character step id (e.g. "trim") or NA_character_
stage_marker_step <- function(stage_name) {
  switch(stage_name,
        pipe_trim_collapse = "trim",
        pipe_align_clean = "aligned",
        NA_character_)
}

#' Total sample count per experiment, from the live pipelines object
#'
#' Not derivable from flags alone -- needs the actual per-run metadata rows.
#' @param pipelines the pipelines list (see pipeline_init_all())
#' @return named integer vector, names are experiment ids
experiment_sample_counts <- function(pipelines) {
  counts <- list()
  for (p in pipelines) {
    for (organism in names(p$organisms)) {
      exp <- p$organisms[[organism]]$conf["exp"]
      counts[[exp]] <- nrow(p$study[ScientificName == organism])
    }
  }
  unlist(counts)
}

#' Nextflow-style stage checklist
#'
#' Read-only: reports on experiment-level flags (existing) and sample-level
#' markers (new, only available for stages with a plain per-sample loop --
#' see stage_marker_step()). Never downloads, mutates config, or blocks.
#'
#' @param pipelines the pipelines list
#' @param config the mNGSp config object
#' @param print logical, default TRUE. If TRUE, render via message() (not
#' cat() -- under MulticoreParam, cat() would write straight to the shared
#' real stdout and interleave with other concurrently-running stage-group
#' workers; message() is what BiocParallel's own logdir= already diverts
#' per-worker, same as the pre-existing "Sleep"/"Stopped sleeping" messages).
#' @return invisible(data.table) with columns: stage, done, total, state
#' ("done"/"running"/"queued"), active_experiment, active_done, active_total
#' (the latter three NA when no marker evidence is available/applicable)
pipeline_checklist <- function(pipelines, config, print = TRUE) {
  exps <- pipelines_names(pipelines)
  sample_totals <- experiment_sample_counts(pipelines)

  rows <- lapply(names(config$flag_steps), function(stage_name) {
    flags <- config$flag_steps[[stage_name]]
    done_vec <- vapply(exps, function(e) all(vapply(flags, step_is_done, logical(1),
                                                     config = config, experiment = e)),
                       logical(1))
    n_done <- sum(done_vec)
    n_total <- length(exps)

    marker_step <- stage_marker_step(stage_name)
    active_experiment <- NA_character_
    active_done <- NA_integer_
    active_total <- NA_integer_
    state <- if (n_done == n_total) "done" else "queued"

    if (n_done < n_total && !is.na(marker_step)) {
      candidates <- exps[!done_vec]
      started <- vapply(candidates, function(e) n_samples_done(config, marker_step, e) > 0, logical(1))
      if (any(started)) {
        active_experiment <- candidates[started][1]
        active_done <- n_samples_done(config, marker_step, active_experiment)
        active_total <- unname(sample_totals[active_experiment])
        state <- "running"
      }
    }

    data.table::data.table(stage = stage_name, done = n_done, total = n_total,
                           state = state, active_experiment = active_experiment,
                           active_done = active_done, active_total = active_total)
  })
  tab <- data.table::rbindlist(rows)
  if (print) message(format_checklist(tab))
  invisible(tab)
}

#' Format a checklist data.table as a compact printable table
#' @param tab data.table, as returned by pipeline_checklist(print = FALSE)
#' @return character, one formatted line per stage
format_checklist <- function(tab) {
  lines <- vapply(seq_len(nrow(tab)), function(i) {
    row <- tab[i]
    mark <- switch(row$state, done = "\u2714", running = ">", "-")
    detail <- if (!is.na(row$active_experiment)) {
      sprintf(" -- active: %s (%d/%d samples)", row$active_experiment, row$active_done, row$active_total)
    } else if (row$state == "running") {
      " -- active (no per-sample detail available for this stage)"
    } else ""
    sprintf("[%s] %-24s %d/%d studies done%s", mark, row$stage, row$done, row$total, detail)
  }, character(1))
  paste(lines, collapse = "\n")
}
