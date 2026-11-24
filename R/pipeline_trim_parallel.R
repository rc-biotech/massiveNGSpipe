#' Trim several libraries of one experiment at once
#'
#' The trim stage does two different kinds of work per library:
#' \itemize{
#'  \item detection (\code{trim_sample_detect()}): adapter scan, 1 M-read detection sample, barcode/UMI calls and the
#'   adapter check on 1 M raw reads. Reads only the start of the file; mostly single-threaded R (~1 min per library).
#'  \item the trim pass (\code{trim_sample_finish()}): one full read of the raw file by the trimcollapse binary (or
#'   fastp), which already streams internally (reader thread -> trim workers -> collapse shards -> pigz output). It is
#'   bound by its single decompression thread and uses ~5-6 cores, whatever the thread count above 8.
#' }
#' \code{method = "serial"} runs detect + trim per library in order, exactly as before. \code{"staged"} (default) forks
#' detection jobs and trim jobs from this (already loaded) process, with separate slot counts: as soon as a library's
#' detection finishes its trim pass is started, while the next libraries are being detected, so the read-, CPU- and
#' write-heavy phases of different libraries overlap. Forking avoids the package loading of new R sessions (~16-35 s).
#' \code{"naive"} forks \code{workers} whole-library jobs (detect + trim) at once, for comparison.
#'
#' Every job writes only its own run's files (trim/report_<run>.json, trim/SINGLE/collapsed_trimmed_<run>.fasta.gz,
#' barcode_detect_sample/<run>, its sample markers); the one study-wide file written during trimming
#' (mate_override_auto.csv) is written under a file lock. On an error no new library is started, the running ones
#' finish (and get their markers), then the stage fails with all errors, so the stage flag is not set, as in the serial
#' path.
#' @param ctx list(all_files, runs, mode, target_dir, trimmed_dir, source_dir, config, experiment)
#' @param method "serial", "staged" or "naive"
#' @param workers NULL (dynamic, \code{trim_workers_dynamic()}) or a fixed number of parallel libraries
#' @return list of barcode_dt rows, one per run
#' @noRd
trim_samples_run <- function(ctx, method = "serial", workers = NULL) {
  n <- nrow(ctx$runs)
  method <- match.arg(method, c("serial", "staged", "naive"))
  if (is.null(workers)) workers <- trim_workers_dynamic(n)
  workers <- max(1L, min(as.integer(workers), n))
  if (method == "serial" || workers == 1L || n == 1L || .Platform$OS.type == "windows")
    return(lapply(seq_len(n), function(i) trim_sample_finish(ctx, i, trim_sample_detect(ctx, i))))
  # one-time setup shared by all jobs, before forking: the adapter candidate file (written on first use) and the
  # trimcollapse binary (compiled on first use) must not be created by several jobs at once
  invisible(adapter_list())
  if (!isFALSE(getOption("massiveNGSpipe.fast_trim", TRUE))) try(trimcollapse_binary(), silent = TRUE)
  slots <- trim_slots(workers, ctx$config)
  # largest raw input first: the longest trim pass then overlaps all smaller libraries instead of running alone at
  # the end (a 14 GB library started last made the parallel run no faster than serial)
  sizes <- vapply(ctx$all_files, function(f) sum(file.size(f[!is.na(f)]), na.rm = TRUE), numeric(1))
  start_order <- order(-sizes, seq_len(n))
  message("Parallel trim (", method, "): ", n, " libraries, ", workers, " at once",
          if (method == "staged") paste0(" (", slots[["detect"]], " detecting, ", slots[["trim"]], " trimming)"))
  if (method == "naive") {
    return(run_fork_pool(n, list(all = function(i) trim_sample_finish(ctx, i, trim_sample_detect(ctx, i))),
                         stages = "all", slots = c(all = workers), run_ids = ctx$runs$Run, start_order = start_order))
  }
  run_fork_pool(n, list(detect = function(i) trim_sample_detect(ctx, i),
                        trim = function(i, det) trim_sample_finish(ctx, i, det)),
                stages = c("detect", "trim"), slots = slots, run_ids = ctx$runs$Run, start_order = start_order,
                # at most 1.5 x workers jobs: trims first, new detections fill the rest
                max_jobs = max(2L, as.integer(ceiling(1.5 * workers))))
}

#' Detection and trim slots of the staged trim pool for a number of parallel libraries
#'
#' Detection is mostly single-threaded R (~1 min per library, whatever its size), a trim pass ~5-6 cores for ~0.2-10
#' min (scales with the file): \code{workers} detection slots and as many trim slots, so a detected library never waits
#' for a trim slot (with half as many trim slots, two small libraries detected first held both slots and the largest
#' library, whose trim pass sets the stage's wall time, waited 3 min for one). A fixed
#' \code{config$trim_stage_slots} or \code{options(massiveNGSpipe.trim_stage_slots = c(detect = d, trim = t))} wins.
#' @noRd
trim_slots <- function(workers, config = list()) {
  fixed <- config[["trim_stage_slots", exact = TRUE]]
  if (is.null(fixed)) fixed <- getOption("massiveNGSpipe.trim_stage_slots", NULL)
  if (!is.null(fixed)) return(c(detect = max(1L, as.integer(fixed[["detect"]])), trim = max(1L, as.integer(fixed[["trim"]]))))
  workers <- max(1L, as.integer(workers))
  c(detect = workers, trim = workers)
}

#' Detection half of one library's trim: everything that only reads the start of the raw file
#' @param ctx see trim_samples_run()
#' @param i row of ctx$runs
#' @return list(kind, ...) for trim_sample_finish()
#' @noRd
trim_sample_detect <- function(ctx, i) {
  study_sample <- ctx$runs[i]
  run <- study_sample$Run
  message(run)
  filenames <- ctx$all_files[[i]]
  assert_base_space_fastq(filenames[!is.na(filenames)], run)
  single_end <- is.na(filenames[2])
  file <- filenames[1]
  file2 <- if (!single_end) filenames[2]
  # Opt-in: trim + collapse only the chosen mate (mate_override_manual.csv, R/mate_override.R)
  if (run %in% mate_overrides(ctx$trimmed_dir)$Run) return(list(kind = "mate_override"))
  # Paired footprint libraries: only the mate with the insert is trimmed, collapsed and aligned, as
  # single-end (read 1, or read 2 when read 1 has no insert); recorded in mate_override_auto.csv
  if (!single_end && study_sample$LIBRARYTYPE %in% footprint_library_types() && !isFALSE(ctx$config$auto_footprint_mate)) {
    decision <- paired_footprint_mate(file, file2)
    message(run, ": footprint mate ", decision$mate, " (", decision$reason, ")")
    return(list(kind = "footprint", decision = decision))
  }
  check_for_barcodes <- study_sample$LIBRARYTYPE == "RFP"
  # RNA inputs built with small-RNA kits may require explicit
  # random-end crops. Honor reviewed manual barcode rows without
  # enabling automatic barcode guesses for ordinary RNA libraries.
  manual_barcodes <- file.path(ctx$trimmed_dir, "barcodes_manual.csv")
  if (!check_for_barcodes && file.exists(manual_barcodes)) {
    manual <- data.table::fread(manual_barcodes)
    check_for_barcodes <- "Run" %in% names(manual) && run %in% manual$Run
  }
  if (check_for_barcodes) {
    # Detects barcode sizes from a cheap subsample; the trim pass then applies front/adapter then post-adapter tail
    # cuts (run_barcode_detection_and_trim()).
    return(list(kind = "barcode",
                barcode_dt = barcode_detection_for_trim(study_sample, ctx$source_dir, ctx$target_dir, ctx$trimmed_dir,
                                                        ctx$mode)))
  }
  list(kind = "adapter", resolved = resolve_adapter_for_trim(file, ctx$target_dir))
}

#' Trim half of one library: the full-file trim (+ collapse) pass, sample markers and read orientation
#' @param det result of trim_sample_detect()
#' @return the library's barcode_dt row (also stored in its "trim" marker)
#' @noRd
trim_sample_finish <- function(ctx, i, det) {
  study_sample <- ctx$runs[i]
  run <- study_sample$Run
  config <- ctx$config; experiment <- ctx$experiment; trimmed_dir <- ctx$trimmed_dir; target_dir <- ctx$target_dir
  filenames <- ctx$all_files[[i]]
  single_end <- is.na(filenames[2])
  file <- filenames[1]
  file2 <- if (!single_end) filenames[2]
  fast_trim_state$collapsed <- NULL
  barcode_dt <- data.table()
  if (det$kind == "mate_override") {
    barcode_dt <- trim_with_mate_override(run, filenames, target_dir, trimmed_dir)
    set_sample_flag(config, "trim", experiment, run, value = barcode_dt)
    set_sample_flag(config, "collapsed", experiment, run)
    return(barcode_dt)
  }
  if (det$kind == "footprint") {
    decision <- det$decision
    barcode_dt <- trim_from_footprint_mate(study_sample, if (decision$mate == 2L) file2 else file,
                                           target_dir, trimmed_dir, ctx$mode, decision)
    set_sample_flag(config, "trim", experiment, run, value = barcode_dt)
    if (!is.null(fast_trim_state$collapsed)) {
      if (decision$mate == 1L) orient_collapsed_fasta(fast_trim_state$collapsed, run, trimmed_dir)
      set_sample_flag(config, "collapsed", experiment, run)
    }
    return(barcode_dt)
  }
  if (det$kind == "barcode") {
    barcode_dt <- barcode_trim_pass(det$barcode_dt, study_sample, target_dir, trimmed_dir, file, file2)
  } else {
    adapter_trim_pass(det$resolved, file, target_dir, file2)
  }
  # Store the actual barcode_dt row (not just TRUE) so it
  # can be reconstructed below even across a resumed run.
  set_sample_flag(config, "trim", experiment, run, value = barcode_dt)
  # The fast path also collapsed the run: finish it like pipeline_collapse() would (orientation, marker)
  if (!is.null(fast_trim_state$collapsed)) {
    orient_collapsed_fasta(fast_trim_state$collapsed, run, trimmed_dir)
    set_sample_flag(config, "collapsed", experiment, run)
  }
  barcode_dt
}

#' Run per-item stages in forked jobs with a slot limit per stage
#'
#' Item i goes through \code{stages} in order; stage k > 1 receives the result of stage k - 1. A free slot of a later
#' stage is filled first (finish started libraries before starting new ones). Jobs are forks of this process: no
#' package loading. Each job's messages are prefixed with its run id. After the first failure no new item is started;
#' running jobs finish, then all errors are raised together.
#' @param n number of items
#' @param funs named list of functions: the first takes (i), later ones (i, previous_result)
#' @param stages names of funs in order
#' @param slots named integer, concurrent jobs per stage
#' @param run_ids labels for messages
#' @param poll seconds between checks
#' @param start_order order in which items enter the first stage (later stages: in order of readiness)
#' @param max_jobs cap on running jobs over all stages (later stages are filled first)
#' @return list of the last stage's results, by item
#' @noRd
run_fork_pool <- function(n, funs, stages = names(funs), slots, run_ids = as.character(seq_len(n)), poll = 0.5,
                          start_order = seq_len(n), max_jobs = Inf) {
  stage_of <- integer(n)                           # last finished stage per item (0 = none)
  value <- vector("list", n)
  running <- list()                                # each: list(job, i, stage)
  errors <- character()
  start_job <- function(i, k) {
    f <- funs[[stages[k]]]; prev <- value[[i]]; label <- run_ids[i]
    job <- parallel::mcparallel({
      withCallingHandlers(
        list(if (k == 1L) f(i) else f(i, prev)),
        message = function(m) {
          m$message <- paste0("[", label, " ", format(Sys.time(), "%H:%M:%S"), "] ", m$message)
          message(m); invokeRestart("muffleMessage")
        })
    }, mc.set.seed = FALSE, silent = FALSE)
    running[[length(running) + 1L]] <<- list(job = job, i = i, stage = k)
  }
  next_item <- 1L
  repeat {
    # collect finished jobs
    if (length(running)) {
      done <- parallel::mccollect(lapply(running, `[[`, "job"), wait = FALSE, timeout = poll)
      if (!is.null(done)) {
        pids <- vapply(running, function(r) r$job$pid, integer(1))
        for (pid in names(done)) {
          idx <- which(pids == as.integer(pid))
          if (!length(idx)) next
          r <- running[[idx]]
          tryCatch(suppressWarnings(sys::exec_status(r$job$pid, wait = TRUE)), error = function(e) NULL)
          res <- done[[pid]]
          if (is.null(res)) {
            errors <- c(errors, paste0(run_ids[r$i], " (", stages[r$stage], "): job died without a result (killed or out of memory?)"))
          } else if (inherits(res, "try-error")) {
            errors <- c(errors, paste0(run_ids[r$i], " (", stages[r$stage], "): ", trimws(as.character(res))))
          } else {
            value[r$i] <- list(res[[1]]); stage_of[r$i] <- r$stage
          }
          running <- running[-idx]; pids <- pids[-idx]
        }
      }
    }
    # start jobs: later stages first, then new items, within each stage's free slots
    if (!length(errors)) {
      for (k in rev(seq_along(stages))) {
        busy <- sum(vapply(running, function(r) r$stage == k, logical(1)))
        free <- min(slots[[stages[k]]] - busy, max_jobs - length(running))
        if (free <= 0) next
        if (k == 1L) {
          while (free > 0 && next_item <= n) { start_job(start_order[next_item], 1L); next_item <- next_item + 1L; free <- free - 1L }
        } else {
          in_flight <- vapply(running, `[[`, integer(1), "i")
          ready <- start_order[stage_of[start_order] == k - 1L & !(start_order %in% in_flight)]
          for (i in utils::head(ready, free)) start_job(i, k)
        }
      }
    }
    # nothing running after the start step = all items done, or a failure stopped new starts
    if (!length(running)) break
  }
  suppressWarnings(reap_zombie_children())   # a child exiting during the /proc scan only warns
  if (length(errors)) stop("Parallel trim failed for ", length(errors), " job(s):\n", paste(errors, collapse = "\n"),
                           call. = FALSE)
  value
}

#' Trim method from the config or options: "staged" (default), "serial" or "naive"
#' @noRd
trim_parallel_method <- function(config = list()) {
  m <- config[["trim_parallel", exact = TRUE]]
  if (is.null(m)) m <- getOption("massiveNGSpipe.trim_parallel", "staged")
  if (isFALSE(m)) m <- "serial"
  if (isTRUE(m)) m <- "staged"
  match.arg(as.character(m), c("serial", "staged", "naive"))
}

#' Fixed number of libraries trimmed at once, from the config or options (NULL = dynamic)
#' @noRd
trim_workers_fixed <- function(config = list()) {
  w <- config[["trim_workers", exact = TRUE]]
  if (is.null(w)) w <- getOption("massiveNGSpipe.trim_workers", NULL)
  if (is.null(w)) return(NULL)
  max(1L, as.integer(w))
}

#' Number of libraries to trim at once on this host now
#'
#' A fixed \code{options(massiveNGSpipe.trim_workers = n)} (or \code{pipeline_config(trim_workers)}, passed in by the
#' caller) wins. Otherwise: idle CPUs (cores minus the 1-min load average) at
#' \code{massiveNGSpipe.trim_cores_per_worker} (default 4: detection ~1.5 cores, a trim pass ~5-6 cores; one library is
#' in one of the two at a time, and short detections overlap long trims), capped
#' by free memory above the reserve (MemAvailable - \code{massiveNGSpipe.mem_reserve_gb}, default 80) at
#' \code{massiveNGSpipe.trim_gb_per_worker} (default 8 GB: peak PSS 26-30 GB for 4 libraries at once, mostly the 1 M-read
#' detection samples and their 8 adapter-search forks), with a floor of \code{massiveNGSpipe.trim_min_workers}
#' (default 2) on a busy host unless memory is short, bounded to \code{massiveNGSpipe.trim_max_workers} (default 8)
#' and the number of libraries. Where /proc is unavailable: 1 (serial).
#' @param n_samples number of libraries left to trim
#' @param cores,load1,mem_gb host state (read from /proc when NULL; arguments for tests)
#' @return integer >= 1
#' @noRd
trim_workers_dynamic <- function(n_samples = Inf, cores = NULL, load1 = NULL, mem_gb = NULL) {
  fixed <- getOption("massiveNGSpipe.trim_workers", NULL)
  snap <- resource_snapshot()
  if (!is.null(cores)) snap$cores <- cores
  if (!is.null(load1)) snap$load1 <- load1
  if (!is.null(mem_gb)) snap$mem_available_gb <- mem_gb
  if (!is.null(fixed))   # a fixed count is not load-based; only an active resource mode caps it
    return(resource_threads("trim_workers", max(1L, as.integer(fixed)), cores_per_worker = 0, per_worker_gb = NA,
                            snapshot = snap))
  if (length(snap$load1) != 1 || length(snap$mem_available_gb) != 1) return(1L)
  # one formula (R/resources.R), parameters as before + the resource mode cap
  as.integer(resource_threads("trim_workers", getOption("massiveNGSpipe.trim_max_workers", 8L), n = n_samples,
                              cores_per_worker = getOption("massiveNGSpipe.trim_cores_per_worker", 4),
                              per_worker_gb = getOption("massiveNGSpipe.trim_gb_per_worker", 8),
                              min = getOption("massiveNGSpipe.trim_min_workers", 2L), snapshot = snap))
}
