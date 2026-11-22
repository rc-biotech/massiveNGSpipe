# Children of this R session in the given state ("Z" = zombie); NULL = any state. The tests compare pid sets, not
# counts: in a full test run, children of earlier test files (BiocParallel workers, callr) can exit at any moment.
child_pids <- function(state = NULL) {
  me <- as.character(Sys.getpid())
  pids <- lapply(Sys.glob("/proc/[0-9]*/stat"), function(f) {
    s <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
    if (!length(s)) return(NULL)
    fl <- strsplit(sub("^.*\\) ", "", s[1]), " ", fixed = TRUE)[[1]]
    if (length(fl) < 2 || fl[2] != me || (!is.null(state) && fl[1] != state)) return(NULL)
    as.integer(sub(" .*", "", s[1]))
  })
  unlist(pids)
}
# mccollect() returns once the result is read, possibly before the child has exited: wait until it is a zombie or gone
wait_exited <- function(pids, timeout = 5) {
  t0 <- Sys.time()
  while (length(setdiff(intersect(pids, child_pids()), child_pids("Z"))) && Sys.time() - t0 < timeout) Sys.sleep(0.05)
}

test_that("forks after a callr call are reaped, live callr children are kept", {
  skip_on_os("windows")
  j <- parallel::mcparallel(0); parallel::mccollect(j)   # parallel installs its SIGCHLD handler once
  long <- callr::r_bg(function() { Sys.sleep(4); 42 })   # processx replaces it
  p <- callr::r_bg(function() 1); p$wait()
  forks <- integer()
  for (i in 1:5) { j <- parallel::mcparallel(i); forks <- c(forks, j$pid); parallel::mccollect(j) }
  wait_exited(forks)
  expect_gt(length(intersect(forks, child_pids("Z"))), 0)   # documents the processx / parallel SIGCHLD conflict this guards against
  massiveNGSpipe:::reap_zombie_children(keep = long$get_pid())
  expect_length(intersect(forks, child_pids()), 0)
  long$wait()
  expect_equal(long$get_result(), 42)
  expect_equal(long$get_exit_status(), 0)
})

test_that("run_forked leaves no zombies after callr use", {
  skip_on_os("windows")
  j <- parallel::mcparallel(0); parallel::mccollect(j)
  p <- callr::r_bg(function() 1); p$wait()
  before <- child_pids()
  res <- massiveNGSpipe:::run_forked(list(a = function() 1, b = function() 2))
  expect_equal(res, list(a = 1, b = 2))
  expect_length(setdiff(child_pids("Z"), before), 0)   # no zombie left by run_forked's own forks
})
