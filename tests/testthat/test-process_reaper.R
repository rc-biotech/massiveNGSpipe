zombie_children <- function() {
  st <- Sys.glob("/proc/[0-9]*/stat")
  sum(vapply(st, function(f) {
    s <- tryCatch(readLines(f, warn = FALSE), error = function(e) "")
    length(s) > 0 && grepl(paste0("^[0-9]+ \\(.*\\) Z ", Sys.getpid(), " "), s[1])
  }, logical(1)))
}

test_that("forks after a callr call are reaped, live callr children are kept", {
  skip_on_os("windows")
  j <- parallel::mcparallel(0); parallel::mccollect(j)   # parallel installs its SIGCHLD handler once
  long <- callr::r_bg(function() { Sys.sleep(4); 42 })   # processx replaces it
  p <- callr::r_bg(function() 1); p$wait()
  for (i in 1:5) { j <- parallel::mcparallel(i); parallel::mccollect(j) }
  expect_gt(zombie_children(), 0)   # documents the processx / parallel SIGCHLD conflict this guards against
  massiveNGSpipe:::reap_zombie_children(keep = long$get_pid())
  expect_equal(zombie_children(), 0)
  long$wait()
  expect_equal(long$get_result(), 42)
  expect_equal(long$get_exit_status(), 0)
})

test_that("run_forked leaves no zombies after callr use", {
  skip_on_os("windows")
  j <- parallel::mcparallel(0); parallel::mccollect(j)
  p <- callr::r_bg(function() 1); p$wait()
  res <- massiveNGSpipe:::run_forked(list(a = function() 1, b = function() 2))
  expect_equal(res, list(a = 1, b = 2))
  expect_equal(zombie_children(), 0)
})
