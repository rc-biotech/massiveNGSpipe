# Unit tests for R/star_index_lock.R -- the cross-process claim
# mechanism that stops one process from removing an organism's STAR
# shared-memory index while a DIFFERENT process is still using it
# (built from a real ~70-minute silent stall this session,
# PRJNA1051101-homo_sapiens, 2026-10-08).

test_that("claim then release leaves no active claims behind", {
  index <- tempfile("star_index_")
  dir.create(index, recursive = TRUE)

  owner <- star_index_claim(index)
  expect_true(file.exists(file.path(index, ".align_claims", paste0(owner, ".rds"))))
  expect_true(star_index_safe_to_remove(index, exclude_owner = owner))

  star_index_release(index, owner)
  expect_false(file.exists(file.path(index, ".align_claims", paste0(owner, ".rds"))))
})

test_that("an active claim from ANOTHER owner blocks removal, and never counts against itself", {
  index <- tempfile("star_index_")
  dir.create(index, recursive = TRUE)

  mine <- star_index_claim(index, owner_id = "mine")
  # A second, genuinely different owner, same live PID (this test process) --
  # simulates a sibling organism in the SAME run rather than a different
  # process, which is enough to exercise the "someone else is active" path.
  other <- star_index_claim(index, owner_id = "other")

  expect_false(star_index_safe_to_remove(index, exclude_owner = mine))
  expect_identical(star_index_active_claims(index, exclude_owner = mine), "other")
  expect_identical(star_index_active_claims(index, exclude_owner = other), "mine")

  star_index_release(index, mine)
  star_index_release(index, other)
})

test_that("a claim left by a dead PID is purged automatically and does not block removal", {
  index <- tempfile("star_index_")
  dir.create(index, recursive = TRUE)

  d <- star_index_lock_dir(index)
  dir.create(d, recursive = TRUE)
  # An implausibly large PID almost certainly not in use.
  saveRDS(2147483000L, file.path(d, "stale_owner.rds"))

  expect_true(star_index_safe_to_remove(index, exclude_owner = "someone_else"))
  expect_false(file.exists(file.path(d, "stale_owner.rds"))) # purged as a side effect
})

test_that("star_index_safe_to_remove() is TRUE when the lock dir doesn't exist yet", {
  index <- tempfile("star_index_") # never created
  expect_true(star_index_safe_to_remove(index, exclude_owner = "anyone"))
})

test_that("star_index_gc() removes an organism's index only when no claim is active, and reports why not otherwise", {
  ref_base <- tempfile("ref_")
  org_dir <- file.path(ref_base, "homo_sapiens")
  genome_dir <- file.path(org_dir, "STAR_index", "genomeDir")
  dir.create(genome_dir, recursive = TRUE)

  config <- fake_config()
  config$config["ref"] <- ref_base

  # "true" is a real system binary that always exits 0 regardless of
  # its arguments -- avoids mocking system2() wholesale, which would
  # also swallow star_index_active_claims()'s own internal
  # system2("kill", ...) PID-liveness calls used by the SAME code path.
  result <- star_index_gc(config, star.path = "true")
  expect_equal(nrow(result), 1)
  expect_true(result$removed)
  expect_identical(result$reason, "ok")

  # Now with an active claim: must NOT attempt removal.
  owner <- star_index_claim(file.path(org_dir, "STAR_index"), owner_id = "busy")
  result2 <- star_index_gc(config, star.path = "true")
  expect_false(result2$removed)
  expect_identical(result2$reason, "other process holds an active claim")
  star_index_release(file.path(org_dir, "STAR_index"), owner)
})
