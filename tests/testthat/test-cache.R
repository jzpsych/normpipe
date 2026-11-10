test_that("a cached result is reused for the same key and recomputed for a new one", {
  out <- file.path(tempdir(), "np_cache_test"); unlink(out, recursive = TRUE)
  n <- 0
  f <- function() { n <<- n + 1; n }
  k1 <- cache_key(list(a = 1), "x"); k2 <- cache_key(list(a = 2), "x")
  expect_false(identical(k1, k2)); expect_identical(k1, cache_key(list(a = 1), "x"))
  expect_equal(cache_get(out, "res", k1, f), 1)
  expect_equal(cache_get(out, "res", k1, f), 1)          # reused
  expect_true(cache_valid(out, "res", k1)); expect_false(cache_valid(out, "res", k2))
  msg <- character(0)
  expect_equal(cache_get(out, "res", k2, f, say = \(...) msg <<- paste0(...)), 2)
  expect_match(msg, "inputs changed")
  expect_length(list.files(file.path(out, "superseded"), pattern = "^res_"), 1)   # old result kept
  expect_equal(readRDS(list.files(file.path(out, "superseded"), full.names = TRUE)), 1)
  st <- cache_status(out)
  expect_equal(st$name, "res"); expect_equal(st$key, k2); expect_true(st$exists)
  # a file without a manifest entry (written by an older version) is not trusted
  saveRDS(99, file.path(out, "legacy.rds"))
  expect_equal(cache_get(out, "legacy", k1, function() 5), 5)
  expect_equal(nrow(cache_status(out)), 2)
})

test_that("run_norming() recomputes only what depends on a change", {
  inst <- read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
  inst$uncertainty$n_boot <- 5
  sim <- simulate_reference(inst, N = 200, seed = 11)
  out <- file.path(tempdir(), "np_cache_run"); unlink(out, recursive = TRUE)
  st <- c("P0", "M3a", "M3b")
  run_norming(inst, sim, out, stages = st, outcomes = c("full_a", "full_b"), verbose = FALSE)
  m1 <- cache_status(out)
  expect_setequal(m1$name, c("p0_reference", paste0(rep(c("frame_", "m3a_norms_", "m3b_continuous_"), 2),
                                                     rep(c("full_a", "full_b"), each = 3))))
  run_norming(inst, sim, out, stages = st, outcomes = c("full_a", "full_b"), verbose = FALSE)
  expect_equal(cache_status(out)$time, m1$time)                 # nothing recomputed
  expect_false(dir.exists(file.path(out, "superseded")))
  # more bootstrap draws: frames and norms are outdated, the reference scores are not
  inst$uncertainty$n_boot <- 6
  run_norming(inst, sim, out, stages = st, outcomes = c("full_a", "full_b"), verbose = FALSE)
  m2 <- cache_status(out)
  ch <- m2$key != m1$key[match(m2$name, m1$name)]
  expect_equal(sort(m2$name[!ch]), "p0_reference")
  expect_length(list.files(file.path(out, "superseded")), 6)
  # a changed answer of one respondent: everything that uses the reference sample
  sim$items$it5[1] <- if (sim$items$it5[1] == 1) 2 else 1
  run_norming(inst, sim, out, stages = "M3a", outcomes = "full_b", verbose = FALSE)
  m3 <- cache_status(out)
  expect_true(all(m3$key[m3$name %in% c("p0_reference", "frame_full_b", "m3a_norms_full_b")] !=
                    m2$key[m2$name %in% c("p0_reference", "frame_full_b", "m3a_norms_full_b")]))
  # stages one by one on a context
  ctx <- norming_context(inst, sim, out, outcomes = "full_b", verbose = FALSE)
  run_stage(ctx, "M3a")
  expect_setequal(ctx$done, c("P0", "M3a"))
  expect_error(run_stage(ctx, "M9"), "Unknown stage")
})
