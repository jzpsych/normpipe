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
