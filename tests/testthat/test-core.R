demo_inst <- function() read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))

test_that("item table, values and grids follow the configuration", {
  inst <- demo_inst()
  it <- inst$items
  expect_equal(nrow(it), 8)
  expect_equal(it$values[[1]], 0:6)
  expect_equal(it$values[[3]], 6:0)          # reverse-keyed
  a <- get_scale(inst, "full_a"); b <- get_scale(inst, "full_b")
  expect_equal(a$grid, seq(0, 6, by = 0.25)); expect_equal(a$step, 0.25)
  expect_equal(b$grid, 0:24); expect_equal(b$step, 1)
  ab <- get_scale(inst, "full_ab")
  expect_false(ab$dense)                     # equal component lengths: finite grid
  expect_equal(ab$grid, seq(0, 6, by = 0.125))  # = item mean over 8 items
  expect_equal(attainable_sums(list(c(0, 2), c(0, 1))), c(0, 1, 2, 3))
  expect_equal(attainable_sums(list(c(0, 2), c(0, 2))), c(0, 2, 4))
})

test_that("snap_to_grid clamps and rounds ties up", {
  s <- list(grid = seq(0, 6, by = 0.25))
  expect_equal(snap_to_grid(c(-1, 0.1, 0.125, 0.2, 7, NA), s), c(0, 0, 0.25, 0.25, 6, NA))
})

test_that("Lord-Wingersky recursion equals brute-force pattern summation", {
  set.seed(3)
  K <- c(3L, 4L, 2L); vals <- list(0:2, 3:0, c(0L, 1L))   # one reversed item
  th <- seq(-2, 2, length.out = 5); Theta <- matrix(th)
  P <- lapply(seq_along(K), \(i) grm_probs(runif(1, 0.8, 2), sort(rnorm(K[i] - 1), decreasing = TRUE), Theta))
  L <- lw_recursion(P, vals)
  pats <- expand.grid(lapply(K, seq_len))
  brute <- matrix(0, length(th), ncol(L))
  for (r in seq_len(nrow(pats))) {
    pr <- Reduce(`*`, lapply(seq_along(K), \(i) P[[i]][, pats[r, i]]))
    S <- sum(mapply(\(v, k) v[k], vals, pats[r, ]))
    brute[, S + 1] <- brute[, S + 1] + pr
  }
  expect_equal(L, brute, tolerance = 1e-12)
  expect_equal(rowSums(L), rep(1, length(th)), tolerance = 1e-12)
})

test_that("summed-score EAP is monotone for a unidimensional GRM", {
  inst <- demo_inst(); s <- get_scale(inst, "full_a")
  # model values are increasing (value-ordered categories, as fit_scale_model builds them)
  m <- new_scale_model("grm", s$items, rep(7L, 4), lapply(s$values, sort),
                       a = c(1.2, 0.9, 1.6, 1.1),
                       d = replicate(4, seq(2.5, -2.5, length.out = 6), simplify = FALSE))
  bl <- make_blocks(s, m)
  th <- seq(-6, 6, length.out = 241); pw <- dnorm(th); pw <- pw / sum(pw)
  tab <- eapsum_table(bl, s, th, pw)
  expect_equal(tab$sum_score, 0:24)
  expect_true(all(diff(tab$theta_eap) > 0))
  expect_true(all(diff(tab$tau_eap) > 0))
  expect_equal(sum(tab$p_marg), 1, tolerance = 1e-8)
  esc <- expected_score_curve(bl, s)
  expect_true(all(diff(esc$mean_score) > 0))
  expect_equal(range(esc$mean_score), c(0, 6), tolerance = 0.05)
})

test_that("a frequency table gives the same point norms as the individual data", {
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 400, seed = 11)
  sc <- score_scales(sim$items, inst)
  s <- get_scale(inst, "full_b")
  ref <- dplyr::mutate(sim$meta, score = sc$full_b)
  fr1 <- direct_frame(ref, s, inst, n_boot = 0, prior_mass = 0)
  freq <- add_age_group(ref, inst) |>
    dplyr::count(gender, age_group, score, wt = weight, name = "count")
  fr2 <- frequency_frame(freq, s, inst, n_boot = 0, prior_mass = 0)
  n1 <- norms_all_cells(fr1, s, inst); n2 <- norms_all_cells(fr2, s, inst)
  expect_equal(n1$overall$T, n2$overall$T, tolerance = 1e-8)
  expect_equal(n1$cells$T, n2$cells$T, tolerance = 1e-8)
})

test_that("the Dirichlet prior adds exactly prior_mass per cell and keeps T monotone", {
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 300, seed = 5)
  sc <- score_scales(sim$items, inst); s <- get_scale(inst, "full_a")
  ref <- dplyr::mutate(sim$meta, score = sc$full_a)
  fr <- direct_frame(ref, s, inst, n_boot = 0, prior_mass = 1)
  pri <- fr[grepl("^prior", fr$id), ]
  expect_true(nrow(pri) > 0)
  tot <- tapply(pri$shape, paste(pri$gender, pri$age_group), sum)
  expect_equal(as.vector(tot), rep(1, length(tot)), tolerance = 1e-8)
  nm <- norms_all_cells(fr, s, inst)
  expect_true(all(tapply(nm$cells$T, paste(nm$cells$gender, nm$cells$age_group), \(x) all(diff(x) >= 0))))
  expect_true(all(nm$cells$pctl_med > 0))
  # PR = P(X <= x): only the top score reaches 1 (clamped on the T scale)
  expect_true(all(nm$cells$pctl_med[nm$cells$raw < max(s$grid)] < 1))
})

test_that("partial item sets score only the scales whose items are present", {
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 50, seed = 2)
  full <- score_scales(sim$items, inst)
  part <- score_scales(sim$items[, c("it1", "it2", "it5", "it6")], inst, partial = TRUE)
  expect_equal(part$short_a, full$short_a)
  expect_true(all(is.na(part$full_a)) && all(is.na(part$full_total)))
  expect_error(score_scales(sim$items[, c("it1", "it2")], inst), "missing")
})
