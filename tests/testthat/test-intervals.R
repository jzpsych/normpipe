iv_inst <- function() read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))

# a graded response model of scale A with fixed parameters
iv_model <- function(s) new_scale_model("grm", s$items, rep(7L, 4), lapply(s$values, sort),
                                        a = c(1.2, 0.9, 1.6, 1.1),
                                        d = replicate(4, seq(2.5, -2.5, length.out = 6), simplify = FALSE))

# one draw of the reference sample (no bootstrap), with the scores of scale A
iv_frame <- function(inst, s, N = 300, seed = 6) {
  sim <- simulate_reference(inst, N = N, seed = seed)
  ref <- dplyr::mutate(sim$meta, score = score_scales(sim$items, inst)$full_a)
  direct_frame(ref, s, inst, n_boot = 0, prior_mass = 1)
}

test_that("the effective N of a norm table is that of the weights, not of one bootstrap draw", {
  inst <- iv_inst(); s <- get_scale(inst, "full_a")
  sim <- simulate_reference(inst, N = 300, seed = 4)
  set.seed(8)
  ref <- dplyr::mutate(sim$meta, weight = runif(dplyr::n(), 0.5, 3),
                       score = score_scales(sim$items, inst)$full_a)
  fr <- direct_frame(ref, s, inst, n_boot = 400, prior_mass = 1, seed = 1)
  nm <- weighted_norms(fr, s, inst, cells = character(0))
  expect_equal(nm$n_eff[1], kish_neff(ref$weight), tolerance = 0.01)
  expect_equal(nm$n[1], 300)
  # a single bootstrap draw has about half of that
  expect_lt(kish_neff(fr$weight[fr$.draw == 1 & fr$shape >= 1]), 0.7 * nm$n_eff[1])
})

test_that("score-conditional intervals take their prior from all draws of the frame", {
  inst <- iv_inst(); s <- get_scale(inst, "full_a"); bl <- make_blocks(s, iv_model(s))
  one <- iv_frame(inst, s)
  two <- dplyr::bind_rows(one, dplyr::mutate(one, .draw = 2L, score = snap_to_grid(.data$score + 1, s)))
  stacked <- dplyr::mutate(two, id = paste(.data$id, .data$.draw), .draw = 1L)   # the same rows as one draw
  nm <- norms_all_cells(one, s, inst)
  a <- conditional_tables(bl, s, two, nm, inst)$overall
  b <- conditional_tables(bl, s, stacked, nm, inst)$overall
  expect_equal(a$T_lo, b$T_lo); expect_equal(a$T_hi, b$T_hi); expect_equal(a$T_eap, b$T_eap)
  first <- conditional_tables(bl, s, one, nm, inst)$overall
  expect_gt(max(abs(a$T_eap - first$T_eap)), 0.5)
})

test_that("integrated intervals mix over the norm table of every draw", {
  inst <- iv_inst(); s <- get_scale(inst, "full_a"); m <- iv_model(s)
  bf <- \(mm) make_blocks(s, mm)
  one <- iv_frame(inst, s)
  stack <- function(shift) purrr::map_dfr(seq_along(shift), \(d)
    dplyr::mutate(one, .draw = d, score = snap_to_grid(.data$score + shift[d], s)))

  # draws that agree: the integrated interval is the score-conditional one
  # (up to the resolution of the latent grid)
  same <- integrated_intervals(bf, list(m, m), s, stack(c(0, 0, 0, 0)), inst)$overall
  cond <- conditional_tables(bf(m), s, one, norms_all_cells(one, s, inst), inst)$overall
  expect_lt(max(abs(same$T_lo - cond$T_lo)), 1)
  expect_lt(max(abs(same$T_hi - cond$T_hi)), 1)

  # draws that disagree strongly: two sets (one per parameter draw) of two
  # draws each, the score distribution shifted down in one draw and up in the
  # other. The definition, written out: the prior pools the draws of a set,
  # and every draw contributes the norm table of its own distribution.
  fr <- stack(c(-1.5, 1.5, -1.5, 1.5))
  wide <- integrated_intervals(bf, list(m, m), s, fr, inst)$overall
  esc <- expected_score_curve(bf(m), s, grid = seq(-6, 6, length.out = 241))
  T_grid <- norm_grid(inst$norm_scale)
  mix <- 0
  for (set in list(1:2, 3:4)) {
    sub <- fr[fr$.draw %in% set, ]
    ep <- eapsum_table(bf(m), s, esc$theta,
                       score_dist_to_theta_prior(sub$score, sub$weight, esc, s$step),
                       return_posterior = TRUE)
    for (d in set) {
      x <- fr[fr$.draw == d, ]
      p <- wecdf_at(x$score, x$weight, s$grid, "midpoint")
      Tm <- make_Tmap(esc, tibble::tibble(raw = s$grid,
                                          T = pctl_to_norm(p, inst$uncertainty$pctl_clamp, inst$norm_scale)))
      mix <- mix + sapply(seq_along(ep$sum_scores), \(j)
        approx(Tm(ep$theta), cumsum(ep$post[, j]), xout = T_grid, yleft = 0, yright = 1, ties = "ordered")$y)
    }
  }
  mix <- mix / 4
  expect_equal(wide$T_lo, T_grid[apply(mix >= inst$uncertainty$probs[1], 2, which.max)])
  expect_equal(wide$T_hi, T_grid[apply(mix >= inst$uncertainty$probs[2], 2, which.max)])

  # pooling the draws of a set into one norm table (normpipe <= 0.4.0) hides
  # the disagreement: both sets then have the same norm table
  pooled <- dplyr::mutate(fr, id = paste(.data$id, .data$.draw), .draw = ifelse(.data$.draw <= 2, 1L, 2L))
  narrow <- integrated_intervals(bf, list(m, m), s, pooled, inst)$overall
  expect_gt(median(wide$T_hi - wide$T_lo), 1.5 * median(narrow$T_hi - narrow$T_lo))
})
