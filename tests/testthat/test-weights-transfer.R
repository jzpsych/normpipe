demo_inst <- function() read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
demo_transfer <- function() read_instrument(system.file("instruments/demo7/demo7_transfer.yml", package = "normpipe"))

rake_data <- function(seed = 1, n = 600) {
  set.seed(seed)
  tibble::tibble(id = seq_len(n), age = sample(18:85, n, TRUE),
                 gender = sample(c("m", "f"), n, TRUE, prob = c(.4, .6)),
                 edu = sample(c("low", "mid", "high"), n, TRUE, prob = c(.2, .3, .5)),
                 weight = rgamma(n, 4, 4))
}
rake_margins <- function(inst) {
  pop <- tidyr::crossing(age = 18:85, gender = c("m", "f")) |> dplyr::mutate(Freq = 1000)
  list(make_margins(pop, inst),
       data.frame(edu = c("low", "mid", "high"), Freq = c(30, 50, 20) * 1e3))
}

test_that("raking hits several margins at once and re-raking restores them after trimming", {
  inst <- demo_inst(); d <- rake_data(); m <- rake_margins(inst)
  r <- rake_to_margins(d, m, inst, trim_quantile = 0.95)
  dg <- attr(r, "diagnostics")
  expect_lt(dg$max_margin_deviation, 2e-3)
  expect_equal(sum(r$weight_raked), nrow(r))
  sh <- tapply(r$weight_raked, r$edu, sum) / sum(r$weight_raked)
  expect_equal(as.vector(sh[c("low", "mid", "high")]), c(.3, .5, .2), tolerance = 2e-3)
  # one-shot trimming (normpipe 0.1.0) leaves the margins off
  r0 <- rake_to_margins(d, m, inst, trim_quantile = 0.95, rerake = FALSE)
  expect_gt(attr(r0, "diagnostics")$max_margin_deviation, dg$max_margin_deviation)
  # the joint cell margin alone still works as before
  r1 <- rake_to_margins(d, m[[1]], inst, rerake = FALSE)
  expect_true(all(c("weight_raked", "raked") %in% names(r1)))
})

test_that("respondents who cannot be raked are dropped or kept with a design-based weight", {
  inst <- demo_inst(); d <- rake_data(2); m <- rake_margins(inst)
  d$age[1:2] <- NA; d$gender[3:5] <- "d"
  rd <- rake_to_margins(d, m, inst)
  expect_equal(nrow(rd), nrow(d) - 5); expect_equal(attr(rd, "diagnostics")$n_not_raked, 5)
  rk <- rake_to_margins(d, m, inst, incomplete = "keep")
  expect_equal(nrow(rk), nrow(d)); expect_equal(sum(!rk$raked), 5)
  expect_equal(sum(rk$weight_raked), nrow(d))
  expect_error(rake_to_margins(dplyr::mutate(d, edu = ifelse(id == 10, "other", edu)), m, inst),
               "without a population count")
})

test_that("transfer specification merges defaults, block and per-outcome overrides", {
  inst <- demo_transfer()
  sp <- transfer_spec(inst, "full_a")
  expect_equal(sp$family, "gaussian"); expect_equal(sp$sigma, "br_mean + s(age)")
  expect_null(sp$exclude)
  inst$transfer$formula <- "{bridge} + s(age, by = gender) + gender + sample"
  inst$transfer$outcomes <- list(full_a = list(family = "skew_normal", adapt_delta = 0.99))
  sp <- transfer_spec(inst, "full_a")
  expect_equal(sp$family, "skew_normal"); expect_equal(sp$adapt_delta, 0.99)
  expect_equal(expand_rhs(sp$formula, c("br_a2", "br_b2")),
               "br_a2 + br_b2 + s(age, by = gender) + gender + sample")
  expect_equal(transfer_spec(inst, "full_b")$family, "gaussian")
  inst$transfer$crosscheck <- list("short_a")
  expect_equal(transfer_spec(inst, "short_a")$exclude, "a2")     # leave-own-domain-out
  inst$transfer$exclude <- list(short_a = c("a2", "b2"))
  expect_equal(transfer_spec(inst, "short_a")$exclude, c("a2", "b2"))
  inst$transfer$family <- "student"
  expect_error(transfer_spec(inst, "full_b"))
})

test_that("projected draws can carry Bayesian-bootstrap weights", {
  inst <- demo_transfer(); s <- get_scale(inst, "full_a")
  pf <- tibble::tibble(id = 1:50, age = 30, gender = "f", weight = rep(c(1, 2), 25))
  pp <- matrix(rnorm(20 * 50, 3, 1), 20, 50)
  f0 <- draws_to_frame(pp, pf, s)
  expect_equal(nrow(f0), 1000); expect_true(all(f0$score %in% s$grid))
  expect_equal(f0$weight, rep(pf$weight, 20))
  set.seed(1); f1 <- draws_to_frame(pp, pf, s, bootstrap_weights = TRUE)
  expect_equal(f1$score, f0$score)
  expect_false(isTRUE(all.equal(f1$weight, f0$weight)))
  expect_equal(mean(f1$weight / f0$weight), 1, tolerance = 0.1)
})

test_that("the cross-check table pairs direct and projected norms by cell and score", {
  inst <- demo_inst(); s <- get_scale(inst, "full_b")
  sim <- simulate_reference(inst, N = 300, seed = 9); sc <- score_scales(sim$items, inst)
  fr <- direct_frame(dplyr::mutate(sim$meta, score = sc$full_b), s, inst, n_boot = 20)
  nm <- norms_all_cells(fr, s, inst)
  cc <- crosscheck_table(nm, nm)
  expect_equal(nrow(cc), nrow(nm$cells) + nrow(nm$overall))
  expect_true(all(cc$dT == 0) && all(cc$inside))
})

test_that("the beta-binomial summed-score representation follows the scale", {
  inst <- demo_transfer()
  a <- bb_scale(get_scale(inst, "full_a"))       # 4 items, values 0..6, item mean
  expect_equal(a$trials, 24L); expect_equal(a$offset, 0L); expect_equal(a$unit_per_sum, 0.25)
  b <- bb_scale(get_scale(inst, "full_b"))       # sum score
  expect_equal(b$trials, 24L); expect_equal(b$unit_per_sum, 1)
  expect_equal(bb_scale(get_scale(inst, "full_ab"))$trials, 48L)   # equal component lengths
  dense <- get_scale(read_instrument(system.file("instruments/pid5/pid5.yml", package = "normpipe")),
                     "full_negative_affectivity")
  expect_error(bb_scale(dense), "count structure")
  sp <- transfer_spec(inst, "full_a")
  expect_equal(sp$prior_b, "normal(0, 1)"); expect_equal(sp$phi, sp$sigma)
  inst$transfer$family <- "beta_binomial"; inst$transfer$phi <- "br_mean"
  sp <- transfer_spec(inst, "full_a")
  expect_equal(sp$prior_b, "normal(0, 1.5)"); expect_equal(sp$phi, "br_mean")
})
