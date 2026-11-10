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
  expect_equal(sp$prior_b_phi, "normal(0, 2.5)")
})

test_that("observed bridge items are not predicted but added (all families)", {
  inst <- demo_transfer()
  s <- get_scale(inst, "full_a")                      # items 1-4 (values 0..6, mean), bridge a2 = items 1, 2
  expect_equal(observed_items(s, inst), c(1L, 2L))
  expect_length(observed_items(s, inst, exclude = "a2"), 0)
  # weights of the item values in the raw score
  expect_equal(unname(item_weights(s, inst)), rep(1 / 4, 4))
  expect_equal(unname(item_weights(get_scale(inst, "full_b"), inst)), rep(1, 4))
  expect_equal(unname(item_weights(get_scale(inst, "full_ab"), inst)), rep(1 / 8, 8))
  set.seed(3)
  it <- as.data.frame(matrix(sample(1:7, 400, TRUE), 50, 8, dimnames = list(NULL, paste0("it", 1:8))))
  bs <- bridge_scores(it, inst)
  expect_true(all(c("bi_1", "bi_2", "bi_5", "bi_6") %in% names(bs)))
  info <- observed_info(s, inst, c(1L, 2L), bs)
  expect_equal(info$range_units, c(min = 0, max = 12))
  expect_equal(info$rest, c(min = 0, max = 3))           # items 3, 4: (0..6 + 0..6) / 4
  bb <- bb_scale(s, known = c(1L, 2L), values = info$range_units)
  expect_equal(bb$trials, 12L); expect_equal(bb$offset, 0L)
  expect_error(bb_scale(get_scale(inst, "short_a"), known = c(1L, 2L), values = info$range_units),
               "all items are bridge items")
  # observed part: raw contribution of items 1 and 2 to the scale mean
  sc <- score_scales(it, inst)
  op <- observed_part(bs, info)
  expect_equal(op, (bs$bi_1 + bs$bi_2) / 4)
  expect_true(all(sc$full_a - op >= info$rest[["min"]] - 1e-9 & sc$full_a - op <= info$rest[["max"]] + 1e-9))
  # a skipped observed item is drawn anew in every draw from its distribution
  # given the other item of the pair; answered items stay fixed
  bs2 <- bs; bs2$bi_1[1] <- NA
  m <- observed_part(bs2, info, D = 200)
  expect_equal(dim(m), c(200, 50))
  expect_true(all(m[, -1] == matrix(op[-1], 200, 49, byrow = TRUE)))
  expect_gt(length(unique(m[, 1])), 1)
  expect_true(all(round((m[, 1] - bs2$bi_2[1] / 4) * 4, 8) %in% 0:6))
})

test_that("extra priors replace defaults of the same class and dpar; init is checked", {
  inst <- demo_transfer()
  inst$transfer$priors <- list(list(prior = "normal(0, 0.5)", class = "b", dpar = "sigma"),
                               list(prior = "normal(0, 2)", class = "b"))
  pr <- transfer_priors(transfer_spec(inst, "full_a"), "gaussian")
  expect_equal(nrow(pr), 2)
  expect_equal(pr$prior[pr$dpar == ""], "normal(0, 2)")
  inst$transfer$init <- "x"
  expect_error(transfer_spec(inst, "full_a"), "random")
  inst$transfer$init <- 0
  expect_equal(transfer_spec(inst, "full_a")$init, 0)
})

test_that("full_init gives every parameter of a censored model a value", {
  skip_if_not_installed("brms")
  set.seed(1)
  d <- data.frame(x = rnorm(50), age = runif(50, 20, 80))
  d$y <- round(pmin(pmax(1 + 0.5 * d$x + rnorm(50, 0, 0.5), 0), 3), 1)
  d$lo <- d$y - 0.05; d$hi <- d$y + 0.05; d$ci <- "interval"
  f <- brms::bf(lo | cens(ci, hi) ~ x + s(age), sigma ~ x)
  ini <- full_init(f, d, stats::gaussian(), list(Intercept = 1, Intercept_sigma = log(0.5)))()
  expect_setequal(names(ini), c("b", "Intercept", "bs", "zs_1_1", "sds_1", "b_sigma", "Intercept_sigma"))
  expect_equal(ini$Intercept, 1, tolerance = 0.3)
  expect_equal(ini$Intercept_sigma, log(0.5), tolerance = 0.3)
  expect_length(ini$b, 1)
  expect_true(all(ini$sds_1 > 0))
})

test_that("several candidate families are accepted; priors follow the family of each fit", {
  inst <- demo_transfer()
  inst$transfer$family <- list("beta_binomial", "gaussian")
  sp <- transfer_spec(inst, "full_a")
  expect_equal(sp$families, c("beta_binomial", "gaussian")); expect_equal(sp$family, "beta_binomial")
  expect_equal(sp$select_rhat, 1.01)
  skip_if_not_installed("brms")
  pb <- function(fam) { p <- transfer_priors(sp, fam); p$prior[p$class == "b" & p$dpar == ""] }
  expect_equal(pb("beta_binomial"), "normal(0, 1.5)"); expect_equal(pb("gaussian"), "normal(0, 1)")
  sp$prior_b_user <- "normal(0, 2)"
  expect_equal(pb("gaussian"), "normal(0, 2)")
  inst$transfer$family <- list("beta_binomial", "poisson")
  expect_error(transfer_spec(inst, "full_a"), "poisson")
})
