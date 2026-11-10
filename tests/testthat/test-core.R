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

test_that("midpoint percentile ranks sit half a score mass below the cumulative ones", {
  x <- c(5, 5, 5, 6, 7, 7, 9); w <- rep(1, 7); g <- 5:9
  cum <- wecdf_at(x, w, g, "cumulative"); mid <- wecdf_at(x, w, g)
  expect_equal(wecdf_at(x, w, g, "midpoint"), mid)
  expect_equal(cum, c(3, 4, 6, 6, 7) / 7)
  expect_equal(mid, cum - c(3, 1, 2, 0, 1) / 14)
  expect_error(wecdf_at(x, w, g, "other"))
  inst <- demo_inst(); inst$uncertainty$pctl_convention <- "midpoint"
  sim <- simulate_reference(inst, N = 200, seed = 7); sc <- score_scales(sim$items, inst)
  s <- get_scale(inst, "full_b"); ref <- dplyr::mutate(sim$meta, score = sc$full_b)
  fr <- direct_frame(ref, s, inst, n_boot = 0, prior_mass = 0)
  nm_mid <- weighted_norms(fr, s, inst, cells = character(0))   # default = midpoint
  nm_cum <- weighted_norms(fr, s, inst, cells = character(0), convention = "cumulative")
  expect_true(all(nm_mid$pctl_med <= nm_cum$pctl_med + 1e-12))
})

test_that("numeric respondent ids work with the Dirichlet prior", {
  inst <- demo_inst(); s <- get_scale(inst, "full_a")
  sim <- simulate_reference(inst, N = 100, seed = 3); sc <- score_scales(sim$items, inst)
  ref <- dplyr::mutate(sim$meta, id = seq_len(dplyr::n()), score = sc$full_a)
  fr <- direct_frame(ref, s, inst, n_boot = 2, prior_mass = 1)
  expect_type(fr$id, "character")
  expect_true(all(c("1", "100") %in% fr$id))
})

test_that("skipped items are imputed within the response categories and reproducibly", {
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 300, seed = 4)
  expect_identical(impute_items(sim$items, inst, M = 2)[[1]], sim$items)   # complete data unchanged
  it <- sim$items
  set.seed(1); holes <- cbind(sample(300, 60), sample(8, 60, TRUE)); it[holes] <- NA
  it[7, ] <- NA                                                           # nothing answered
  imp <- impute_items(it, inst, M = 5, seed = 2)
  expect_length(imp, 5)
  expect_equal(sum(attr(imp, "n_imputed")), sum(is.na(it)))
  for (d in imp) {
    expect_false(anyNA(d))
    expect_true(all(as.matrix(d) %in% 1:7))
    expect_identical(d[!is.na(it)], it[!is.na(it)])                       # answered items stay
  }
  expect_identical(impute_items(it, inst, M = 5, seed = 2), imp)
  expect_false(identical(imp[[1]], imp[[2]]))
  # the draw follows the answered items of the block (here item 2)
  it2 <- sim$items; hi <- rowMeans(it2[, 2:4]) >= 5.5; lo <- rowMeans(it2[, 2:4]) <= 2.5
  it2$it1[which(hi)[1:5]] <- NA; it2$it1[which(lo)[1:5]] <- NA
  d <- do.call(cbind, lapply(impute_items(it2, inst, M = 40), \(x) x$it1))
  expect_gt(mean(d[which(hi)[1:5], ]), mean(d[which(lo)[1:5], ]) + 1)
  # the bridge pair is the block when the instrument has bridge predictors
  tr <- read_instrument(system.file("instruments/demo7/demo7_transfer.yml", package = "normpipe"))
  expect_equal(imputation_blocks(tr, 1:8)[["1"]], c(1L, 2L))
  expect_equal(imputation_blocks(inst, 1:8)[["1"]], 1:2)       # shortest scale: short_a
  expect_equal(imputation_blocks(inst, 1:8)[["3"]], 1:4)
  expect_equal(imputation_blocks(tr, c(1, 2, 5, 6))[["5"]], c(5L, 6L))
})

test_that("direct frames take the scores of the completed data sets draw by draw", {
  inst <- demo_inst(); inst$versions$full$scales$b$min_answered <- 0.5
  s <- get_scale(inst, "full_b")
  sim <- simulate_reference(inst, N = 120, seed = 5)
  it <- sim$items; it$it5[1:10] <- NA
  sc <- score_scales(it, inst)
  imp <- impute_items(it, inst, M = 3)
  si <- sapply(imp, \(d) score_scales(d, inst)$full_b)
  ref <- dplyr::mutate(sim$meta, score = sc$full_b)
  ref$score[11] <- NA                                    # excluded respondent stays out
  f0 <- direct_frame(ref, s, inst, n_boot = 7)
  f1 <- direct_frame(ref, s, inst, n_boot = 7, scores_imp = si)
  expect_equal(f1$weight, f0$weight); expect_equal(f1$id, f0$id)
  id <- as.character(sim$meta$id)
  expect_false(id[11] %in% f1$id)
  for (r in 1:7) {
    x <- f1[f1$.draw == r, ]
    expect_equal(x$score[match(id[-11], x$id)], si[-11, (r - 1) %% 3 + 1])
  }
  complete <- f1$id %in% id[12:120]
  expect_equal(f1$score[complete], f0$score[complete])
  expect_equal(f1$score[!f1$id %in% id], f0$score[!f0$id %in% id])   # prior rows untouched
  expect_error(direct_frame(ref, s, inst, scores_imp = si[-1, ]))
})

test_that("run_norming() norms a reference sample with skipped items from completed data", {
  inst <- demo_inst(); inst$uncertainty$n_boot <- 10
  sim <- simulate_reference(inst, N = 300, seed = 6)
  set.seed(2); sim$items[cbind(sample(300, 40), sample(8, 40, TRUE))] <- NA
  inst$missing$reference <- "impute"; inst$missing$n_imputations <- 4
  for (v in names(inst$versions)) for (k in names(inst$versions[[v]]$scales))
    inst$versions[[v]]$scales[[k]]$min_answered <- 0.5
  out <- file.path(tempdir(), "np_impute"); unlink(out, recursive = TRUE)
  suppressMessages(run_norming(inst, sim, out, stages = c("P0", "M3a"), verbose = FALSE))
  expect_true(file.exists(file.path(out, "p0_imputed_items.rds")))
  fr <- frame_long(readRDS(file.path(out, "frame_full_b.rds")))
  expect_false(anyNA(fr$score))
  expect_true(all(fr$score %in% get_scale(inst, "full_b")$grid))
  expect_gt(length(unique(table(fr$score[fr$.draw == 1]) - table(factor(fr$score[fr$.draw == 2], levels = sort(unique(fr$score[fr$.draw == 1])))))), 1)
  inst$missing$reference <- "listwise"
  expect_error(suppressMessages(run_norming(inst, sim, file.path(tempdir(), "np_bad"), stages = "P0", verbose = FALSE)), "prorate")
})
