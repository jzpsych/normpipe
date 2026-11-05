# ==============================================================================
# intervals.R -- individual intervals on the T metric
#
# Score-conditional intervals (Supplement S5.2): per reporting cell the
# weighted score distribution is carried onto the latent metric as the
# population prior, the summed-score EAP gives the posterior of the trait per
# attainable score, and the expected-score curve together with the norm table
# of the cell maps it onto T.
#
# Integrated intervals (Supplement S5.4): the same computation repeated under
# joint draws of item parameters (draw_model_params) and norm-level draws
# (the draw frame: Bayesian bootstrap of the direct route, or transfer
# posterior draws), with the per-draw T posteriors mixed on a common T grid.
# ==============================================================================

#' Population prior on the latent grid from a weighted score distribution
#'
#' Every attainable score s stands for a latent interval (s - h/2, s + h/2]
#' on the expected-score curve; its weighted mass is spread uniformly over the
#' theta interval that the curve assigns to it (floor and ceiling intervals
#' closed at the ends of the grid).
score_dist_to_theta_prior <- function(scores, w, esc, resolution) {
  th <- esc$theta; n <- length(th)
  edges <- c(th[1] - (th[2] - th[1]) / 2, (th[-1] + th[-n]) / 2, th[n] + (th[n] - th[n - 1]) / 2)
  h <- resolution
  keep <- !is.na(scores) & !is.na(w)
  s <- round(scores[keep] / h) * h
  mass <- tapply(w[keep], round(s, 8), sum)
  vals <- as.numeric(names(mass))
  wt <- numeric(n)
  for (i in seq_along(vals)) {
    lo_s <- vals[i] - h / 2; hi_s <- vals[i] + h / 2
    th_lo <- if (lo_s <= min(esc$mean_score)) edges[1] else approx(esc$mean_score, esc$theta, xout = lo_s, rule = 2)$y
    th_hi <- if (hi_s >= max(esc$mean_score)) edges[n + 1] else approx(esc$mean_score, esc$theta, xout = hi_s, rule = 2)$y
    overlap <- pmax(0, pmin(edges[-1], th_hi) - pmax(edges[-(n + 1)], th_lo))
    if (sum(overlap) <= 0) { k <- which.min(abs(th - (th_lo + th_hi) / 2)); overlap[k] <- 1 }
    wt <- wt + mass[[i]] * overlap / sum(overlap)
  }
  wt / sum(wt)
}

#' theta -> T through the expected-score curve and a norm table (raw, T)
make_Tmap <- function(esc, norm_tab) {
  function(theta) {
    m <- approx(esc$theta, esc$mean_score, xout = theta, rule = 2)$y
    approx(norm_tab$raw, norm_tab$T, xout = m, rule = 2)$y
  }
}

#' Split a frame into the reporting cells (list: cell key -> subset), with the
#' overall table last
split_cells <- function(fr, inst) {
  fr <- add_age_group(fr, inst)
  cv <- cell_vars(inst)
  out <- list()
  if (length(cv)) {
    sub <- fr
    if ("gender" %in% cv) sub <- filter(sub, .data$gender %in% inst$cells$gender)
    if ("age_group" %in% cv) sub <- filter(sub, !is.na(.data$age_group))
    cells <- distinct(sub, across(all_of(cv))) %>% arrange(across(all_of(cv)))
    for (r in seq_len(nrow(cells))) {
      key <- cell_key(cells[r, , drop = FALSE])
      sel <- rep(TRUE, nrow(sub))
      for (v in cv) sel <- sel & as.character(sub[[v]]) == as.character(cells[[v]][r])
      out[[key]] <- sub[sel, ]
    }
  }
  if (isTRUE(inst$cells$overall) || !length(cv)) out$overall <- fr
  out
}

cell_key <- function(row) paste(map_chr(row, as.character), collapse = " ")

#' Norm table of one cell from a norms_all_cells() result
norm_table_of_cell <- function(nm, key, inst) {
  if (key == "overall") return(nm$overall)
  cv <- cell_vars(inst)
  keys <- apply(as.data.frame(nm$cells[cv]), 1, \(r) paste(r, collapse = " "))
  nm$cells[keys == key, ]
}

#' Score-conditional EAP tables of a scale for every reporting cell
#'
#' @param blocks block list of the scale (point estimates)
#' @param fr point frame (or any frame; the point weights are used)
#' @param nm norm tables of the scale (norms_all_cells)
#' @return named list (cell key) of tibbles with sum_score, raw, T_eap, T_lo,
#'   T_hi, tau_eap, tau_sd, p_marg
conditional_tables <- function(blocks, s, fr, nm, inst, sum_scores = NULL,
                               probs = inst$uncertainty$probs, ...) {
  fr <- point_frame(fr)
  esc <- expected_score_curve(blocks, s)
  gh <- statmod::gauss.quad.prob(21, "normal")
  prep <- map(blocks, block_prep, grid = esc$theta, gh = gh)
  sum_scores <- sum_scores %||% eap_sum_rows(s)
  cells <- split_cells(fr, inst)
  imap(cells, \(sub, key) {
    ntab <- norm_table_of_cell(nm, key, inst)
    prw <- score_dist_to_theta_prior(sub$score, sub$weight, esc, s$step)
    tab <- eapsum_table(blocks, s, esc$theta, prw, sum_scores = sum_scores,
                        probs = probs, prep = prep, ...)
    Tmap <- make_Tmap(esc, ntab)
    tab %>% mutate(T_eap = Tmap(.data$theta_eap), T_lo = Tmap(.data$theta_lo),
                   T_hi = Tmap(.data$theta_hi))
  })
}

#' Integrated intervals of a scale for every reporting cell
#'
#' @param block_fun function(models) -> block list, where models is one
#'   parameter draw of the scale's model(s) (see draws_of_blocks)
#' @param model_draws list of n_param_draws model (or model list) draws
#' @param fr draw frame with several draws (Bayesian bootstrap or projection)
integrated_intervals <- function(blocks_fun, model_draws, s, fr, inst,
                                 sum_scores = NULL, probs = inst$uncertainty$probs,
                                 T_grid = seq(10, 90, by = 0.1),
                                 clamp = inst$uncertainty$pctl_clamp,
                                 convention = inst$uncertainty$pctl_convention %||% "midpoint", ...) {
  n_param <- length(model_draws)
  sum_scores <- sum_scores %||% eap_sum_rows(s)
  draw_ids <- sort(unique(fr$.draw))
  chunks <- split(draw_ids, cut(seq_along(draw_ids), n_param, labels = FALSE))
  gh <- statmod::gauss.quad.prob(21, "normal")
  theta <- seq(-6, 6, length.out = 241)
  cells <- split_cells(fr, inst)
  acc <- map(cells, \(x) matrix(0, length(T_grid), length(sum_scores)))
  for (r in seq_len(n_param)) {
    bl <- blocks_fun(model_draws[[r]])
    prep <- map(bl, block_prep, grid = theta, gh = gh)
    esc <- expected_score_curve(bl, s, grid = theta, prep = prep)
    for (key in names(cells)) {
      sub <- cells[[key]] %>% filter(.data$.draw %in% chunks[[r]])
      pw <- score_dist_to_theta_prior(sub$score, sub$weight, esc, s$step)
      p <- wecdf_at(sub$score, sub$weight, s$grid, convention)
      Tmap <- make_Tmap(esc, tibble(raw = s$grid, T = pctl_to_T(p, clamp)))
      ep <- eapsum_table(bl, s, esc$theta, pw, sum_scores = sum_scores,
                         return_posterior = TRUE, prep = prep, ...)
      T_theta <- Tmap(ep$theta)
      for (j in seq_along(sum_scores))
        acc[[key]][, j] <- acc[[key]][, j] +
          approx(x = T_theta, y = cumsum(ep$post[, j]), xout = T_grid,
                 yleft = 0, yright = 1, ties = "ordered")$y
    }
  }
  map(acc, \(A) map_dfr(seq_along(sum_scores), \(j) {
    mix <- A[, j] / n_param
    tibble(sum_score = sum_scores[j], raw = sum_to_raw(sum_scores[j], s),
           T_eap = sum(head(T_grid, -1) * diff(mix)) + T_grid[1] * mix[1],
           T_lo = T_grid[which.max(mix >= probs[1])],
           T_hi = T_grid[which.max(mix >= probs[2])])
  }))
}

#' Marginal reliability of the summed score in the reference population
#' (from the overall conditional table; Supplement S7.1)
marginal_reliability <- function(eap_overall) {
  d <- eap_overall %>% filter(!is.na(.data$tau_eap), !is.na(.data$tau_sd), .data$p_marg > 0)
  w <- d$p_marg / sum(d$p_marg)
  m <- sum(w * d$tau_eap)
  vb <- sum(w * (d$tau_eap - m)^2); vw <- sum(w * d$tau_sd^2)
  vb / (vb + vw)
}
