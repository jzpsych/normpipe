# ==============================================================================
# transfer-project.R -- projection of the reference sample through a transfer
# model, and the leave-own-domain-out cross-check
# ==============================================================================

#' Project the reference sample through a transfer model
#'
#' Posterior predictive draws of the unobserved part of the outcome for every
#' reference respondent, plus the respondent's observed part, snapped to the
#' attainable scores. The result is the draw frame of the outcome.
#'
#' @param fit Transfer model ([transfer_model()]).
#' @param pred_frame Prediction frame of the reference sample: `id`, `age`,
#'   `gender`, `weight`, the bridge columns and `sample`; or a list of such
#'   frames, one per completed data set of the reference sample
#'   ([impute_items()]), over which the posterior draws are spread evenly.
#' @param s Scale specification ([get_scale()]).
#' @param n_draws Number of posterior predictive draws per respondent.
#' @param level Level of `sample` at which the reference sample is projected:
#'   a donor sample, or `"average"` (the mean of the sample effects over all
#'   donor samples). Default: the level already set in `pred_frame`.
#' @param bootstrap_weights If `TRUE`, the reference weights are replaced in
#'   every posterior draw by `weight x Exp(1)` (Bayesian bootstrap), so that
#'   the norm band contains the sampling uncertainty of the reference sample
#'   as well as the uncertainty of the transfer model.
#' @return Draw frame (see [direct_frame()]).
#' @export
project_outcome <- function(fit, pred_frame, s, n_draws = 500, level = NULL,
                            bootstrap_weights = FALSE) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  lv_all <- levels(fit$data$sample) %||% unique(as.character(fit$data$sample))
  fam <- stats::family(fit)$family
  is_count <- fam %in% count_families
  info <- fit$normpipe
  if (is.null(info))
    stop("This transfer fit has no record of the observed items (made before normpipe 0.2.4); refit it.")
  known <- info$known_items
  bb <- if (is_count) bb_scale(s, known = known, values = info$range_units)
  if (length(lv_all) < 2) level <- NULL            # no sample effect in the model
  if (!is.null(level) && !identical(level, "average") && !level %in% lv_all)
    stop("Projection level '", level, "' is not a donor sample (",
         paste(lv_all, collapse = ", "), ").")
  frames <- if (inherits(pred_frame, "data.frame")) list(pred_frame) else pred_frame
  ids <- sample(seq_len(brms::ndraws(fit)), min(n_draws, brms::ndraws(fit)))
  M <- min(length(frames), length(ids))
  chunk <- split(ids, rep(seq_len(M), length.out = length(ids)))
  pp <- do.call(rbind, map(seq_len(M), \(m) {
    pf <- frames[[m]]
    if (fam == "beta_binomial") pf$y_trials <- bb$trials
    if (identical(level, "average")) {
      x <- predict_average_level(fit, pf, levels = lv_all, ids = chunk[[m]])
    } else {
      if (!is.null(level)) pf$sample <- factor(level, levels = lv_all)
      x <- brms::posterior_predict(fit, newdata = pf, draw_ids = chunk[[m]])  # draws x persons
    }
    if (is_count) {
      if (fam == "cumulative") x <- x - 1L                      # category 1..K+1 -> sum 0..K
      x <- x + bb$offset
      if (length(known)) x <- x + observed_part(pf, info, wt = rep(1, length(known)), D = nrow(x))
      x * bb$unit_per_sum                                       # summed units -> raw metric
    } else if (length(known)) {
      # latent rest, censored at its bounds, plus the observed part
      pmin(pmax(x, info$rest[["min"]]), info$rest[["max"]]) + observed_part(pf, info, D = nrow(x))
    } else x
  }))
  draws_to_frame(pp, frames[[1]], s, bootstrap_weights)
}

#' Draw frame from a draws x persons matrix of latent predictions
#'
#' Snaps the predictions to the attainable scores and, optionally, draws a
#' Bayesian-bootstrap weight per person and draw.
#' @noRd
draws_to_frame <- function(pp, pred_frame, s, bootstrap_weights = FALSE) {
  n <- nrow(pred_frame); D <- nrow(pp)
  stopifnot(ncol(pp) == n)
  w <- rep(pred_frame$weight, times = D)
  if (isTRUE(bootstrap_weights)) w <- w * stats::rexp(n * D)
  tibble(id = rep(pred_frame$id, times = D),
         age = rep(pred_frame$age, times = D),
         gender = rep(as.character(pred_frame$gender), times = D),
         weight = w,
         shape = 1,
         .draw = rep(seq_len(D), each = n),
         score = snap_to_grid(as.vector(t(pp)), s))
}

#' Posterior predictions at the average of the sample effects
#'
#' The linear predictor of the mean is averaged over all donor levels of
#' `sample` (equal weights) in every posterior draw; the other distributional
#' parameters (and the thresholds of the cumulative model) are taken from the
#' same draw. Returns draws on the scale of `posterior_predict()`.
#' @noRd
predict_average_level <- function(fit, pred_frame, n_draws = NULL, levels, ids = NULL) {
  ids <- ids %||% sample(seq_len(brms::ndraws(fit)), min(n_draws, brms::ndraws(fit)))
  mu <- Reduce(`+`, map(levels, \(lv) {
    nd <- pred_frame; nd$sample <- factor(lv, levels = levels)
    brms::posterior_linpred(fit, newdata = nd, dpar = "mu", draw_ids = ids)
  })) / length(levels)
  nd <- pred_frame; nd$sample <- factor(levels[1], levels = levels)
  fam <- stats::family(fit)$family
  if (fam == "cumulative") {
    # category (1..K+1, as posterior_predict) by inversion of
    # P(y <= k) = logistic(disc * (threshold_k - mu))
    disc <- brms::posterior_linpred(fit, newdata = nd, dpar = "disc", transform = TRUE, draw_ids = ids)
    thr <- unclass(brms::as_draws_matrix(fit, variable = "^b_Intercept\\[", regex = TRUE))[ids, , drop = FALSE]
    u <- matrix(stats::runif(length(mu)), nrow(mu), ncol(mu))
    out <- matrix(1L, nrow(mu), ncol(mu))
    for (k in seq_len(ncol(thr))) out <- out + (u > stats::plogis(disc * (thr[, k] - mu)))
    return(out)
  }
  if (fam == "beta_binomial") {
    mu <- stats::plogis(mu)
    phi <- brms::posterior_linpred(fit, newdata = nd, dpar = "phi", transform = TRUE, draw_ids = ids)
    if (!requireNamespace("extraDistr", quietly = TRUE)) stop("extraDistr is not installed.")
    out <- extraDistr::rbbinom(length(mu), size = pred_frame$y_trials[1], alpha = mu * phi, beta = (1 - mu) * phi)
    return(matrix(out, nrow(mu), ncol(mu)))
  }
  sig <- brms::posterior_linpred(fit, newdata = nd, dpar = "sigma", transform = TRUE, draw_ids = ids)
  if (fam == "gaussian") {
    out <- stats::rnorm(length(mu), mu, sig)
  } else if (fam == "skew_normal") {
    alpha <- brms::as_draws_matrix(fit, variable = "alpha")[ids, 1]
    alpha <- matrix(alpha, nrow(mu), ncol(mu))
    out <- brms::rskew_normal(length(mu), mu = mu, sigma = sig, alpha = alpha)
  } else stop("Family '", fam, "' is not supported for the average projection level.")
  matrix(out, nrow(mu), ncol(mu))
}

#' Prediction frames of the reference sample for a fitted transfer model
#'
#' One frame per completed item data set of the reference sample, built from
#' the reference scores and items alone; the factor levels of `gender` and
#' `sample` come from the fit, `sample` is set to its reference (first)
#' level. Respondents without age, or with a gender the fit does not know,
#' are left out (they cannot be projected).
#'
#' @param ref_scores Scores of the reference sample (`id`, `age`, `gender`,
#'   `weight`).
#' @param ref_ids Ids of the rows of the item data sets.
#' @param items_list List of item data sets (one per completed data set, or a
#'   single data frame).
#' @param fit Transfer model.
#' @param inst Instrument.
#' @return List of prediction frames.
#' @noRd
prediction_frames <- function(ref_scores, ref_ids, items_list, fit, inst) {
  if (inherits(items_list, "data.frame")) items_list <- list(items_list)
  g_lv <- levels(fit$data$gender) %||% unique(as.character(fit$data$gender))
  s_lv <- levels(fit$data$sample) %||% unique(as.character(fit$data$sample))
  ok <- !is.na(ref_scores$age) & !is.na(ref_scores$gender) &
    (!length(g_lv) | as.character(ref_scores$gender) %in% g_lv)
  map(items_list, \(items) {
    pf <- bind_cols(ref_scores %>% select("id", "age", "gender", "weight"),
                    bridge_scores(items, inst)[match(ref_scores$id, ref_ids), ])[ok, ]
    if (length(g_lv)) pf$gender <- factor(as.character(pf$gender), levels = g_lv)
    if (length(s_lv)) pf$sample <- factor(s_lv[1], levels = s_lv)
    pf
  })
}

#' Leave-own-domain-out cross-check of a directly observed scale
#'
#' Compares the norms of a directly observed scale with the norms obtained
#' when the same scale is projected through the transfer module without its
#' own bridge items (`transfer$crosscheck`). [run_norming()] writes the
#' result of all cross-check outcomes to `c1_crosscheck.csv`.
#'
#' @param direct,lodo Outputs of [norms_all_cells()] for the direct and the
#'   projected norms of the same scale.
#' @return Tibble with the cell variables (`NA` for the overall table),
#'   `raw`, `T_direct`, `T_lo_direct`, `T_hi_direct`, `T_lodo`, `T_lo_lodo`,
#'   `T_hi_lodo`, `dT` (projected minus direct) and `inside` (the projected
#'   point value lies in the band of the direct norms).
#' @export
crosscheck_table <- function(direct, lodo) {
  flat <- function(x) {
    bind_rows(if (!is.null(x$cells)) mutate(x$cells, gender = as.character(.data$gender),
                                            age_group = as.character(.data$age_group)),
              if (!is.null(x$overall)) mutate(x$overall, gender = NA_character_,
                                              age_group = NA_character_)) %>%
      mutate(raw_key = round(.data$raw, 8)) %>%
      select("gender", "age_group", "raw_key", "raw", "T", "T_lo", "T_hi")
  }
  inner_join(flat(direct), select(flat(lodo), -"raw"),
             by = c("gender", "age_group", "raw_key"), suffix = c("_direct", "_lodo")) %>%
    select(-"raw_key") %>%
    mutate(dT = .data$T_lodo - .data$T_direct,
           inside = .data$T_lodo >= .data$T_lo_direct & .data$T_lodo <= .data$T_hi_direct)
}
