# ==============================================================================
# transfer.R -- optional stage M2: conditional transfer models
#
# Needed only for scales whose items the reference sample did not answer. A
# BRIDGE (items answered by the reference sample and by one or more donor
# samples that answered the full instrument) predicts every such scale score
# through a Bayesian distributional regression (brms). By default: bridge
# scores + s(age) + gender + sample, an interval-censored normal likelihood on
# the latent score scale and sigma ~ bridge severity + s(age). The posterior
# predictive draws for every reference respondent, discretized to the
# attainable scores, form the draw frame of the scale.
#
# Configuration (instrument$transfer):
#   enabled: true
#   bridge:                 # bridge predictors: named item sets
#     - {id: a2, items: [1, 2]}
#     ...
#   formula: "{bridge} + s(age) + gender + sample"   # mean model (right-hand side)
#   sigma: "br_mean + s(age)"                        # sigma model (right-hand side)
#   family: gaussian | skew_normal | beta_binomial   # likelihood (see below)
#   phi: "br_mean + s(age)"                          # precision model (beta_binomial)
#   reference_level: quota       # first level of `sample` in the fit
#   projection_level: quota      # level used for the projection; a donor
#                                # sample or "average" [reference_level]
#   bootstrap_weights: true      # Bayesian bootstrap of the reference weights
#                                # in every posterior draw
#   exclude: {<outcome>: [bridge ids]}   # bridge predictors left out per outcome
#   crosscheck: [<outcome>, ...]         # directly observed outcomes that are
#                                        # also projected (leave-own-domain-out)
#   outcomes: {<outcome>: {formula, sigma, family, adapt_delta, ...}}  # overrides
#   n_draws: 500                 # posterior predictive draws per respondent
#   adapt_delta: 0.95
#   max_treedepth: 10
#   chains: 4, iter: 2000
#   backend: cmdstanr | rstan    # default: cmdstanr when installed, else rstan
# Donor samples are attached to the instrument before run_norming() as
# inst$donors_data <- list(meta, items), with meta$sample naming the sample.
# ==============================================================================

#' Transfer module: bridge scores, transfer models and projection
#'
#' The transfer module is needed only for scales whose items the reference
#' sample did not answer. A bridge (items answered by the reference sample
#' and by one or more donor samples that answered the full instrument)
#' predicts every such scale score through a Bayesian distributional
#' regression fitted with brms. By default the bridge scores, `s(age)`,
#' `gender` and `sample` predict the score with an interval-censored normal
#' likelihood, and `sigma` depends on the mean bridge score and `s(age)`. The
#' mean model, the sigma model and the likelihood family are configurable
#' (`transfer$formula`, `transfer$sigma`, `transfer$family`, per outcome
#' under `transfer$outcomes`), e.g.
#' `"{bridge} + s(age, by = gender) + gender + sample"` for gender-specific
#' age trajectories; `{bridge}` stands for the bridge predictors.
#'
#' Families: `gaussian` (default) and `skew_normal` model the score as
#' interval-censored on a latent continuous scale, with `sigma` (and, for
#' the skew normal, `alpha`) as distributional parameters. `beta_binomial`
#' models the summed score in integer units directly,
#' `sum | trials(max sum) ~ ...` with a logit link on the mean and the
#' precision `phi` modelled by `transfer$phi` (default: the sigma formula).
#' It needs no censoring, predicts on the attainable scores and has the
#' sign-changing skewness of a bounded scale built in (right-skewed near the
#' floor, left-skewed near the ceiling). It is available for scales whose
#' attainable scores form an arithmetic sequence (`sum`, `mean`, and
#' `mean_of_means` with equal component lengths), not for dense scales.
#' The posterior predictive draws for every reference respondent, snapped to
#' the attainable scores, form the draw frame of the scale. [run_norming()] runs
#' these steps in stage `M2`; the donor samples are supplied as
#' `inst$donors_data <- list(meta, items)` with a column `meta$sample`.
#'
#' The Stan backend is `transfer$backend` of the instrument, by default
#' cmdstanr when installed and rstan otherwise. cmdstanr needs an installed
#' CmdStan (`cmdstanr::install_cmdstan()`).
#'
#' @param items_df Data frame with the bridge item columns.
#' @param inst Instrument ([read_instrument()]) with a `transfer` block.
#' @return `bridge_scores()`: tibble with one column `br_<id>` per bridge
#'   predictor (mean item value) and `br_mean` (mean over all bridge items).
#' @export
bridge_scores <- function(items_df, inst) {
  br <- inst$transfer$bridge
  stopifnot(!is.null(br))
  idx <- item_indices(items_df, inst, items = unique(unlist(map(br, "items"))))
  vals <- item_values(idx, inst, items = unique(unlist(map(br, "items"))))
  unit <- attr(inst$items, "unit")
  out <- as_tibble(setNames(map(br, \(b) {
    cols <- inst$items$col[match(b$items, inst$items$item)]
    rowMeans(as.matrix(vals[cols]), na.rm = TRUE) * unit
  }), paste0("br_", map_chr(br, "id"))))
  all_items <- unique(unlist(map(br, "items")))
  out$br_mean <- rowMeans(as.matrix(vals[inst$items$col[match(all_items, inst$items$item)]]), na.rm = TRUE) * unit
  # single bridge items (scored, integer units): the beta-binomial adds the
  # bridge items an outcome contains as known values (see fit_transfer())
  for (it in all_items) out[[paste0("bi_", it)]] <- vals[[inst$items$col[match(it, inst$items$item)]]]
  out
}

#' Transfer specification of one outcome
#'
#' Merges the defaults, the `transfer` block of the instrument and the
#' outcome's entry under `transfer$outcomes`. For an outcome listed under
#' `transfer$crosscheck`, the bridge predictors that share items with the
#' outcome are excluded (leave-own-domain-out), unless `exclude` names them
#' explicitly.
#'
#' @param inst Instrument ([read_instrument()]).
#' @param outcome Outcome name (`<version>_<scale>`).
#' @return List with `formula`, `sigma`, `phi`, `family`, `exclude`,
#'   `adapt_delta`, `max_treedepth`, `chains`, `iter`, `prior_b`,
#'   `prior_b_phi`, `priors` (further priors as `{prior, class, dpar, coef}`)
#'   and `init` (`NULL` = family default, `"random"` or 0).
#' @export
transfer_spec <- function(inst, outcome) {
  tr <- inst$transfer %||% list()
  def <- list(formula = "{bridge} + s(age) + gender + sample",
              sigma = "br_mean + s(age)", phi = NULL,
              family = "gaussian",
              adapt_delta = 0.95, max_treedepth = 10,
              chains = 4, iter = 2000, prior_b = NULL, prior_b_phi = NULL,
              priors = NULL, init = NULL, select_rhat = 1.01)
  top <- tr[intersect(names(def), names(tr))]
  spec <- modifyList(def, top)
  ov <- tr$outcomes[[outcome]]
  if (!is.null(ov)) spec <- modifyList(spec, ov)
  # family: one family, or several candidates among which run_transfer_stage()
  # chooses by LOO (see select_transfer_family()); spec$family is the first
  fams <- unique(unlist(spec$family))
  bad <- setdiff(fams, c("gaussian", "skew_normal", "beta_binomial"))
  if (length(bad) || !length(fams))
    stop("transfer$family must be gaussian, skew_normal or beta_binomial (or a list of them); got ",
         paste(bad, collapse = ", "))
  spec$families <- fams
  spec$family <- fams[1]
  spec$phi <- spec$phi %||% spec$sigma
  # regression coefficients: the logit scale of the beta-binomial needs a wider
  # prior; the default follows the family of each fit (prior_b_user: as configured)
  spec$prior_b_user <- spec$prior_b
  spec$prior_b <- spec$prior_b %||% default_prior_b(spec$family)
  # coefficients of the log precision of the beta-binomial (brms default: flat)
  spec$prior_b_phi <- spec$prior_b_phi %||% "normal(0, 2.5)"
  # further priors: list of {prior, class, dpar, coef}; they replace a default
  # of the same class/dpar/coef (e.g. {prior: "normal(0, 1)", class: b, dpar: sigma})
  spec$priors <- spec$priors %||% list()
  # initial values: NULL = family default (beta-binomial 0; censored families:
  # intercepts at the outcome's location and log scale), "random", or 0
  if (!is.null(spec$init) && !(identical(spec$init, "random") || (is.numeric(spec$init) && spec$init == 0)))
    stop("transfer$init must be \"random\" or 0 (or left out for the family default).")
  excl <- ov$exclude %||% tr$exclude[[outcome]]
  if (is.null(excl) && outcome %in% unlist(tr$crosscheck))
    excl <- own_bridge(inst, outcome)
  spec$exclude <- unlist(excl)
  spec
}

#' Bridge predictors that share items with an outcome
#' @noRd
own_bridge <- function(inst, outcome) {
  s <- get_scale(inst, outcome)
  br <- inst$transfer$bridge
  keep <- map_lgl(br, \(b) length(intersect(b$items, s$items)) > 0)
  map_chr(br, "id")[keep]
}

#' Right-hand side of the mean model with the bridge placeholder filled in
#' @noRd
expand_rhs <- function(rhs, preds) {
  rhs <- gsub("{bridge}", paste(preds, collapse = " + "), rhs, fixed = TRUE)
  rhs <- gsub("\\+\\s*\\+", "+", rhs)
  trimws(gsub("^\\s*\\+|\\+\\s*$", "", rhs))
}

#' @rdname bridge_scores
#' @param outcome Outcome name (`<version>_<scale>`), a column of `dat`.
#' @param dat Donor frame with the outcome column, the bridge columns
#'   (`bridge_scores()`), `age`, `gender` and `sample` (factor; first level =
#'   reference level).
#' @param s Scale specification ([get_scale()]).
#' @param exclude Bridge ids to leave out of the predictors (e.g. the bridge
#'   items of the outcome itself); default from [transfer_spec()].
#' @param formula,sigma Right-hand sides of the mean and the sigma model;
#'   `{bridge}` in `formula` stands for the bridge predictors. Default from
#'   [transfer_spec()].
#' @param family `"gaussian"`, `"skew_normal"` (latent distribution of the
#'   interval-censored likelihood) or `"beta_binomial"` (summed score).
#' @param phi Right-hand side of the precision model of the beta-binomial;
#'   default from [transfer_spec()].
#' @param adapt_delta,max_treedepth Sampler settings passed to brms.
#' @param file Optional file path for caching the brms fit.
#' @param ... Further arguments passed to [brms::brm()] (e.g. `chains`,
#'   `iter`, `seed`, `refresh`).
#' @return `fit_transfer()`: a `brmsfit`.
#' @export
fit_transfer <- function(outcome, dat, s, inst, exclude = NULL,
                         formula = NULL, sigma = NULL, family = NULL, phi = NULL,
                         adapt_delta = NULL, max_treedepth = NULL,
                         file = NULL, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  spec <- transfer_spec(inst, outcome)
  exclude <- exclude %||% spec$exclude
  formula <- formula %||% spec$formula
  sigma <- sigma %||% spec$sigma
  phi <- phi %||% spec$phi
  family <- match.arg(family %||% spec$family, c("gaussian", "skew_normal", "beta_binomial"))
  adapt_delta <- adapt_delta %||% spec$adapt_delta
  max_treedepth <- max_treedepth %||% spec$max_treedepth
  preds <- setdiff(grep("^br_", names(dat), value = TRUE), c("br_mean", paste0("br_", exclude)))
  r2 <- summary(stats::lm(stats::reformulate(preds, outcome), data = dat))$r.squared
  if (r2 > 0.999)
    stop("Outcome '", outcome, "' is a (near-)exact linear combination of the ",
         "bridge predictors (R^2 = ", round(r2, 5), "); exclude its own bridge items.")
  # observed items: only the rest of the score is modelled (see observed_items())
  known <- observed_items(s, inst, exclude)
  info <- observed_info(s, inst, known, dat)
  h <- s$step; lo <- s$range[1]; hi <- s$range[2]
  y <- dat[[outcome]]
  y_full <- y
  if (length(known) && family != "beta_binomial") {
    y <- y - observed_part(dat, info)
    lo <- info$rest[["min"]]; hi <- info$rest[["max"]]
  }
  eps <- 1e-9
  dat <- dat %>%
    mutate(cens_ind = case_when(is.na(y) ~ NA_character_,
                                y <= lo + eps ~ "left", y >= hi - eps ~ "right", TRUE ~ "interval"),
           cens_lo = case_when(.data$cens_ind == "left" ~ lo + h / 2,
                               .data$cens_ind == "right" ~ hi - h / 2,
                               TRUE ~ y - h / 2),
           cens_hi = case_when(.data$cens_ind == "interval" ~ y + h / 2, TRUE ~ .data$cens_lo))
  rhs <- expand_rhs(formula, preds)
  if (family == "beta_binomial") {
    # Only the sum of the outcome's unobserved items is modelled. Modelling the
    # full sum would also be wrong numerically: given the bridge, a sum that
    # contains observed items is underdispersed relative to the binomial, which
    # a beta-binomial cannot represent (phi drifts to infinity, no convergence).
    bb <- bb_scale(s, known = known, values = info$range_units)
    kn <- observed_part(dat, info, wt = rep(1, length(known)))
    dat$y_sum <- as.integer(round(y_full / bb$unit_per_sum - kn)) - bb$offset
    dat$y_trials <- bb$trials
    f <- brms::bf(stats::as.formula(paste0("y_sum | trials(y_trials) ~ ", rhs)),
                  stats::as.formula(paste0("phi ~ ", phi)))
    fam <- brms::beta_binomial()
    # All parameters start at zero on the unconstrained scale (phi = 1, mean at
    # the logit of 0.5 for centred predictors). From the default random values
    # a chain can start on the plateau of large phi, where the likelihood is
    # flat in phi, and never leave it (every iteration at the maximum tree
    # depth); see also the proper prior on the phi coefficients below.
    init_fun <- 0
  } else {
    f <- brms::bf(stats::as.formula(paste0("cens_lo | cens(cens_ind, cens_hi) ~ ", rhs)),
                  stats::as.formula(paste0("sigma ~ ", sigma)))
    fam <- switch(family, gaussian = stats::gaussian(), skew_normal = brms::skew_normal())
    # Initial values: the intercepts start at the location and log scale of the
    # outcome, all other parameters at zero (lower-bounded ones at one).
    # Without them, the censored likelihood of an outcome far from zero (e.g. a
    # sum score) can be -Inf at every random initial value, notably for the
    # skew normal, whose lcdf underflows early. The list is complete: rstan
    # fails on some models when only part of the parameters is given
    # ("no more scalars to read").
    y_mid <- (dat$cens_lo + dat$cens_hi) / 2
    init_fun <- full_init(f, dat, fam,
                          list(Intercept = mean(y_mid, na.rm = TRUE),
                               Intercept_sigma = log(stats::sd(y_mid, na.rm = TRUE))))
  }
  if (!is.null(spec$init)) init_fun <- if (identical(spec$init, "random")) NULL else 0
  dots <- list(...)
  args <- c(list(formula = f, data = dat, family = fam,
                 prior = transfer_priors(spec, family),
                 chains = spec$chains, cores = min(spec$chains, 4), iter = spec$iter,
                 backend = inst$transfer$backend %||%
                   if (requireNamespace("cmdstanr", quietly = TRUE) &&
                       !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                                         error = function(e) NULL))) "cmdstanr" else "rstan",
                 control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                 file = file),
            if (!is.null(init_fun)) list(init = init_fun),
            dots)
  args <- args[!duplicated(names(args), fromLast = TRUE)]
  rds <- if (!is.null(file)) { if (grepl("\\.rds$", file)) file else paste0(file, ".rds") }
  cached <- !is.null(rds) && file.exists(rds)
  fit <- do.call(brms::brm, args)
  # record the observed items with a new fit (a cached fit without this
  # record modelled the full score and keeps that meaning in the projection)
  if (!cached) {
    fit$normpipe <- info
    if (!is.null(rds)) saveRDS(fit, rds)
  }
  fit
}

#' @rdname bridge_scores
#' @param fit Transfer model (`fit_transfer()`).
#' @param pred_frame Prediction frame of the reference sample: `id`, `age`,
#'   `gender`, `weight`, the bridge columns and `sample`; or a list of such
#'   frames, one per completed data set of the reference sample
#'   ([impute_items()]), over which the posterior draws are spread evenly.
#' @param n_draws Number of posterior predictive draws per respondent.
#' @param level Level of `sample` at which the reference sample is projected:
#'   a donor sample, or `"average"` (the mean of the sample effects over all
#'   donor samples). Default: the level already set in `pred_frame`.
#' @param bootstrap_weights If `TRUE`, the reference weights are replaced in
#'   every posterior draw by `weight x Exp(1)` (Bayesian bootstrap), so that
#'   the norm band contains the sampling uncertainty of the reference sample
#'   as well as the uncertainty of the transfer model.
#' @return `project_outcome()`: draw frame (see [direct_frame()]).
#' @export
project_outcome <- function(fit, pred_frame, s, n_draws = 500, level = NULL,
                            bootstrap_weights = FALSE) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  lv_all <- levels(fit$data$sample) %||% unique(as.character(fit$data$sample))
  is_bb <- stats::family(fit)$family == "beta_binomial"
  info <- fit$normpipe                     # NULL for fits made before normpipe 0.2.3
  if (!is.null(info) && length(info$known_items) && is.null(info$cond))
    stop("This transfer fit was made by normpipe 0.2.3 (observed items without imputation ",
         "model). Move the cached fit away and refit it with normpipe >= 0.2.4.")
  known <- info$known_items
  bb <- if (is_bb) bb_scale(s, known = known, values = info$range_units)
  if (!is.null(level) && !identical(level, "average") && !level %in% lv_all)
    stop("Projection level '", level, "' is not a donor sample (",
         paste(lv_all, collapse = ", "), ").")
  # one prediction frame, or one per completed data set of the reference
  # sample (impute_items()): the posterior draws are spread over the frames
  frames <- if (inherits(pred_frame, "data.frame")) list(pred_frame) else pred_frame
  ids <- sample(seq_len(brms::ndraws(fit)), min(n_draws, brms::ndraws(fit)))
  M <- min(length(frames), length(ids))
  chunk <- split(ids, rep(seq_len(M), length.out = length(ids)))
  pp <- do.call(rbind, map(seq_len(M), \(m) {
    pf <- frames[[m]]
    if (is_bb) pf$y_trials <- bb$trials
    if (identical(level, "average")) {
      x <- predict_average_level(fit, pf, levels = lv_all, ids = chunk[[m]])
    } else {
      if (!is.null(level)) pf$sample <- factor(level, levels = lv_all)
      x <- brms::posterior_predict(fit, newdata = pf, draw_ids = chunk[[m]])  # draws x persons
    }
    if (is_bb) {
      x <- x + bb$offset
      if (length(known)) x <- x + observed_part(pf, info, wt = rep(1, length(known)), D = nrow(x))
      x * bb$unit_per_sum                                       # summed units -> raw metric
    } else if (length(known)) {
      # latent rest, censored at its bounds, plus the observed part
      pmin(pmax(x, info$rest[["min"]]), info$rest[["max"]]) +
        observed_part(pf, info, D = nrow(x))
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
#' `sample` (equal weights) in every posterior draw; sigma (and alpha of the
#' skew normal) are taken from the same draw. Censoring plays no role in the
#' prediction of the latent value.
#' @noRd
predict_average_level <- function(fit, pred_frame, n_draws = NULL, levels, ids = NULL) {
  ids <- ids %||% sample(seq_len(brms::ndraws(fit)), min(n_draws, brms::ndraws(fit)))
  mu <- Reduce(`+`, map(levels, \(lv) {
    nd <- pred_frame; nd$sample <- factor(lv, levels = levels)
    brms::posterior_linpred(fit, newdata = nd, dpar = "mu", draw_ids = ids)
  })) / length(levels)
  nd <- pred_frame; nd$sample <- factor(levels[1], levels = levels)
  fam <- stats::family(fit)$family
  sig <- if (fam != "beta_binomial")
    brms::posterior_linpred(fit, newdata = nd, dpar = "sigma", transform = TRUE, draw_ids = ids)
  if (fam == "beta_binomial") {
    mu <- stats::plogis(mu)
    phi <- brms::posterior_linpred(fit, newdata = nd, dpar = "phi", transform = TRUE, draw_ids = ids)
    if (!requireNamespace("extraDistr", quietly = TRUE)) stop("extraDistr is not installed.")
    out <- extraDistr::rbbinom(length(mu), size = pred_frame$y_trials[1], alpha = mu * phi, beta = (1 - mu) * phi)
    return(matrix(out, nrow(mu), ncol(mu)))
  }
  if (fam == "gaussian") {
    out <- stats::rnorm(length(mu), mu, sig)
  } else if (fam == "skew_normal") {
    alpha <- brms::as_draws_matrix(fit, variable = "alpha")[ids, 1]
    alpha <- matrix(alpha, nrow(mu), ncol(mu))
    out <- brms::rskew_normal(length(mu), mu = mu, sigma = sig, alpha = alpha)
  } else stop("Family '", fam, "' is not supported for the average projection level.")
  matrix(out, nrow(mu), ncol(mu))
}

#' Draw frame of a transfer outcome from the cached fit (used by run_norming)
#' @noRd
project_frame <- function(o, s, inst, out_dir, suffix = "") {
  fit <- readRDS(file.path(out_dir, paste0("transfer_", o, suffix, ".rds")))
  f_imp <- file.path(out_dir, "m2_prediction_frames_imputed.rds")
  pred <- if (file.exists(f_imp)) readRDS(f_imp) else
    list(readRDS(file.path(out_dir, "m2_prediction_frame.rds")))
  ok <- !is.na(pred[[1]]$age) & !is.na(pred[[1]]$gender)
  pred <- map(pred, \(p) p[ok, ])
  if (length(pred) == 1) pred <- pred[[1]]
  set.seed(sum(utf8ToInt(o)))
  project_outcome(fit, pred, s, n_draws = inst$transfer$n_draws %||% 500,
                  level = inst$transfer$projection_level %||% NULL,
                  bootstrap_weights = inst$transfer$bootstrap_weights %||% TRUE)
}

#' Stage M2 of run_norming(): donor frame, transfer fits, prediction frame
#' @noRd
run_transfer_stage <- function(inst, reg, ref, ref_scores, out_dir, say = message, ref_imp = NULL) {
  donors <- inst$donors_data
  if (is.null(donors)) stop("Transfer is enabled but no donor samples were supplied ",
                            "(inst$donors_data <- list(meta, items)).")
  need <- reg$outcome[reg$route == "transfer"]
  cross <- intersect(as.character(unlist(inst$transfer$crosscheck)), reg$outcome[reg$route != "transfer"])
  if (!length(need) && !length(cross)) return(invisible(NULL))
  dm <- as_tibble(donors$meta)
  lv <- inst$transfer$reference_level %||% levels(factor(dm$sample))[1]
  dm$sample <- factor(dm$sample, levels = c(lv, setdiff(unique(as.character(dm$sample)), lv)))
  g_don <- unique(stats::na.omit(as.character(dm$gender)))
  g_lv <- c(intersect(inst$cells$gender, g_don), sort(setdiff(g_don, inst$cells$gender)))
  dm$gender <- factor(as.character(dm$gender), levels = g_lv)
  donor_frame <- bind_cols(dm, bridge_scores(donors$items, inst), score_scales(donors$items, inst, prorate = FALSE))
  if (!is.null(ref_scores))
    donor_frame <- donor_frame %>% filter(between(.data$age, min(ref_scores$age, na.rm = TRUE), max(ref_scores$age, na.rm = TRUE)))
  saveRDS(donor_frame, file.path(out_dir, "m2_donor_frame.rds"))
  pred_of <- function(items) bind_cols(ref_scores %>% select("id", "age", "gender", "weight"),
                    bridge_scores(items, inst)[match(ref_scores$id, ref$meta$id), ]) %>%
    mutate(gender = factor(as.character(.data$gender), levels = g_lv),
           sample = factor(lv, levels = levels(dm$sample)))
  pred <- pred_of(ref$items)
  # completed data sets of the reference sample: one prediction frame each
  f_imp <- file.path(out_dir, "m2_prediction_frames_imputed.rds")
  if (!is.null(ref_imp)) saveRDS(map(ref_imp, pred_of), f_imp) else if (file.exists(f_imp)) file.remove(f_imp)
  n_na <- sum(is.na(pred$age) | is.na(pred$gender))
  if (n_na) say("M2: ", n_na, " reference respondent(s) without age or with a gender ",
                "category absent from the donor samples are not projected")
  saveRDS(pred, file.path(out_dir, "m2_prediction_frame.rds"))
  todo <- c(set_names(need, need), set_names(cross, paste0(cross, "_lodo")))
  for (i in seq_along(todo)) {
    o <- todo[[i]]
    f <- file.path(out_dir, paste0("transfer_", names(todo)[i]))
    if (file.exists(paste0(f, ".rds"))) next
    s <- get_scale(inst, o)
    sp <- transfer_spec(inst, o)
    fams <- sp$families
    # a dense scale (no count structure) cannot be beta-binomial
    if ("beta_binomial" %in% fams && length(fams) > 1 &&
        inherits(try(bb_scale(s), silent = TRUE), "try-error"))
      fams <- setdiff(fams, "beta_binomial")
    excl <- if (length(sp$exclude)) paste0(", without ", paste(sp$exclude, collapse = "/")) else ""
    if (length(fams) == 1) {
      say("M2 ", names(todo)[i], ": fitting transfer model (", fams, excl, ")")
      fit_transfer(o, donor_frame, s, inst, family = fams, file = f)
      next
    }
    fits <- list()
    for (fam in fams) {
      ff <- paste0(f, "__", fam)
      if (!file.exists(paste0(ff, ".rds")))
        say("M2 ", names(todo)[i], ": fitting candidate ", fam, excl)
      fits[[fam]] <- fit_transfer(o, donor_frame, s, inst, family = fam, file = ff)
    }
    sel <- select_transfer_family(fits, rhat_max = sp$select_rhat %||% 1.01)
    sel$outcome <- names(todo)[i]
    saveRDS(sel, paste0(f, "_selection.rds"))
    best <- sel$family[sel$chosen]
    say("M2 ", names(todo)[i], ": ", best, " chosen by LOO (",
        paste(sprintf("%s %.1f", sel$family, sel$elpd_loo), collapse = ", "), ")")
    file.copy(paste0(f, "__", best, ".rds"), paste0(f, ".rds"), overwrite = TRUE)
  }
  sel_files <- file.path(out_dir, paste0("transfer_", names(todo), "_selection.rds"))
  sel_files <- sel_files[file.exists(sel_files)]
  if (length(sel_files))
    readr::write_csv(bind_rows(lapply(sel_files, readRDS)) %>% select("outcome", everything()),
                     file.path(out_dir, "m2_family_selection.csv"))
  invisible(NULL)
}

#' Choose the family of a transfer model by LOO
#'
#' Compares candidate fits of the same outcome (same response events: the
#' censored likelihood gives the probability of the attainable score, the
#' beta-binomial that of the same summed score) by their expected log
#' predictive density (PSIS-LOO). Only candidates whose largest R-hat is at
#' most `rhat_max` and that did not saturate the tree depth are eligible; if
#' none is, all are, and the table says so.
#'
#' @param fits Named list of brms fits (names: families).
#' @param rhat_max Largest R-hat of an eligible fit.
#' @return Tibble with `family`, `elpd_loo`, `se_elpd_loo`, `elpd_diff`,
#'   `se_diff` (to the chosen), `rhat_max`, `treedepth_hits`, `eligible`,
#'   `chosen`, `pareto_k_bad` (observations with k > 0.7).
#' @export
select_transfer_family <- function(fits, rhat_max = 1.01) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  rows <- lapply(names(fits), function(fam) {
    fit <- fits[[fam]]
    rh <- suppressWarnings(max(brms::rhat(fit), na.rm = TRUE))
    np <- brms::nuts_params(fit)
    td <- sum(np$Value[np$Parameter == "treedepth__"] >= (fit$fit@stan_args[[1]]$control$max_treedepth %||% 10))
    lo <- suppressWarnings(brms::loo(fit))
    list(row = tibble(family = fam, elpd_loo = lo$estimates["elpd_loo", 1],
                      se_elpd_loo = lo$estimates["elpd_loo", 2], rhat_max = rh,
                      treedepth_hits = td,
                      pareto_k_bad = sum(lo$diagnostics$pareto_k > 0.7)),
         pw = lo$pointwise[, "elpd_loo"])
  })
  tab <- bind_rows(lapply(rows, `[[`, "row"))
  tab$eligible <- tab$rhat_max <= rhat_max & tab$treedepth_hits == 0
  pool <- if (any(tab$eligible)) which(tab$eligible) else seq_len(nrow(tab))
  best <- pool[which.max(tab$elpd_loo[pool])]
  tab$chosen <- seq_len(nrow(tab)) == best
  d <- lapply(rows, function(r) r$pw - rows[[best]]$pw)
  tab$elpd_diff <- map_dbl(d, sum)
  tab$se_diff <- map_dbl(d, \(x) sqrt(length(x)) * stats::sd(x))
  tab
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

#' Summed-score representation of a scale for the beta-binomial family
#'
#' The attainable raw scores must form an arithmetic sequence: the summed
#' score in integer units is `raw / unit_per_sum`, the trials of the
#' beta-binomial are the range of the attainable sums, and `offset` is the
#' smallest attainable sum (non-zero when item values do not start at 0).
#'
#' A dense scale (mean of facet means of unequal length, tabulated on a
#' reporting grid) has no count structure: modelled as a count out of the
#' grid steps, the beta-binomial loses the floor and widens the upper tail
#' (PID-5 domain pilot, 2026-11). Such scales keep the censored families.
#' @noRd
bb_scale <- function(s, known = NULL, values = NULL) {
  if (isTRUE(s$dense))
    stop("Scale '", s$version, "_", s$id, "': a dense scale (mean of means of unequal ",
         "length, tabulated on the reporting grid) has no count structure; use ",
         "family 'gaussian' or 'skew_normal' for it (transfer$outcomes).")
  if (!isTRUE(all.equal(diff(s$grid), rep(s$step, length(s$grid) - 1))))
    stop("Scale '", s$version, "_", s$id, "': the beta-binomial family needs an ",
         "arithmetic sequence of attainable scores (sum, mean, or mean_of_means ",
         "with equal component lengths).")
  lo <- min(s$sums); hi <- max(s$sums)
  if (!isTRUE(all.equal(sort(unique(s$sums)), lo:hi)))
    stop("Scale '", s$version, "_", s$id, "': the attainable sums have gaps; ",
         "the beta-binomial family is not available.")
  # observed (known) items: their value range is removed from the modelled sum
  if (length(known)) {
    lo <- lo - values[["min"]]; hi <- hi - values[["max"]]
    if (hi <= lo) stop("Scale '", s$version, "_", s$id, "': all items are bridge items; ",
                       "the scale should be observed directly.")
  }
  list(unit_per_sum = s$step, offset = as.integer(round(lo)), trials = as.integer(round(hi - lo)))
}

default_prior_b <- function(family) if (family == "beta_binomial") "normal(0, 1.5)" else "normal(0, 1)"

#' Priors of a transfer model: defaults plus transfer$priors
#' @noRd
transfer_priors <- function(spec, family) {
  key <- function(class, dpar = "", coef = "") paste(class, dpar %||% "", coef %||% "", sep = "|")
  pr <- list(list(prior = spec$prior_b_user %||% default_prior_b(family), class = "b", dpar = "", coef = ""))
  if (family == "beta_binomial") pr <- c(pr, list(list(prior = spec$prior_b_phi, class = "b", dpar = "phi", coef = "")))
  for (x in spec$priors) {
    if (is.null(x$prior) || is.null(x$class)) stop("Each entry of transfer$priors needs 'prior' and 'class'.")
    k <- key(x$class, x$dpar, x$coef)
    pr <- c(Filter(\(p) key(p$class, p$dpar, p$coef) != k, pr),
            list(list(prior = x$prior, class = x$class, dpar = x$dpar %||% "", coef = x$coef %||% "")))
  }
  do.call(c, map(pr, \(p) brms::set_prior(p$prior, class = p$class, dpar = p$dpar, coef = p$coef)))
}

#' Complete initial values for a brms model
#'
#' Reads the parameters block of the Stan code and returns a function that
#' gives every parameter a value: `values` where named, otherwise zero (one for
#' parameters with a lower bound of zero), each with a small jitter.
#' @noRd
full_init <- function(f, dat, fam, values = list()) {
  code <- strsplit(brms::make_stancode(f, data = dat, family = fam), "\n")[[1]]
  sd <- brms::make_standata(f, data = dat, family = fam)
  i0 <- grep("^parameters \\{", code); i1 <- i0 + which(grepl("^\\}", code[-(1:i0)]))[1]
  decl <- trimws(sub("//.*$", "", code[(i0 + 1):(i1 - 1)]))
  decl <- decl[nzchar(decl)]
  pars <- lapply(decl, function(d) {
    m <- regmatches(d, regexec("^(real|vector)(<[^>]*>)?(\\[(.*)\\])?\\s+([A-Za-z0-9_]+);$", d))[[1]]
    if (!length(m)) stop("normpipe cannot set initial values for the declaration '", d, "'.")
    n <- if (nzchar(m[5])) as.integer(eval(parse(text = m[5]), envir = sd)) else 1L
    pos <- grepl("lower\\s*=\\s*0", m[3])
    list(name = m[6], vector = m[2] == "vector", n = n, base = if (pos) 1 else 0)
  })
  function() {
    out <- lapply(pars, function(p) {
      v <- values[[p$name]] %||% p$base
      v <- rep_len(v, p$n) + stats::rnorm(p$n, 0, 0.05)
      if (p$base == 1) v <- abs(v)
      if (p$vector) array(v, dim = p$n) else v
    })
    stats::setNames(out, map_chr(pars, "name"))
  }
}

#' Observed items of a transfer outcome
#'
#' The bridge items that belong to an outcome are answered in the reference
#' sample. The transfer model therefore predicts only the part of the score
#' that the other items contribute ("rest"), and the projection adds the
#' respondent's observed part. This holds for every family; it is not a
#' configuration option. Items of excluded bridge predictors (e.g. in the
#' leave-own-domain-out cross-check) count as unobserved.
#' @noRd
observed_items <- function(s, inst, exclude = NULL) {
  br <- inst$transfer$bridge
  used <- br[!map_chr(br, "id") %in% exclude]
  intersect(s$items, unique(unlist(map(used, "items"))))
}

#' Weight of each item value (integer units) in the raw score of a scale
#'
#' The raw score is linear in the item values: sum (weight = unit), mean
#' (unit / n) or mean of component means (unit / (n_component * n_components)).
#' @noRd
item_weights <- function(s, inst) {
  unit <- attr(inst$items, "unit")
  w <- switch(s$metric,
    sum = setNames(rep(unit, length(s$items)), s$items),
    mean = setNames(rep(unit / length(s$items), length(s$items)), s$items),
    mean_of_means = {
      C <- length(s$components)
      unlist(unname(map(s$components, \(it) setNames(rep(unit / (length(it) * C), length(it)), it))))
    })
  w[as.character(s$items)]
}

#' Everything the projection needs about the observed items (stored with the fit)
#'
#' Weights, the range of the observed part (integer units, for the
#' beta-binomial) and of the rest (raw metric, for the censored families),
#' and, per observed item, the distribution of its value given the other
#' item(s) of its bridge pair in the donor data. The latter is used to draw a
#' value for a reference respondent who skipped the item (a new value in
#' every posterior draw, so that the norm band carries this uncertainty).
#' @noRd
observed_info <- function(s, inst, known, dat = NULL) {
  w <- item_weights(s, inst)
  v <- setNames(inst$items$values[match(s$items, inst$items$item)], s$items)
  vmin <- map_dbl(v, min); vmax <- map_dbl(v, max)
  kc <- as.character(known); uc <- setdiff(names(w), kc)
  out <- list(known_items = known, w = w[kc],
              range_units = c(min = sum(vmin[kc]), max = sum(vmax[kc])),
              rest = c(min = sum(w[uc] * vmin[uc]), max = sum(w[uc] * vmax[uc])),
              cond = list())
  if (!length(known)) return(out)
  br <- inst$transfer$bridge
  for (it in known) {
    b <- br[[which(map_lgl(br, \(x) it %in% x$items))[1]]]
    partners <- setdiff(b$items, it)
    vals <- sort(unique(v[[as.character(it)]]))
    cnt <- NULL; keys <- numeric(0)
    if (!is.null(dat) && length(partners)) {
      pk <- round(rowMeans(as.matrix(dat[paste0("bi_", partners)])))
      x <- dat[[paste0("bi_", it)]]
      ok <- !is.na(pk) & !is.na(x)
      keys <- sort(unique(pk[ok]))
      cnt <- vapply(keys, \(k) tabulate(match(x[ok][pk[ok] == k], vals), length(vals)), numeric(length(vals)))
      cnt <- matrix(cnt, nrow = length(vals))          # values x keys
    }
    marg <- if (!is.null(cnt)) rowSums(cnt) else rep(1, length(vals))
    out$cond[[as.character(it)]] <- list(partners = partners, values = vals, keys = keys,
                                         prob = if (!is.null(cnt)) sweep(cnt + 0.5, 2, colSums(cnt + 0.5), `/`),
                                         marginal = (marg + 0.5) / sum(marg + 0.5))
  }
  out
}

#' Observed part of the score per draw (draws x persons)
#'
#' `wt`: weight per observed item (raw metric for the censored families,
#' 1 for the integer-unit sums of the beta-binomial). With `D = 1` and
#' complete data (fitting) the result is a vector.
#' @noRd
observed_part <- function(d, info, wt = info$w, D = 1) {
  known <- info$known_items
  if (!length(known)) return(if (D == 1) rep(0, nrow(d)) else matrix(0, D, nrow(d)))
  cols <- paste0("bi_", known)
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop("Bridge item column(s) ", paste(miss, collapse = ", "), " missing; ",
                         "build the frame with bridge_scores() of normpipe >= 0.2.3.")
  n <- nrow(d)
  tot <- matrix(0, D, n)
  for (j in seq_along(known)) {
    x <- d[[cols[j]]]
    vj <- matrix(x, D, n, byrow = TRUE)
    na <- which(is.na(x))
    if (length(na) && D > 1) {
      cj <- info$cond[[as.character(known[j])]]
      pk <- if (length(cj$partners)) round(rowMeans(as.matrix(d[paste0("bi_", cj$partners)])[na, , drop = FALSE])) else rep(NA, length(na))
      for (k in seq_along(na)) {
        col <- match(pk[k], cj$keys)
        pr <- if (!is.na(col)) cj$prob[, col] else cj$marginal
        vj[, na[k]] <- cj$values[sample.int(length(cj$values), D, TRUE, prob = pr)]
      }
    }
    tot <- tot + wt[[j]] * vj
  }
  if (D == 1) drop(tot) else tot
}
