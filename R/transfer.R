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
#   family: gaussian | skew_normal                   # latent distribution
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
#' mean model, the sigma model and the latent distribution (`gaussian` or
#' `skew_normal`) are configurable (`transfer$formula`, `transfer$sigma`,
#' `transfer$family`, per outcome under `transfer$outcomes`), e.g.
#' `"{bridge} + s(age, by = gender) + gender + sample"` for gender-specific
#' age trajectories; `{bridge}` stands for the bridge predictors. The
#' posterior predictive draws for every reference respondent, snapped to the
#' attainable scores, form the draw frame of the scale. [run_norming()] runs
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
#' @return List with `formula`, `sigma`, `family`, `exclude`, `adapt_delta`,
#'   `max_treedepth`, `chains`, `iter`, `prior_b`.
#' @export
transfer_spec <- function(inst, outcome) {
  tr <- inst$transfer %||% list()
  def <- list(formula = "{bridge} + s(age) + gender + sample",
              sigma = "br_mean + s(age)",
              family = "gaussian",
              adapt_delta = 0.95, max_treedepth = 10,
              chains = 4, iter = 2000, prior_b = "normal(0, 1)")
  top <- tr[intersect(names(def), names(tr))]
  spec <- modifyList(def, top)
  ov <- tr$outcomes[[outcome]]
  if (!is.null(ov)) spec <- modifyList(spec, ov)
  spec$family <- match.arg(spec$family, c("gaussian", "skew_normal"))
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
#' @param family `"gaussian"` or `"skew_normal"` (the latent distribution of
#'   the interval-censored likelihood).
#' @param adapt_delta,max_treedepth Sampler settings passed to brms.
#' @param file Optional file path for caching the brms fit.
#' @param ... Further arguments passed to [brms::brm()] (e.g. `chains`,
#'   `iter`, `seed`, `refresh`).
#' @return `fit_transfer()`: a `brmsfit`.
#' @export
fit_transfer <- function(outcome, dat, s, inst, exclude = NULL,
                         formula = NULL, sigma = NULL, family = NULL,
                         adapt_delta = NULL, max_treedepth = NULL,
                         file = NULL, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  spec <- transfer_spec(inst, outcome)
  exclude <- exclude %||% spec$exclude
  formula <- formula %||% spec$formula
  sigma <- sigma %||% spec$sigma
  family <- match.arg(family %||% spec$family, c("gaussian", "skew_normal"))
  adapt_delta <- adapt_delta %||% spec$adapt_delta
  max_treedepth <- max_treedepth %||% spec$max_treedepth
  preds <- setdiff(grep("^br_", names(dat), value = TRUE), c("br_mean", paste0("br_", exclude)))
  r2 <- summary(stats::lm(stats::reformulate(preds, outcome), data = dat))$r.squared
  if (r2 > 0.999)
    stop("Outcome '", outcome, "' is a (near-)exact linear combination of the ",
         "bridge predictors (R^2 = ", round(r2, 5), "); exclude its own bridge items.")
  h <- s$step; lo <- s$range[1]; hi <- s$range[2]
  y <- dat[[outcome]]
  dat <- dat %>%
    mutate(cens_ind = case_when(is.na(y) ~ NA_character_,
                                y <= lo ~ "left", y >= hi ~ "right", TRUE ~ "interval"),
           cens_lo = case_when(.data$cens_ind == "left" ~ lo + h / 2,
                               .data$cens_ind == "right" ~ hi - h / 2,
                               TRUE ~ y - h / 2),
           cens_hi = case_when(.data$cens_ind == "interval" ~ y + h / 2, TRUE ~ .data$cens_lo))
  f <- brms::bf(stats::as.formula(paste0("cens_lo | cens(cens_ind, cens_hi) ~ ",
                                         expand_rhs(formula, preds))),
                stats::as.formula(paste0("sigma ~ ", sigma)))
  fam <- switch(family, gaussian = stats::gaussian(), skew_normal = brms::skew_normal())
  # Initial values: the intercepts start at the location and log scale of the
  # outcome (other parameters random). Without them, the censored likelihood
  # of an outcome far from zero (e.g. a sum score) can be -Inf at every random
  # initial value, notably for the skew normal, whose lcdf underflows early.
  y_mid <- (dat$cens_lo + dat$cens_hi) / 2
  init_fun <- function() list(Intercept = mean(y_mid, na.rm = TRUE) + stats::rnorm(1, 0, 0.05),
                              Intercept_sigma = log(stats::sd(y_mid, na.rm = TRUE)) + stats::rnorm(1, 0, 0.05))
  dots <- list(...)
  args <- c(list(formula = f, data = dat, family = fam,
                 prior = brms::set_prior(spec$prior_b, class = "b"),
                 chains = spec$chains, cores = min(spec$chains, 4), iter = spec$iter,
                 backend = inst$transfer$backend %||%
                   if (requireNamespace("cmdstanr", quietly = TRUE) &&
                       !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                                         error = function(e) NULL))) "cmdstanr" else "rstan",
                 control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                 init = init_fun, file = file),
            dots)
  args <- args[!duplicated(names(args), fromLast = TRUE)]
  do.call(brms::brm, args)
}

#' @rdname bridge_scores
#' @param fit Transfer model (`fit_transfer()`).
#' @param pred_frame Prediction frame of the reference sample: `id`, `age`,
#'   `gender`, `weight`, the bridge columns and `sample`.
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
  if (identical(level, "average")) {
    pp <- predict_average_level(fit, pred_frame, n_draws, lv_all)
  } else {
    if (!is.null(level)) {
      if (!level %in% lv_all) stop("Projection level '", level, "' is not a donor sample (",
                                   paste(lv_all, collapse = ", "), ").")
      pred_frame$sample <- factor(level, levels = lv_all)
    }
    pp <- brms::posterior_predict(fit, newdata = pred_frame, ndraws = n_draws)  # draws x persons
  }
  draws_to_frame(pp, pred_frame, s, bootstrap_weights)
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
predict_average_level <- function(fit, pred_frame, n_draws, levels) {
  ids <- sample(seq_len(brms::ndraws(fit)), min(n_draws, brms::ndraws(fit)))
  mu <- Reduce(`+`, map(levels, \(lv) {
    nd <- pred_frame; nd$sample <- factor(lv, levels = levels)
    brms::posterior_linpred(fit, newdata = nd, dpar = "mu", draw_ids = ids)
  })) / length(levels)
  nd <- pred_frame; nd$sample <- factor(levels[1], levels = levels)
  sig <- brms::posterior_linpred(fit, newdata = nd, dpar = "sigma", transform = TRUE, draw_ids = ids)
  fam <- stats::family(fit)$family
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
  pred <- readRDS(file.path(out_dir, "m2_prediction_frame.rds"))
  pred <- pred[!is.na(pred$age) & !is.na(pred$gender), ]
  set.seed(sum(utf8ToInt(o)))
  project_outcome(fit, pred, s, n_draws = inst$transfer$n_draws %||% 500,
                  level = inst$transfer$projection_level %||% NULL,
                  bootstrap_weights = inst$transfer$bootstrap_weights %||% TRUE)
}

#' Stage M2 of run_norming(): donor frame, transfer fits, prediction frame
#' @noRd
run_transfer_stage <- function(inst, reg, ref, ref_scores, out_dir, say = message) {
  donors <- inst$donors_data
  if (is.null(donors)) stop("Transfer is enabled but no donor samples were supplied ",
                            "(inst$donors_data <- list(meta, items)).")
  need <- reg$outcome[reg$route == "transfer"]
  cross <- intersect(unlist(inst$transfer$crosscheck), reg$outcome[reg$route != "transfer"])
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
  pred <- bind_cols(ref_scores %>% select("id", "age", "gender", "weight"),
                    bridge_scores(ref$items, inst)[match(ref_scores$id, ref$meta$id), ]) %>%
    mutate(gender = factor(as.character(.data$gender), levels = g_lv),
           sample = factor(lv, levels = levels(dm$sample)))
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
    say("M2 ", names(todo)[i], ": fitting transfer model (", sp$family,
        if (length(sp$exclude)) paste0(", without ", paste(sp$exclude, collapse = "/")), ")")
    fit_transfer(o, donor_frame, s, inst, file = f)
  }
  invisible(NULL)
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
