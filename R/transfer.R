# ==============================================================================
# transfer.R -- optional stage M2: conditional transfer models
#
# Needed only for scales whose items the reference sample did not answer. A
# BRIDGE (items answered by the reference sample and by one or more donor
# samples that answered the full instrument) predicts every such scale score
# through a Bayesian distributional regression (brms): bridge scores + s(age)
# + gender + sample, with an interval-censored normal likelihood on the latent
# score scale and sigma ~ bridge severity + s(age). The posterior predictive
# draws for every reference respondent, discretized to the attainable scores,
# form the draw frame of the scale.
#
# Configuration (instrument$transfer):
#   enabled: true
#   bridge:                 # bridge predictors: named item sets
#     - {id: a2, items: [1, 2]}
#     ...
#   reference_level: quota       # donor sample used as the reference level
#   n_draws: 500                 # posterior predictive draws per respondent
#   adapt_delta: 0.95
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
#' regression fitted with brms: the bridge scores, `s(age)`, `gender` and
#' `sample` predict the score with an interval-censored normal likelihood,
#' and `sigma` depends on the mean bridge score and `s(age)`. The posterior
#' predictive draws for every reference respondent, snapped to the attainable
#' scores, form the draw frame of the scale. [run_norming()] runs these steps
#' in stage `M2`; the donor samples are supplied as
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

#' @rdname bridge_scores
#' @param outcome Outcome name (`<version>_<scale>`), a column of `dat`.
#' @param dat Donor frame with the outcome column, the bridge columns
#'   (`bridge_scores()`), `age`, `gender` and `sample` (factor; first level =
#'   reference level).
#' @param s Scale specification ([get_scale()]).
#' @param exclude Bridge ids to leave out of the predictors (e.g. the bridge
#'   items of the outcome itself).
#' @param adapt_delta,max_treedepth Sampler settings passed to brms.
#' @param file Optional file path for caching the brms fit.
#' @return `fit_transfer()`: a `brmsfit`.
#' @export
fit_transfer <- function(outcome, dat, s, inst, exclude = NULL,
                         adapt_delta = inst$transfer$adapt_delta %||% 0.95,
                         max_treedepth = 10, file = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
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
  rhs <- paste(c(preds, "s(age)", "gender", "sample"), collapse = " + ")
  f <- brms::bf(stats::as.formula(paste0("cens_lo | cens(cens_ind, cens_hi) ~ ", rhs)),
                sigma ~ br_mean + s(age))
  brms::brm(f, data = dat, family = stats::gaussian(),
            prior = brms::prior(normal(0, 1), class = b),
            chains = 4, cores = 4, iter = 2000,
            backend = inst$transfer$backend %||%
              if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan",
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
            file = file)
}

#' @rdname bridge_scores
#' @param fit Transfer model (`fit_transfer()`).
#' @param pred_frame Prediction frame of the reference sample: `id`, `age`,
#'   `gender`, `weight`, the bridge columns and `sample` (set to the
#'   reference level).
#' @param n_draws Number of posterior predictive draws per respondent.
#' @return `project_outcome()`: draw frame (see [direct_frame()]).
#' @export
project_outcome <- function(fit, pred_frame, s, n_draws = 500) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  pp <- brms::posterior_predict(fit, newdata = pred_frame, ndraws = n_draws)  # draws x persons
  n <- nrow(pred_frame)
  tibble(id = rep(pred_frame$id, times = n_draws),
         age = rep(pred_frame$age, times = n_draws),
         gender = rep(as.character(pred_frame$gender), times = n_draws),
         weight = rep(pred_frame$weight, times = n_draws),
         shape = 1,
         .draw = rep(seq_len(n_draws), each = n),
         score = snap_to_grid(as.vector(t(pp)), s))
}

#' Draw frame of a transfer outcome from the cached fit (used by run_norming)
#' @noRd
project_frame <- function(o, s, inst, out_dir) {
  fit <- readRDS(file.path(out_dir, paste0("transfer_", o, ".rds")))
  pred <- readRDS(file.path(out_dir, "m2_prediction_frame.rds"))
  set.seed(sum(utf8ToInt(o)))
  project_outcome(fit, pred, s, n_draws = inst$transfer$n_draws %||% 500)
}

#' Stage M2 of run_norming(): donor frame, transfer fits, prediction frame
#' @noRd
run_transfer_stage <- function(inst, reg, ref, ref_scores, out_dir, say = message) {
  donors <- inst$donors_data
  if (is.null(donors)) stop("Transfer is enabled but no donor samples were supplied ",
                            "(inst$donors_data <- list(meta, items)).")
  need <- reg$outcome[reg$route == "transfer"]
  if (!length(need)) return(invisible(NULL))
  dm <- as_tibble(donors$meta)
  lv <- inst$transfer$reference_level %||% levels(factor(dm$sample))[1]
  dm$sample <- factor(dm$sample, levels = c(lv, setdiff(unique(dm$sample), lv)))
  donor_frame <- bind_cols(dm, bridge_scores(donors$items, inst), score_scales(donors$items, inst, prorate = FALSE))
  if (!is.null(ref_scores))
    donor_frame <- donor_frame %>% filter(between(.data$age, min(ref_scores$age, na.rm = TRUE), max(ref_scores$age, na.rm = TRUE)))
  saveRDS(donor_frame, file.path(out_dir, "m2_donor_frame.rds"))
  pred <- bind_cols(ref_scores %>% select("id", "age", "gender", "weight"),
                    bridge_scores(ref$items, inst)) %>%
    mutate(sample = factor(lv, levels = levels(dm$sample)))
  saveRDS(pred, file.path(out_dir, "m2_prediction_frame.rds"))
  for (o in need) {
    f <- file.path(out_dir, paste0("transfer_", o))
    if (file.exists(paste0(f, ".rds"))) next
    s <- get_scale(inst, o)
    say("M2 ", o, ": fitting transfer model")
    fit_transfer(o, donor_frame, s, inst, file = f)
  }
  invisible(NULL)
}
