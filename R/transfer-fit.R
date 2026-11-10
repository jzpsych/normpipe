# ==============================================================================
# transfer-fit.R -- fitting transfer models (brms) and choosing the family
# ==============================================================================

#' Stan backend of the transfer models
#'
#' `transfer$backend` of the instrument; by default cmdstanr when it and
#' CmdStan are installed, otherwise rstan.
#' @noRd
transfer_backend <- function(inst) {
  inst$transfer$backend %||%
    if (requireNamespace("cmdstanr", quietly = TRUE) &&
        !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                          error = function(e) NULL))) "cmdstanr" else "rstan"
}

#' Fit a transfer model
#'
#' `transfer_model()` fits the model defined by [transfer_data()] with brms.
#' `fit_transfer()` does both steps for one outcome; its arguments override
#' the transfer specification of the outcome.
#'
#' The fit records the observed items in `fit$normpipe`, which
#' [project_outcome()] needs.
#'
#' @param td Output of [transfer_data()].
#' @param backend Stan backend (`"cmdstanr"` or `"rstan"`).
#' @param file Optional file path for caching the brms fit (brms `file`).
#' @param ... Further arguments passed to [brms::brm()] (e.g. `seed`,
#'   `refresh`, `chains`, `iter`, `cores`).
#' @return A `brmsfit`.
#' @export
transfer_model <- function(td, backend = "rstan", file = NULL, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  stopifnot(inherits(td, "normpipe_transfer_data"))
  sp <- td$spec
  args <- c(list(formula = td$formula, data = td$data, family = td$bfamily, prior = td$prior,
                 chains = sp$chains, cores = min(sp$chains, 4), iter = sp$iter, backend = backend,
                 control = list(adapt_delta = sp$adapt_delta, max_treedepth = sp$max_treedepth),
                 file = file),
            if (!is.null(td$init)) list(init = td$init),
            list(...))
  args <- args[!duplicated(names(args), fromLast = TRUE)]
  rds <- if (!is.null(file)) { if (grepl("\\.rds$", file)) file else paste0(file, ".rds") }
  cached <- !is.null(rds) && file.exists(rds)
  fit <- do.call(brms::brm, args)
  if (!cached) {
    fit$normpipe <- td$info
    if (!is.null(rds)) saveRDS(fit, rds)
  }
  fit
}

#' @rdname transfer_model
#' @param outcome Outcome name (`<version>_<scale>`), a column of `dat`.
#' @param dat Donor frame ([transfer_donor_frame()]).
#' @param s Scale specification ([get_scale()]).
#' @param inst Instrument.
#' @param exclude Bridge block ids to leave out of the predictors.
#' @param formula,sigma,phi Right-hand sides of the mean, sigma and precision
#'   models; `{bridge}` in `formula` stands for the bridge predictors.
#' @param family `"gaussian"`, `"skew_normal"` or `"beta_binomial"`.
#' @param adapt_delta,max_treedepth Sampler settings.
#' @export
fit_transfer <- function(outcome, dat, s, inst, exclude = NULL,
                         formula = NULL, sigma = NULL, family = NULL, phi = NULL,
                         adapt_delta = NULL, max_treedepth = NULL,
                         file = NULL, ...) {
  spec <- transfer_spec(inst, outcome)
  over <- list(exclude = exclude, formula = formula, sigma = sigma, phi = phi,
               adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  for (k in names(over)) if (!is.null(over[[k]])) spec[[k]] <- over[[k]]
  td <- transfer_data(outcome, dat, inst, family = family, spec = spec, s = s)
  transfer_model(td, backend = transfer_backend(inst), file = file, ...)
}

#' Convergence summary of a brms fit
#' @noRd
fit_convergence <- function(fit) {
  rh <- suppressWarnings(max(brms::rhat(fit), na.rm = TRUE))
  np <- brms::nuts_params(fit)
  mtd <- tryCatch(fit$fit@stan_args[[1]]$control$max_treedepth, error = function(e) NULL) %||% 10
  list(rhat_max = rh,
       divergent = sum(np$Value[np$Parameter == "divergent__"]),
       treedepth_hits = sum(np$Value[np$Parameter == "treedepth__"] >= mtd))
}

#' Choose the family of a transfer model by LOO
#'
#' Compares candidate fits of the same outcome by their expected log
#' predictive density (PSIS-LOO). The candidates describe the same events
#' (the censored likelihood gives the probability of the attainable score,
#' the beta-binomial that of the same summed score), so their elpd values are
#' comparable. Only candidates whose largest R-hat is at most `rhat_max` and
#' that did not saturate the tree depth are eligible; if none is, all are.
#'
#' @param fits Named list of brms fits (names: families).
#' @param rhat_max Largest R-hat of an eligible fit.
#' @return Tibble with `family`, `elpd_loo`, `se_elpd_loo`, `rhat_max`,
#'   `divergent`, `treedepth_hits`, `pareto_k_bad` (observations with
#'   k > 0.7), `eligible`, `chosen`, `elpd_diff` and `se_diff` (to the chosen).
#' @export
select_transfer_family <- function(fits, rhat_max = 1.01) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  rows <- lapply(names(fits), function(fam) {
    fit <- fits[[fam]]
    cv <- fit_convergence(fit)
    lo <- suppressWarnings(brms::loo(fit))
    list(row = tibble(family = fam, elpd_loo = lo$estimates["elpd_loo", 1],
                      se_elpd_loo = lo$estimates["elpd_loo", 2], rhat_max = cv$rhat_max,
                      divergent = cv$divergent, treedepth_hits = cv$treedepth_hits,
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
