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

#' Transfer model on several completed donor data sets
#'
#' When skipped items of the donor samples were imputed
#' (`missing: {method: impute}`), the model is fitted once per completed data
#' set with one chain each and the chains are pooled, so the posterior carries
#' the imputation uncertainty at the cost of a single fit. R-hat over the
#' pooled chains also counts differences between the completed data sets, so
#' convergence of such a fit is judged by the largest split R-hat within a
#' single chain (`fit$normpipe_rhat_within`). A single chain gives a noisier
#' R-hat than four (healthy fits reach about 1.01 to 1.04, a stuck chain 1.5
#' and more), hence its own threshold in [select_transfer_family()]. The
#' pointwise log-likelihood of each chain is evaluated on its own data set;
#' the resulting PSIS-LOO object is stored in `fit$normpipe_loo` and used by
#' [select_transfer_family()].
#'
#' @param tds List of [transfer_data()] outputs, one per completed data set,
#'   for the same outcome and family.
#' @param backend Stan backend.
#' @param ... Further arguments passed to [brms::brm_multiple()].
#' @return A `brmsfit` with pooled chains.
#' @export
transfer_model_multiple <- function(tds, backend = "rstan", ...) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  stopifnot(length(tds) >= 1, all(map_lgl(tds, inherits, "normpipe_transfer_data")))
  if (length(tds) == 1) return(transfer_model(tds[[1]], backend = backend, ...))
  td <- tds[[1]]; sp <- td$spec
  # compile once, then one chain per completed data set, in parallel
  args0 <- c(list(formula = td$formula, data = td$data, family = td$bfamily, prior = td$prior,
                  chains = 0, backend = backend), list(...))
  args0 <- args0[!duplicated(names(args0), fromLast = TRUE)]
  fit0 <- suppressMessages(do.call(brms::brm, args0))
  datas <- map(tds, "data")
  # initial values are drawn here, one set per chain (0: brms default "0")
  init <- if (is.function(td$init)) map(seq_along(tds), \(i) list(td$init())) else
    if (!is.null(td$init)) rep(list(td$init), length(tds))
  control <- list(adapt_delta = sp$adapt_delta, max_treedepth = sp$max_treedepth)
  iter <- list(...)$iter %||% sp$iter
  seed <- list(...)$seed %||% sample.int(1e6, 1)
  # self-contained, so that it runs on the workers of a cluster
  run_one <- function(i, fit0, datas, iter, control, seed, init) {
    a <- list(fit0, newdata = datas[[i]], chains = 1, cores = 1, iter = iter, control = control,
              recompile = FALSE, seed = seed + i, refresh = 0)
    if (!is.null(init)) a$init <- init[[i]]
    do.call(stats::update, a)
  }
  environment(run_one) <- globalenv()
  workers <- min(length(tds), getOption("mc.cores", 4L))
  fits <- NULL
  if (workers > 1) fits <- tryCatch({
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function() { requireNamespace("brms"); NULL })
    parallel::parLapply(cl, seq_along(tds), run_one, fit0 = fit0, datas = datas, iter = iter,
                        control = control, seed = seed, init = init)
  }, error = function(e) { message("transfer_model_multiple: parallel fits failed (", conditionMessage(e),
                                   "), fitting one after the other"); NULL })
  if (is.null(fits)) fits <- map(seq_along(tds), run_one, fit0 = fit0, datas = datas, iter = iter,
                                 control = control, seed = seed, init = init)
  ll <- do.call(rbind, map(fits, brms::log_lik))
  chain <- rep(seq_along(fits), times = map_int(fits, brms::ndraws))
  loo <- suppressWarnings(loo::loo(ll, r_eff = loo::relative_eff(exp(ll), chain_id = chain)))
  # convergence within each chain (split R-hat): pooled R-hat would also
  # count the differences between the completed data sets
  rhat_within <- max(map_dbl(fits, \(f) suppressWarnings(max(brms::rhat(f), na.rm = TRUE))))
  fit <- brms::combine_models(mlist = fits, check_data = FALSE)
  fit$normpipe_rhat_within <- rhat_within
  fit$normpipe <- td$info
  fit$normpipe_loo <- loo
  fit$normpipe_n_data <- length(tds)
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
#' @param family `"gaussian"`, `"skew_normal"`, `"beta_binomial"` or `"cumulative"`.
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
  list(rhat_max = rh, rhat_within = fit$normpipe_rhat_within %||% NA_real_,
       divergent = sum(np$Value[np$Parameter == "divergent__"]),
       treedepth_hits = sum(np$Value[np$Parameter == "treedepth__"] >= mtd),
       phi_plateau = phi_plateau(fit))
}

#' Share of draws of a beta-binomial fit on the plateau of large phi
#'
#' Where the rest is not overdispersed, phi can drift to astronomically large
#' values (the beta-binomial then equals the binomial, and the likelihood is
#' flat in phi). From about phi = 1e12 on, the beta-binomial probabilities
#' that brms evaluates for the log-likelihood (extraDistr) lose precision; at
#' 1e14 and above they are near one for every value, so that LOO of such
#' draws is meaningless (family pilot 2026-11-12: chains stuck at phi ~ 1e19
#' had a log-likelihood of 0 for every donor). Counts draws in which phi
#' exceeds 1e10 for any donor; 0 for other families.
#' @noRd
phi_plateau <- function(fit) {
  if (stats::family(fit)$family != "beta_binomial") return(0)
  lp <- brms::posterior_linpred(fit, dpar = "phi")     # log(phi), draws x donors
  mean(apply(lp, 1, max) > log(1e10))
}

#' Choose the family of a transfer model by LOO
#'
#' Compares candidate fits of the same outcome by their expected log
#' predictive density (PSIS-LOO). The candidates describe the same events
#' (the censored likelihood gives the probability of the attainable score,
#' the beta-binomial and the cumulative model that of the same summed score),
#' so their elpd values are comparable. Only candidates whose largest R-hat
#' is at most `rhat_max`, that did not saturate the tree depth and that have
#' no draw on the plateau of large phi (beta-binomial, see `phi_plateau`)
#' are eligible; if none is, all are.
#'
#' @param fits Named list of brms fits (names: families).
#' @param rhat_max Largest R-hat of an eligible fit.
#' @param rhat_chain Largest within-chain split R-hat of an eligible fit that
#'   was pooled over completed donor data sets ([transfer_model_multiple()]);
#'   for such fits `rhat_max` is reported but not used, because it also counts
#'   the differences between the data sets.
#' @return Tibble with `family`, `elpd_loo`, `se_elpd_loo`, `rhat_max`,
#'   `rhat_within` (`NA` for an ordinary fit),
#'   `divergent`, `treedepth_hits`, `phi_plateau` (share of draws with
#'   phi > 1e10; their elpd is not interpretable), `pareto_k_bad`
#'   (observations with k > 0.7), `eligible`, `chosen`, `elpd_diff` and
#'   `se_diff` (to the chosen).
#' @export
select_transfer_family <- function(fits, rhat_max = 1.01, rhat_chain = 1.05) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  rows <- lapply(names(fits), function(fam) {
    fit <- fits[[fam]]
    cv <- fit_convergence(fit)
    lo <- fit$normpipe_loo %||% suppressWarnings(brms::loo(fit))
    list(row = tibble(family = fam, elpd_loo = lo$estimates["elpd_loo", 1],
                      se_elpd_loo = lo$estimates["elpd_loo", 2], rhat_max = cv$rhat_max,
                      rhat_within = cv$rhat_within,
                      divergent = cv$divergent, treedepth_hits = cv$treedepth_hits,
                      phi_plateau = cv$phi_plateau,
                      pareto_k_bad = sum(lo$diagnostics$pareto_k > 0.7)),
         pw = lo$pointwise[, "elpd_loo"])
  })
  tab <- bind_rows(lapply(rows, `[[`, "row"))
  # fits pooled over completed donor data sets: R-hat within the chains
  tab$eligible <- ifelse(is.na(tab$rhat_within), tab$rhat_max <= rhat_max, tab$rhat_within <= rhat_chain) &
    tab$treedepth_hits == 0 & tab$phi_plateau == 0
  pool <- if (any(tab$eligible)) which(tab$eligible) else seq_len(nrow(tab))
  best <- pool[which.max(tab$elpd_loo[pool])]
  tab$chosen <- seq_len(nrow(tab)) == best
  d <- lapply(rows, function(r) r$pw - rows[[best]]$pw)
  tab$elpd_diff <- map_dbl(d, sum)
  tab$se_diff <- map_dbl(d, \(x) sqrt(length(x)) * stats::sd(x))
  tab
}
