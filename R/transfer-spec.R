# ==============================================================================
# transfer-spec.R -- configuration of the transfer module (stage M2)
#
# The transfer module is needed only for scales whose items the reference
# sample did not answer. Files:
#   transfer-spec.R      configuration: transfer_spec(), priors, formulas
#   transfer-data.R      donor and prediction frames, bridge scores,
#                        transfer_data() (response, censoring, formula; no Stan)
#   transfer-observed.R  observed (bridge) items of an outcome, beta-binomial
#                        count structure
#   transfer-fit.R       transfer_model(), fit_transfer(), select_transfer_family()
#   transfer-project.R   projection of the reference sample, cross-check
#   transfer-stage.R     run_transfer_stage(): stage M2 of run_norming()
#
# Configuration (instrument$transfer):
#   enabled: true
#   bridge:                 # bridge blocks: item sets answered by the reference
#     - {id: a2, items: [1, 2]}   # sample and the donor samples; any length
#     ...
#   formula: "{bridge} + s(age) + gender + sample"   # mean model (right-hand side)
#   sigma: "br_mean + s(age)"                        # sigma model (censored families)
#   phi: "br_mean + s(age)"                          # precision model (beta_binomial)
#   family: gaussian | skew_normal | beta_binomial | [candidates]
#   select_rhat: 1.01            # largest R-hat of a candidate that may be chosen
#   select_rhat_chain: 1.05      # the same within single chains (imputed donors)
#   prior_b, prior_b_phi, priors # priors (see transfer_spec())
#   init: random | 0             # initial values [family default]
#   reference_level: quota       # first level of `sample` in the fit
#   projection_level: quota      # a donor sample or "average" [reference_level]
#   bootstrap_weights: true      # Bayesian bootstrap of the reference weights
#   exclude: {<outcome>: [bridge ids]}   # bridge blocks left out per outcome
#   crosscheck: [<outcome>, ...]         # directly observed outcomes that are
#                                        # also projected (leave-own-domain-out)
#   outcomes: {<outcome>: {formula, family, adapt_delta, ...}}  # overrides
#   n_draws: 500                 # posterior predictive draws per respondent
#   adapt_delta: 0.95, max_treedepth: 10, chains: 4, iter: 2000
#   backend: cmdstanr | rstan    # default: cmdstanr when installed, else rstan
# ==============================================================================

#' Families available for transfer models
#' @noRd
transfer_families <- c("gaussian", "skew_normal", "beta_binomial")

#' Keys of a transfer specification (global or per outcome)
#' @noRd
transfer_spec_defaults <- function() {
  list(formula = "{bridge} + s(age) + gender + sample",
       sigma = "br_mean + s(age)", phi = NULL,
       family = "gaussian", select_rhat = 1.01, select_rhat_chain = 1.05,
       adapt_delta = 0.95, max_treedepth = 10, chains = 4, iter = 2000,
       prior_b = NULL, prior_b_phi = "normal(0, 2.5)", priors = NULL, init = NULL)
}

#' Transfer specification of one outcome
#'
#' Merges the defaults, the `transfer` block of the instrument and the
#' outcome's entry under `transfer$outcomes` (a key given there replaces the
#' global one as a whole). For an outcome listed under `transfer$crosscheck`,
#' the bridge blocks that share items with the outcome are excluded
#' (leave-own-domain-out), unless `exclude` names them explicitly.
#'
#' Defaults that depend on the family (the prior of the mean coefficients)
#' are resolved per fit by [transfer_data()], because `family` may list
#' several candidates.
#'
#' @param inst Instrument ([read_instrument()]).
#' @param outcome Outcome name (`<version>_<scale>`).
#' @return List with `formula`, `sigma`, `phi`, `family` (one family or
#'   several candidates), `select_rhat`, `adapt_delta`, `max_treedepth`,
#'   `chains`, `iter`, `prior_b` (`NULL`: family default), `prior_b_phi`,
#'   `priors` (further priors as `{prior, class, dpar, coef}`), `init`
#'   (`NULL` = family default, `"random"` or 0) and `exclude`.
#' @export
transfer_spec <- function(inst, outcome) {
  tr <- inst$transfer %||% list()
  def <- transfer_spec_defaults()
  ov <- tr$outcomes[[outcome]] %||% list()
  bad <- setdiff(names(ov), c(names(def), "exclude"))
  if (length(bad))
    stop("transfer$outcomes$", outcome, ": unknown key(s) ", paste(bad, collapse = ", "),
         " (allowed: ", paste(c(names(def), "exclude"), collapse = ", "), ").")
  spec <- def
  for (k in intersect(names(def), names(tr))) spec[k] <- list(tr[[k]])
  for (k in intersect(names(def), names(ov))) spec[k] <- list(ov[[k]])
  fams <- unique(as.character(unlist(spec$family)))
  bad <- setdiff(fams, transfer_families)
  if (length(bad) || !length(fams))
    stop("transfer$family must be one of ", paste(transfer_families, collapse = ", "),
         " (or a list of them); got ", paste(bad, collapse = ", "))
  spec$family <- fams
  spec$phi <- spec$phi %||% spec$sigma
  spec$priors <- spec$priors %||% list()
  for (p in spec$priors)
    if (is.null(p$prior) || is.null(p$class)) stop("Each entry of transfer$priors needs 'prior' and 'class'.")
  if (!is.null(spec$init) && !(identical(spec$init, "random") || (is.numeric(spec$init) && spec$init == 0)))
    stop("transfer$init must be \"random\" or 0 (or left out for the family default).")
  excl <- ov$exclude %||% tr$exclude[[outcome]]
  if (is.null(excl) && outcome %in% unlist(tr$crosscheck))
    excl <- own_bridge(inst, outcome)
  spec$exclude <- as.character(unlist(excl))
  spec
}

#' Check the keys of the transfer block
#'
#' Unknown keys at the top of `transfer` are reported as a warning (a typo
#' would otherwise silently fall back to the default); unknown keys under
#' `transfer$outcomes` are an error in [transfer_spec()].
#' @noRd
check_transfer_config <- function(inst) {
  tr <- inst$transfer %||% list()
  known <- c(names(transfer_spec_defaults()), "enabled", "bridge", "reference_level",
             "projection_level", "bootstrap_weights", "exclude", "crosscheck", "outcomes",
             "n_draws", "backend")
  bad <- setdiff(names(tr), known)
  if (length(bad)) warning("Unknown key(s) in transfer: ", paste(bad, collapse = ", "), call. = FALSE)
  for (o in names(tr$outcomes)) transfer_spec(inst, o)
  invisible(TRUE)
}

#' Bridge blocks that share items with an outcome
#' @noRd
own_bridge <- function(inst, outcome) {
  s <- get_scale(inst, outcome)
  br <- inst$transfer$bridge
  keep <- map_lgl(br, \(b) length(intersect(b$items, s$items)) > 0)
  map_chr(br, "id")[keep]
}

#' Right-hand side of a model with the bridge placeholder filled in
#'
#' Terms of a factor with a single level in `dat` (e.g. `sample` with one
#' donor sample, `gender` in a single-gender study) are dropped, and
#' `by = <factor>` in a smooth is removed for such a factor.
#' @noRd
expand_rhs <- function(rhs, preds, dat = NULL) {
  rhs <- gsub("{bridge}", paste(preds, collapse = " + "), rhs, fixed = TRUE)
  if (!is.null(dat)) {
    for (v in c("sample", "gender")) {
      if (!v %in% names(dat) || length(unique(stats::na.omit(as.character(dat[[v]])))) > 1) next
      rhs <- gsub(paste0("(^|\\+)\\s*", v, "\\s*(?=\\+|$)"), "\\1", rhs, perl = TRUE)
      rhs <- gsub(paste0(",\\s*by\\s*=\\s*", v, "\\s*\\)"), ")", rhs)
    }
  }
  rhs <- gsub("\\+\\s*\\+", "+", rhs)
  rhs <- trimws(gsub("^\\s*\\+|\\+\\s*$", "", rhs))
  if (!nzchar(rhs)) "1" else rhs
}

#' Default prior of the mean coefficients for a family
#'
#' The logit scale of the beta-binomial needs a wider prior.
#' @noRd
default_prior_b <- function(family) if (family == "beta_binomial") "normal(0, 1.5)" else "normal(0, 1)"

#' Priors of a transfer model: defaults plus transfer$priors
#'
#' Defaults: `prior_b` (family default if `NULL`) on the mean coefficients and,
#' for the beta-binomial, `prior_b_phi` on the coefficients of log(phi) (brms
#' default: flat). An entry of `priors` replaces a default of the same class,
#' dpar and coef.
#' @noRd
transfer_priors <- function(spec, family) {
  key <- function(class, dpar = "", coef = "") paste(class, dpar %||% "", coef %||% "", sep = "|")
  pr <- list(list(prior = spec$prior_b %||% default_prior_b(family), class = "b", dpar = "", coef = ""))
  if (family == "beta_binomial")
    pr <- c(pr, list(list(prior = spec$prior_b_phi, class = "b", dpar = "phi", coef = "")))
  for (x in spec$priors) {
    k <- key(x$class, x$dpar, x$coef)
    pr <- c(Filter(\(p) key(p$class, p$dpar, p$coef) != k, pr),
            list(list(prior = x$prior, class = x$class, dpar = x$dpar %||% "", coef = x$coef %||% "")))
  }
  do.call(c, map(pr, \(p) brms::set_prior(p$prior, class = p$class, dpar = p$dpar, coef = p$coef)))
}
