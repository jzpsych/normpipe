# normpipe 0.2.5

* `transfer$family` may list several candidates, e.g.
  `[beta_binomial, gaussian]`. `run_transfer_stage()` then fits each
  (`transfer_<outcome>__<family>.rds`), and `select_transfer_family()` keeps
  the candidate with the highest PSIS-LOO expected log predictive density
  among those that converged (largest R-hat at most `transfer$select_rhat`,
  default 1.01, and no fit at the maximum tree depth). The comparison is on
  the same events (the probability of the attainable score). The choice is
  written to `m2_family_selection.csv`. Dense scales keep the censored
  candidates only. In the PID-5 donor data, the beta-binomial fails or loses
  for short outcomes whose unobserved rest is not overdispersed relative to
  the binomial (e.g. two unobserved items): the precision parameter drifts to
  infinity, and censored normal fits better (Delta ELPD +17 to +61).
* The default prior of the mean coefficients (`prior_b`) follows the family
  of each fit rather than the first configured family.
* Censored transfer fits (`gaussian`, `skew_normal`) get complete initial
  values: the intercepts at the location and log scale of the outcome, all
  other parameters at zero (lower-bounded ones at one). With only the two
  intercepts given, rstan failed for some models before sampling ("no more
  scalars to read").

# normpipe 0.2.4

* Observed bridge items are no longer predicted in any family. 0.2.3 applied
  this to the beta-binomial only; the censored families (`gaussian`,
  `skew_normal`) now also model the outcome minus its weighted observed items
  (censored at the attainable range of that rest), and `project_outcome()`
  adds the observed part to every draw. Before, the censored domain models
  predicted the full score, so items predicted themselves and 0.3-1.3 % of
  the PID-5 domain draws fell below the respondent's own observed part.
  Item weights follow the scale metric (`sum`, `mean`, `mean_of_means`).
  The rule is fixed and has no switch.
* An observed item the respondent skipped is no longer set to the value of
  its bridge predictor: it is drawn in every posterior draw from its
  distribution in the donor samples given the answered item of the pair
  (smoothed; marginal distribution when the partner is missing too), so the
  imputation uncertainty enters the norm band.
* New `transfer` keys (global or per outcome): `priors`, a list of further
  brms priors `{prior, class, dpar, coef}` that replace defaults of the same
  class/dpar/coef, and `init` (`random` or `0`; left out: family default).
* Fits record the observed-item information in `fit$normpipe`; censored fits
  cached by 0.2.3 or earlier keep their full-score meaning but should be
  refitted.

# normpipe 0.2.3

* Beta-binomial transfer models no longer model the bridge items an outcome
  contains. These items are observed in the reference sample, so the sum of
  the outcome given the bridge is partly known and underdispersed relative to
  the binomial; the beta-binomial (overdispersion only) then lets the
  precision parameter drift to infinity, and in the PID-5 production run 21
  of the 46 fits of such outcomes did not converge (every iteration at the
  maximum tree depth, R-hat up to infinity). `fit_transfer()` now models the
  sum of the unobserved items only (`trials` = their range), and
  `project_outcome()` adds the respondent's observed bridge items back. A
  missing observed item takes the value of its bridge predictor (the answered
  item of the pair). Outcomes without bridge items and the
  leave-own-domain-out fits are unchanged.
* Beta-binomial fits start all parameters at zero (`init = 0`; from random
  values a chain could stay on the flat plateau of large phi) and use a
  proper prior on the coefficients of log(phi), `transfer$prior_b_phi`
  (default `normal(0, 2.5)`; brms default: flat).
* `bridge_scores()` also returns the scored single bridge items (`bi_<item>`,
  integer units).
* The fit records the observed items (`fit$normpipe`); cached fits of 0.2.2
  keep their meaning (full sum) when projected.

# normpipe 0.2.2

## Beta-binomial transfer models

* `transfer$family: beta_binomial` models the summed score of an outcome as
  a count out of its attainable range (`y_sum | trials(y_trials)`), with a
  logit mean model (`transfer$formula`) and a precision model
  (`transfer$phi`, default: the sigma model). It replaces the
  interval-censored normal for bounded scales: the residual skew of a bounded
  scale changes sign with severity, which a latent normal or skew normal
  with constant shape cannot represent. In the PID-5 donor pilot the
  beta-binomial improved the expected log predictive density of the long
  facets by 30 to 240 points, matched the floor shares and the upper
  quantiles of the donors within one summed point, and ran as fast as the
  censored model.
* The family needs a count structure (sum or mean of items, or means of
  equally long components). A dense scale (mean of facet means of unequal
  length, tabulated on the reporting grid, e.g. the PID-5 domains) keeps
  the censored normal via a per-outcome `family: gaussian`: modelled as a
  count out of the grid steps the beta-binomial loses the floor and widens
  the upper tail (PID-5 domain pilot).
* `project_outcome()` converts the predicted counts back to the raw metric;
  `projection_level = "average"` draws from the beta-binomial via
  `extraDistr::rbbinom()` (new in Suggests).
* `transfer_spec()` returns `phi`, and the default `prior_b` is
  `normal(0, 1.5)` for the beta-binomial (logit scale).

# normpipe 0.2.1

* `direct_frame()` coerces respondent ids to character, so that numeric ids
  (e.g. GESIS panelist ids) can be combined with the character ids of the
  Dirichlet prior rows (this failed in `bind_rows()` before).
* PID-5 instrument: `n_boot` raised from 200 to 1000. With 200 draws the
  Monte Carlo error of the point norm reached 1 to 2 T in the sparse tails
  of the directly observed scales.

# normpipe 0.2.0

## Transfer module

* The mean model, the sigma model and the latent distribution of the transfer
  models are configurable (`transfer$formula`, `transfer$sigma`,
  `transfer$family` with `gaussian` or `skew_normal`), globally and per
  outcome (`transfer$outcomes`); `{bridge}` in a formula stands for the bridge
  predictors. Example: `"{bridge} + s(age, by = gender) + gender + sample"`.
  New `transfer_spec()` shows the specification used for an outcome.
* Bayesian bootstrap of the reference weights in every posterior draw
  (`transfer$bootstrap_weights`, default `true`): the norm bands of projected
  scales now contain the sampling error of the reference sample, as those of
  directly observed scales do. Set `false` to reproduce normpipe 0.1.0.
* Projection level as a parameter (`transfer$projection_level`, argument
  `level` of `project_outcome()`): any donor sample, or `"average"` (mean of
  the sample effects).
* Leave-own-domain-out cross-check: `transfer$exclude` passes excluded bridge
  predictors per outcome to `fit_transfer()`; outcomes listed in
  `transfer$crosscheck` are directly observed but additionally projected
  without their own bridge items, and `run_norming()` writes the comparison
  with the direct norms to `c1_crosscheck.csv` (`crosscheck_table()`).
* Reference respondents without age, or with a gender category absent from
  the donor samples, are left out of the projection with a message instead
  of failing inside brms.
* The default backend is cmdstanr only when CmdStan is actually installed.

## Weighting

* `rake_to_margins()` rakes to several margins at once (a list of margin data
  frames, any variables of `meta`, e.g. education or region next to age band
  x gender).
* After trimming, the weights are raked again and trimmed again until the
  margins hold (`rerake = TRUE`, default); `rerake = FALSE` reproduces the
  one-shot trimming of 0.1.0. New diagnostics: rounds, share of capped
  weights, largest margin deviation, number of rows not raked.
* `incomplete = "keep"` keeps respondents who cannot be raked (missing age,
  gender outside the cells) with a design-based weight, so that they enter
  the overall norms.
* New `weighting` block of the instrument (`trim_quantile`, `rerake`,
  `incomplete`), used by `run_norming()`.

## Other

* `scripts/smoke_transfer.R` compares the norms with a true percentile rank of
  the same convention (it compared midpoint norms with P(X <= x) before,
  which biased the reported deviation and coverage).
