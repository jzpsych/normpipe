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
