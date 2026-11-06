# normpipe 0.2.1

* `direct_frame()` coerces respondent ids to character, so that numeric ids
  (e.g. GESIS panelist ids) can be combined with the character ids of the
  Dirichlet prior rows (this failed in `bind_rows()` before).

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
