# normpipe 0.4.1

Three corrections to what a draw frame contributes to the intervals and to
the effective sample size. The norm values and their bands are unchanged.
The cached results of the stages `norms`, `intervals_conditional` and
`intervals_integrated` are recomputed on the next run; calibrated models,
transfer models and draw frames are kept.

* **Integrated intervals carry the full uncertainty of the norm values.**
  `integrated_intervals()` divides the draws of the draw frame into as many
  sets as there are item-parameter draws. So far the norm table of a set
  was computed from its pooled draws (10 of 500 with the defaults), which
  averaged most of the variability of the norm values out: the integrated
  intervals contained about a tenth of its variance and were hardly wider
  than the score-conditional ones. The norm table is now computed from
  every draw in turn, so that the mixture has one component per draw. The
  prior of a set is still pooled over its draws. Integrated intervals become
  wider where the norm band is wide relative to the measurement error, that
  is, for reliable scales in small reporting cells.
* **Score-conditional intervals take their prior from all draws.**
  `conditional_tables()` used the first draw of the frame only, which is one
  posterior predictive draw with one set of bootstrap weights (transfer
  route) or one bootstrap draw (direct route). It now pools all draws, the
  distribution to which the printed norm values refer.
* **Effective sample size.** `n_eff` in the norm tables (and in the JSON
  export) was the Kish effective N of the first draw, whose bootstrap
  weights about halve it. It is now computed from the weights of the
  respondents averaged over the draws, as in the continuous norms.

# normpipe 0.4.0

The package no longer assumes an instrument normed by gender and age band on
the T metric: the groups for which norms are tabulated and the norm scale
are part of the configuration, and stages, files and export functions are
named after what they do. Instruments and scripts written for earlier
versions keep working; the norms of an unchanged instrument are the same.

* **Stage names.** The stages are `prepare`, `calibrate`, `transfer`,
  `norms`, `continuous_norms`, `intervals_conditional`,
  `intervals_integrated` and `export`; `intervals` stands for both interval
  stages. The former names (`P0`, `M1`, `M2`, `M3a`, `M3b`, `M3c`, `M3d`,
  `X1`) are still understood by `run_norming()` and `run_stage()`.
* **Output folder.** The result tables are in `out_dir`, every cached
  intermediate result in `out_dir/cache` (with `cache_manifest.csv` and
  `superseded/`). Results are named after their content:
  `norm_tables.csv` (was `norm_tables_long.csv`),
  `norm_tables_continuous.csv` (`norm_tables_age_continuous.csv`),
  `scale_registry.csv` (`p0_scale_registry.csv`), `reliability.csv`
  (`x1_reliability.csv`), `reference_scores.rds`, `rake_diagnostics.rds`,
  `transfer_family_selection.csv` (`m2_family_selection.csv`),
  `transfer_crosscheck.csv` (`c1_crosscheck.csv`), `compact/` (`platform/`)
  and `json/` (`portal/`). Cached results: `prepare_reference`,
  `calibrate_model_<outcome>`, `norms_<outcome>`,
  `continuous_norms_<outcome>`, `intervals_conditional_<outcome>`,
  `intervals_integrated_<outcome>`, `transfer_<outcome>_crosscheck`.
* **`migrate_output()`** moves an output folder of an earlier version to
  this layout. The moved results are adopted: the next run accepts them as
  current, so that nothing is computed again.
* **Norm groups.** The new `groups` block names the variables whose groups
  the norms are tabulated for: any number of categorical variables (listed
  codes) and banded numeric variables (`breaks`), e.g.
  `by: {gender: [m, f], age: {breaks: [18, 40, 65, Inf]}}` or
  `by: {education: [low, mid, high]}`. `continuous` asks for kernel norms
  along one numeric variable (`variable`, `range`, `by`, `kernel`). The
  `cells` block of earlier versions (`gender`, `age_breaks`, `ages`,
  `kernel`) is read as before. The reference data need the variables named
  in `groups`; `gender` and `age` are no longer required otherwise.
  `add_groups()` adds the band columns (`add_age_group()` is an alias),
  `group_sources()` lists the variables. Raking, frequency input, the
  transfer module and the invariance analyses follow the configured groups.
* **Norm scale.** `norm_scale` sets the standard-score metric: a preset
  (`T`, `z`, `IQ`, `stanine`, `sten`) or `{name, mean, sd, limits, digits}`.
  The score columns of every table carry the name of the scale (`IQ`,
  `IQ_lo`, `IQ_int_hi`, ...). The default is `T`.
* **Exports.** `export: {formats: [long, compact, json], language: en}`
  selects what the `export` stage writes. `export_compact()` and
  `export_json()` replace `export_platform_tables()` and
  `export_portal_json()`, which remain as aliases. The JSON file has English
  keys by default; `language: de` (and `export_portal_json()`) writes the
  German keys of earlier versions. The compact continuous tables are
  `norms_<version>_continuous.rds`. The argument `portal_citation` of
  `run_norming()` is now `citation`.
* `weighted_norms_continuous()`: the argument `ages` is now `values`
  (`ages` is still accepted); the result has one column per `by` variable
  and a column named after the continuous variable.
* `crosscheck_table()` works for any groups; its columns are
  `<score>_direct`, `<score>_lo_direct`, `<score>_hi_direct`,
  `<score>_projected`, `<score>_lo_projected`, `<score>_hi_projected`,
  `difference` and `inside`. The cross-check fits are cached as
  `transfer_<outcome>_crosscheck`.
* Transfer module: the default formulas are built from the norm groups
  (`{bridge} + s(<numeric variables>) + <categorical variables> + sample`;
  for groups by gender and age the same formulas as before). An instrument
  without groups gets `{bridge} + sample`. The donor data need the
  variables of the norm groups.
* `invariance_screen()` screens by default every scale of the first version
  that is scored from items directly (was: scales with `level: facet`).
* `simulate_reference()` also draws the variables of the norm groups other
  than `age` and `gender`.

# normpipe 0.3.1

* **New transfer family `cumulative`.** The unobserved rest of a scale (the
  summed score in integer units, as for the beta-binomial) is modelled as an
  ordinal outcome with one threshold per attainable value (brms
  `cumulative("logit")`, `y_ord | thres(K) ~ <mean model>`). The thresholds
  give the rest distribution any shape, including a rest that is
  underdispersed relative to the binomial, where the beta-binomial drifts to
  a plateau of huge `phi`. The scale model `disc ~ br_mean + s(age)` (key
  `disc`, default: the `phi` model) has its intercept fixed at zero
  (`constant(0)`; not identified next to free thresholds) and the prior
  `prior_b_disc` (default `normal(0, 1)`). The family can be a candidate in
  `transfer$family`; its PSIS-LOO describes the same events as the
  beta-binomial and the censored families. Like the beta-binomial it needs
  the count structure of a scale, so dense scales keep the censored
  candidates.
* **Family selection: beta-binomial draws on the plateau of large `phi`.**
  Above `phi` of about 1e12 the beta-binomial probabilities that brms uses
  for the log-likelihood (extraDistr) lose precision; at 1e14 and more they
  are near one for every value, so that chains stuck on the plateau look far
  better by LOO than they are. A fit with any draw at
  `phi > 1e10` is now never eligible; `m2_family_selection.csv` gains the
  column `phi_plateau` (share of such draws), and a cached selection is
  recomputed under the new rule.

# normpipe 0.3.0

Consolidation release: one handling of missing responses, caches that know
what they depend on, and a driver whose stages can be run one at a time.

* **Caches with dependency keys.** Every cached result is recorded in
  `cache_manifest.csv` with a hash of the configuration, data and upstream
  results it depends on (`cache_key()`, `cache_get()`, `cache_status()`). A
  rerun recomputes exactly what a change affects; outdated files are moved
  to `superseded/`, not deleted. Files cached by earlier versions have no
  key and are recomputed. Item models are cached per scale
  (`m1_model_<outcome>.rds` instead of `m1_models.rds`), the reference
  scores, weights and completed data sets in `p0_reference.rds`.
* **Measurement invariance.** `run_invariance()` tests the items answered by
  the reference sample across the samples and across the reporting cells and
  screens every facet in the donor samples (`invariance_sequence()`,
  `invariance_screen()`, `invariance_structure()`, `invariance_syntax()`).
  The models (configural, scalar, strict; higher-order structure written out
  by hand, WLSMV, theta parameterization, item intercepts 0) and the rule of
  Chen (2007) for the changes in fit (`invariance_flags()`,
  `invariance_decision()`) are described in `?invariance_syntax`. Needs
  lavaan (suggested); only fit indices are written.
* **Stages.** `norming_context()` collects the inputs of a run and
  `run_stage()` runs one stage; `run_norming()` is a loop over them.
* **Donor samples are an argument** (`run_norming(donors = )`) instead of
  `inst$donors_data`, which is still accepted.
* **Missing responses.** The `missing` block now reads
  `{method: impute, n_imputations, min_answered, covariates, neighbours}`
  (`reference:` of 0.2.6 is accepted for `method:`). With `impute`,
  respondents who answered less than `min_answered` (default one half) of
  the presented items are unit non-respondents and leave the sample before
  raking; this replaces inclusion rules in analysis scripts and the
  per-scale `min_answered`, which applies to `prorate` only. `impute_items()`
  is now a nearest-neighbour hot deck within the largest scale below the
  total (typically the domain), optionally with covariates, and copies all
  skipped items of a block from the same neighbour. `answered_share()` is new.
* Item models are calibrated on all respondents with at least one response
  to the scale (skipped items missing in the likelihood) instead of complete
  cases.
* **Donor samples are completed like the reference sample.** With
  `missing: {method: impute}`, donors below `min_answered` are dropped and
  the skipped items of the others are imputed. Every transfer model is then
  fitted with one chain per completed data set (`transfer_model_multiple()`,
  fitted in parallel after one compilation), so the posterior carries the
  imputation uncertainty at the cost of one fit and no donor is lost
  scale by scale. PSIS-LOO for the family selection uses the log-likelihood
  of each chain on its own data set. Convergence of such fits is judged within
  the single chains (`select_rhat_chain`, default 1.05), because the R-hat
  over the pooled chains also counts the differences between the data sets.
* The projection builds the prediction frames of the reference sample from
  the run itself; `m2_prediction_frame*.rds` are no longer written, and the
  norm stages no longer depend on a rerun of stage `M2`.
* R CMD check runs on pull requests (GitHub Actions).

## Transfer module

* `transfer_data()` stops before any fit when the beta-binomial family is
  requested and `extraDistr` is missing (brms needs it for the
  log-likelihood and the predictions of that family).

* `R/transfer.R` is split by task: `transfer-spec.R` (configuration),
  `transfer-data.R` (donor and prediction frames, `transfer_data()`),
  `transfer-observed.R` (observed bridge items, beta-binomial count
  structure), `transfer-fit.R` (`transfer_model()`, `fit_transfer()`,
  `select_transfer_family()`), `transfer-project.R` (projection,
  cross-check) and `transfer-stage.R` (stage M2).
* New `transfer_data()`: builds one transfer model (response, censoring or
  trials, formula, brms family, priors, initial values) without Stan, so it
  can be inspected and tested. `transfer_model()` fits it; `fit_transfer()`
  does both as before. New `transfer_donor_frame()`.
* `transfer_spec()`: `family` is always the vector of candidate families
  (one or more); the family-dependent default of `prior_b` is resolved per
  fit (`prior_b` stays `NULL` unless configured). A key under
  `transfer$outcomes` replaces the global key as a whole (also lists such as
  `family` or `priors`), and unknown keys are an error.
* Bridge blocks may have any number of items; a block whose items are all
  missing gives `NA`, and block ids must be unique.
* `sample` and `gender` terms are dropped from the formulas when the donor
  data have a single level (one donor sample, single-gender data); the
  projection then ignores `projection_level`.
* Stage M2 takes the donor samples as an argument
  (`run_transfer_stage(..., donors)`) and caches every fit, candidate and
  family choice through a cache function with a fingerprint of the model
  (formula, priors, data, sampler settings, backend).

# normpipe 0.2.6

* One handling of skipped items in the reference sample. New instrument block
  `missing: {reference: impute, n_imputations: 20}` (default `prorate`, the
  behaviour so far): `impute_items()` draws every skipped item from its
  distribution in the reference sample given the answered items of its block
  (bridge predictor, otherwise shortest scale) in `n_imputations` completed
  data sets. `run_norming()` computes the directly normed scales, the bridge
  predictors and the observed part of the transfer outcomes from these data
  sets; draw `r` of a draw frame uses data set `((r - 1) mod M) + 1`, so the
  imputation uncertainty enters the norm bands. Before, directly normed
  scales were prorated and bridge predictors were means of the answered
  items, while only the observed items of transfer outcomes were imputed.
* `direct_frame()` gains `scores_imp`; `project_outcome()` accepts a list of
  prediction frames (one per completed data set) and spreads the posterior
  draws over them.
* The inclusion rule (`min_answered`) is unchanged and is applied to the
  answered items. Transfer fits are not affected; caches `frame_*`, `m3*` and
  `m2_prediction_frame` of an earlier run have to be recomputed. The PID-5
  instrument switches to `impute`.

# normpipe 0.2.5

* `transfer$family` may list several candidates, e.g.
  `[beta_binomial, gaussian]`. `run_transfer_stage()` then fits each
  (`transfer_<outcome>__<family>.rds`), and `select_transfer_family()` keeps
  the candidate with the highest PSIS-LOO expected log predictive density
  among those that converged (largest R-hat at most `transfer$select_rhat`,
  default 1.01, and no fit at the maximum tree depth). The comparison is on
  the same events (the probability of the attainable score). The choice is
  written to `m2_family_selection.csv`. Dense scales keep the censored
  candidates only. The beta-binomial can fail or lose for short outcomes
  whose unobserved rest is not overdispersed relative to the binomial (e.g.
  two unobserved items): the precision parameter drifts to infinity, and a
  censored normal fits better.
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
  adds the observed part to every draw. Before, the censored models
  predicted the full score, so items predicted themselves and a draw could
  fall below the respondent's own observed part.
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
  precision parameter drift to infinity, and such fits may not converge.
  `fit_transfer()` now models the
  sum of the unobserved items only (`trials` = their range), and
  `project_outcome()` adds the respondent's observed bridge items back. A
  missing observed item takes the value of its bridge predictor (the answered
  item of the pair). Outcomes without bridge items and the cross-check fits
  are unchanged.
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
  with constant shape cannot represent.
* The family needs a count structure (sum or mean of items, or means of
  equally long components). A dense scale (mean of component means of unequal
  length, tabulated on the reporting grid) keeps
  the censored normal via a per-outcome `family: gaussian`: modelled as a
  count out of the grid steps the beta-binomial loses the floor and widens
  the upper tail.
* `project_outcome()` converts the predicted counts back to the raw metric;
  `projection_level = "average"` draws from the beta-binomial via
  `extraDistr::rbbinom()` (new in Suggests).
* `transfer_spec()` returns `phi`, and the default `prior_b` is
  `normal(0, 1.5)` for the beta-binomial (logit scale).

# normpipe 0.2.1

* `direct_frame()` coerces respondent ids to character, so that numeric ids
  can be combined with the character ids of the
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
* Cross-check of the projection: `transfer$exclude` passes excluded bridge
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
