# normpipe

Bayesian norming pipeline for self-report questionnaires, generalised from the
PID-5 norming pipeline (Zimmermann et al., preprint; Zenodo code). The PID-5
family is the first use case (`inst/instruments/pid5/`); any instrument is
described by a YAML file and normed with the same code.

What the pipeline does (one instrument, several versions and scales):

| Stage | Content | Original module |
|------|---------|-----------------|
| P0 | scoring of the reference sample, raking to population margins, route per scale (direct or transfer) | data preparation |
| M1 | calibration of the item models in the norm (or a calibration) sample: GRM per facet, bifactor per domain, composite route for totals | Module 1 |
| M2 | optional: transfer models from bridge items to unobserved scales (brms), projection into the reference sample | Module 2 |
| M3a | band norms (percentile ranks and T) per reporting cell and overall, with norm bands from Bayesian bootstrap draws | Module 3 |
| M3b | continuous-age norms (Epanechnikov kernel, adaptive half width) | Module 3 |
| M3c | score-conditional EAP tables (Lord–Wingersky recursion, population prior of the cell) | Module 3 |
| M3d | integrated intervals (item-parameter draws × norm-level draws) | Module 3 |
| X1 | exports: long norm table, platform CSV per version, portal JSON | export_platform_tables |

Design decisions (agreed with J. Zimmermann, 2026-11-04):

- **No fixed response format.** Categories, values (including reverse keying,
  non-uniform values such as 0/1/2/1/0 and half-unit values) and the number
  of items come from the instrument YAML. Attainable scores follow from the
  value sets by convolution.
- **Direct route with Dirichlet prior as default.** Every cell gets
  `prior_mass` pseudo-observations spread uniformly over the attainable
  scores (uniform base measure, as in the preprint), then Bayesian bootstrap
  draws give norm bands. The transfer module is optional and only used for
  scales whose items the reference sample did not see.
- **Frequency-only input.** `frequency_frame()` builds the same draw frame
  from a table of counts per cell and score (weighted counts keep their
  uncertainty through an optional `n_obs`), so published frequency tables can
  be normed without raw data; individual intervals then need item parameters
  (`models`).
- **Measurement model calibrated in the norm sample** by default (`calib`
  overrides). The engine is `mirt` when installed, otherwise a built-in
  penalised marginal-ML estimator for GRM and bifactor models (slower, for
  short scales and tests).
- Item models run on value-ordered categories, so reverse-keyed items keep a
  positive slope and the category indices of the data never enter the models.

## Install and run

```r
# install.packages("remotes")
remotes::install_github("jzpsych/normpipe")
library(normpipe)
inst <- read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
scale_registry(inst)                      # every version_scale outcome
sim <- simulate_reference(inst, N = 1500) # or your own list(meta, items)
res <- run_norming(inst, ref = list(meta = sim$meta, items = sim$items),
                   out_dir = "demo-out",
                   stages = c("P0", "M1", "M3a", "M3b", "M3c", "M3d", "X1"))
```

`ref$meta` needs `id`, `age`, `gender` and optionally `weight`; `ref$items`
holds the item columns named `<prefix><item>`. Results and per-scale caches
land in `out_dir` (`norm_tables_long.csv`, `norm_tables_age_continuous.csv`,
`platform/norms_<version>.csv`, `portal/<instrument>.json`,
`x1_reliability.csv`). Stages are skipped when their cache exists.

Frequency tables instead of raw data:

```r
freq <- list(full_a = tibble(gender = "f", age_group = "[18,40)", score = 0:24, count = ...))
run_norming(inst, freq = freq, out_dir = "out", stages = c("P0", "M3a", "M3b", "X1"))
```

## Instrument YAML

See `inst/instruments/demo7/demo7.yml` (small, all scale types) and
`inst/instruments/pid5/pid5.yml` (the PID-5 family, 75 outcomes). Keys:

- `items`: `n_items` (or `ids`), `prefix`, `categories`, optional `values`,
  `reverse`, per-item `overrides`.
- `versions[].scales[]`: `id`, `label`, `level`, `items` or `components`
  (earlier scales of the version or inline `{id, items}` lists), `metric`
  (`sum`, `mean`, `mean_of_means`), `model` (`grm`, `bifactor`, `composite`,
  `none`), optional `calibrate_from: <version>_<scale>` (short forms use the
  long form's item parameters), `equal_block_slopes` (doublet bifactor
  models), `min_answered`, `observed`.
- `cells`: `gender`, `age_breaks`, `ages`, `kernel`; `uncertainty`:
  `prior_mass`, `n_boot`, `n_param_draws`, `probs`, `pctl_clamp`,
  `pctl_convention` (`cumulative` = P(X ≤ x), the default and the convention
  of the PID-5, IPO-16 and LPFS-BF norms; `midpoint` = P(X < x) + ½ P(X = x),
  the classical percentile rank, which keeps floor scores below T 50; chosen
  by J. Zimmermann on 2026-11-05 as a switch so that the published norms stay
  reproducible and a revision can move to midpoint); `measurement.engine`: `auto`,
  `mirt`, `internal`.

## Status (2026-11-04)

- Tested end to end on the simulated demo instrument (`scripts/smoke.R`):
  GRM facets, sum and mean scales, a bifactor domain, a composite total and a
  short form calibrated from the long form; 27 unit tests
  (`tests/testthat`: Lord–Wingersky vs brute force, frequency vs individual
  frames, prior mass per cell, grids, EAP monotonicity).
- PID-5 configuration checked against the original pipeline without data
  (`scripts/check_pid5_against_original.R path/to/pid5_norming_pipeline.R`):
  item keys and item order of all 75 outcomes (31 full, 31 SF, 6 BF, 7
  PID5BF+M), reverse keying, scoring of random responses (identical to
  machine precision), the 69 transfer outcomes and the 6 direct outcomes
  under the GESIS bridge all agree with the original. The 17 two-item facet
  doublets of the PID5BF+ are bridge predictors only, not normed scales.
- Transfer module (stage M2, brms with rstan) tested end to end on
  simulated data (`scripts/smoke_transfer.R`): with only four bridge items
  observed in the reference sample, the projected norms of the full-form
  scales differ from the norms of the true full responses by 0.6 T (facet)
  and 0.7 T (total) on average, and 97.5 % / 98.5 % of the true values lie
  inside the norm bands (5th to 95th percentile).
- mirt engine tested on the demo instrument (`scripts/smoke.R` with
  `measurement.engine: mirt`): 2.5 min instead of 12 min with the internal
  estimator; the GRM scales agree with the internal estimator to within
  0.5 T in every column of the norm table. The demo's bifactor domain (two
  specific factors of four items) is weakly identified: mirt's Oakes
  information matrix is singular at the default EM tolerance, so
  `fit_with_mirt()` refits with `TOL = 1e-5` (then, as a last resort, with
  the cross-product information; on the demo the EM still reports
  non-convergence after 5000 cycles at that tolerance, with the specific
  factor of one block collapsing towards zero); the resulting maximum-likelihood solution
  differs from the penalised internal one in the conditional and integrated
  intervals of that domain (about 1 T on average, up to 14 T in sparse
  tails). The PID-5 domain models (three facets of 6 to 14 items) are not
  expected to show this.
- cmdstanr backend of the transfer module tested on the same simulated
  case (CmdStan 2.36.0): 37 min including model compilation, T values
  within 0.1 to 0.4 T of the rstan run on average (MCMC noise); both
  backends work, cmdstanr stays the default when installed.
- Not yet run on real data: the PID-5 configuration.
- Raking uses `survey::rake`; margins are given with `make_margins()`.
- Licence: GPL-3, see LICENSE.md.
