# normpipe

Population-referenced norms with two layers of uncertainty for self-report
questionnaires.

normpipe turns a weighted reference sample (or published frequency tables)
into norm tables for any questionnaire that is described by a YAML
configuration file. It reports two kinds of uncertainty:

- **uncertainty of the norm values**: every percentile rank and T value comes
  with a credible band from a Bayesian bootstrap under a Dirichlet prior over
  the attainable scores;
- **measurement error of the individual score**: score-conditional and
  integrated intervals on the T metric from calibrated item response models
  (graded response, bifactor and composite models). The integrated interval
  also contains the uncertainty of the norm values.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jzpsych/normpipe")
```

Optional packages:

| Package | Needed for |
|---|---|
| `mirt` | fast calibration of the item models (otherwise a built-in, slower estimator is used) |
| `brms` with `cmdstanr` or `rstan` | the transfer module (stage `M2`); cmdstanr is used when installed and needs CmdStan (`cmdstanr::install_cmdstan()`) |
| `jsonlite` | the JSON export |

## Quick start

```r
library(normpipe)
inst <- read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
scale_registry(inst)                       # every <version>_<scale> outcome
sim <- simulate_reference(inst, N = 1500)  # or your own list(meta, items)
res <- run_norming(inst, ref = list(meta = sim$meta, items = sim$items),
                   out_dir = "demo-out",
                   stages = c("P0", "M1", "M3a", "M3b", "M3c", "M3d", "X1"))
```

`ref$meta` needs `id`, `age`, `gender` and optionally `weight`; `ref$items`
holds the item columns named `<prefix><item>`. Results and per-scale caches
are written to `out_dir`. Every cached result knows what it was computed
from: a rerun skips what is unchanged and recomputes what depends on a
changed setting or data set (`cache_status()`). Single stages can be run
with `norming_context()` and `run_stage()`.

Norms from frequency tables instead of raw data:

```r
freq <- list(full_b = data.frame(gender = "f", age_group = "[18,40)",
                                 score = 0:24, count = c(...)))
run_norming(inst, freq = freq, out_dir = "out", stages = c("P0", "M3a", "X1"))
```

The getting-started vignette (`vignette("normpipe")`) walks through the
pipeline step by step; `vignette("instrument-config")` documents every key of
the configuration file.

## Pipeline stages

| Stage | Content |
|---|---|
| `P0` | scoring of the reference sample, optional raking to population margins, route per scale (direct, frequency or transfer) |
| `M1` | calibration of the item models in the reference (or a separate calibration) sample: GRM per scale, bifactor per domain, composite route for totals |
| `M2` | optional transfer module: Bayesian distributional regression (brms) from bridge items to scales the reference sample did not answer, projection into the reference sample |
| `M3a` | band norms (percentile ranks and T) per reporting cell and overall, with norm bands |
| `M3b` | continuous-age norms (Epanechnikov kernel with adaptive half-width) |
| `M3c` | score-conditional intervals (summed-score EAP, Lord–Wingersky recursion, population prior of the cell) |
| `M3d` | integrated intervals (item-parameter draws × norm-value draws) |
| `X1` | exports: long norm table, continuous-age table, compact per-version tables, JSON norm file, marginal reliabilities |

## Design

- **Any response format.** Categories, item values (including reverse keying,
  non-uniform values such as 0/1/2/1/0 and half-unit values) and the number
  of items come from the configuration. Attainable scores follow from the
  value sets by convolution; scales can be sums, means or means of component
  means.
- **Direct route with a Dirichlet prior.** Every reporting cell receives
  `prior_mass` pseudo-observations spread uniformly over the attainable
  scores; Bayesian-bootstrap draws then give the norm bands.
- **Frequency input.** `frequency_frame()` builds the same draw frame from
  counts per cell and score (weighted counts keep their uncertainty through
  an optional `n_obs` column), so published frequency tables can be normed
  without raw data. Individual intervals then need item parameters
  (`models`).
- **Measurement model calibrated in the norm sample** by default; a separate
  calibration sample (`calib`) or pre-fitted models (`models`) can be given.
  Short forms can use the item parameters of the long form
  (`calibrate_from`).
- **Percentile convention.** By default the percentile rank is the midpoint
  rank PR = P(X < x) + ½ P(X = x), whose mean over the norm sample is exactly
  50. `uncertainty: pctl_convention: cumulative` switches to the relative
  cumulative frequency PR = P(X ≤ x). T = 50 + 10 Φ⁻¹(PR).
- **Optional transfer module** for scales whose items the reference sample did
  not answer: a bridge of items answered by the reference sample and by donor
  samples predicts the missing scales.

## Output files

| File | Content |
|---|---|
| `norm_tables_long.csv` | one row per version, scale, reporting cell and raw score: percentile rank with band, T with band, score-conditional and integrated interval, extrapolation flag, effective N |
| `norm_tables_age_continuous.csv` | continuous-age norms per gender and year of age |
| `platform/norms_<version>.csv` | compact per-version tables (raw score, T, integrated interval, percentile rank, extrapolation flag) |
| `portal/<instrument>.json` | JSON norm file per instrument (German keys) |
| `x1_reliability.csv` | marginal reliability of every scale with an item model |
| `p0_scale_registry.csv` | scales with their route |

Example output of the demo instrument is in `inst/examples/`.

## Example configurations

- `inst/instruments/demo7/demo7.yml`: a small instrument with a seven-point
  format that uses every scale type (sum, mean, mean of means; GRM,
  bifactor, composite; a short form calibrated from the long form).
- `inst/instruments/demo7/demo7_transfer.yml`: the same with a transfer block.
- `inst/instruments/pid5/pid5.yml`: the Personality Inventory for DSM-5 and
  its short forms (75 scales in four versions), including a transfer block
  for a reference sample that answered only 34 bridge items.

## Testing

`tests/testthat` holds the unit tests (Lord–Wingersky recursion against brute
force, frequency versus individual input, prior mass per cell, grids, EAP
monotonicity, percentile conventions). `scripts/smoke.R` runs the whole
pipeline on simulated data for the demo instrument; `scripts/smoke_transfer.R`
runs the transfer module on simulated data and compares the projected norms
with the norms of the known full responses.

## Licence

GPL-3, see `LICENSE.md`.
