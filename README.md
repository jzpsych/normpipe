# normpipe

Population-referenced norms with two layers of uncertainty for self-report
questionnaires.

normpipe turns a weighted reference sample (or published frequency tables)
into norm tables for any questionnaire that is described by a YAML
configuration file. The configuration names the items and scales, the
groups for which norms are tabulated (for example gender by age band, or
education) and the norm scale (T by default; z, IQ, stanine, sten or any
mean and standard deviation). The package reports two kinds of uncertainty:

- **uncertainty of the norm values**: every percentile rank and standard
  score comes with a credible band from a Bayesian bootstrap under a
  Dirichlet prior over the attainable scores;
- **measurement error of the individual score**: score-conditional and
  integrated intervals on the norm scale from calibrated item response models
  (graded response, bifactor and composite models). The integrated interval
  also contains the uncertainty of the norm values.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jzpsych/normpipe")          # current version
remotes::install_github("jzpsych/normpipe@v1.0.0")   # a fixed version, e.g. for a paper
```

`citation("normpipe")` gives the reference for the installed version.

Optional packages:

| Package | Needed for |
|---|---|
| `mirt` | fast calibration of the item models (otherwise a built-in, slower estimator is used) |
| `brms` with `cmdstanr` or `rstan` | the transfer module (stage `transfer`); cmdstanr is used when installed and needs CmdStan (`cmdstanr::install_cmdstan()`) |
| `jsonlite` | the JSON export |

## Quick start

```r
library(normpipe)
inst <- read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
scale_registry(inst)                       # every <version>_<scale> outcome
sim <- simulate_reference(inst, N = 1500)  # or your own list(meta, items)
res <- run_norming(inst, ref = list(meta = sim$meta, items = sim$items),
                   out_dir = "demo-out",
                   stages = c("calibrate", "norms", "continuous_norms",
                              "intervals", "export"))
```

`ref$meta` needs `id` (unique), the variables of the norm groups (here `age`
and `gender`) and optionally `weight`; `ref$items` holds the item columns
named `<prefix><item>`, in the row order of `ref$meta` or with an `id`
column by which the rows are matched. The result tables are written to `out_dir`, per-scale
caches to `out_dir/cache`. Every cached result knows what it was computed
from: a rerun skips what is unchanged and recomputes what depends on a
changed setting or data set (`cache_status()`). Single stages can be run
with `norming_context()` and `run_stage()`.

Norms from frequency tables instead of raw data:

```r
freq <- list(full_b = data.frame(gender = "f", age_group = "[18,40)",
                                 score = 0:24, count = c(...)))
run_norming(inst, freq = freq, out_dir = "out", stages = c("norms", "export"))
```

The getting-started vignette (`vignette("normpipe")`) walks through the
pipeline step by step; `vignette("instrument-config")` documents every key of
the configuration file.

## Pipeline stages

| Stage | Content |
|---|---|
| `prepare` | scoring of the reference sample, optional raking to population margins, route per scale (direct, frequency or transfer) |
| `calibrate` | calibration of the item models in the reference (or a separate calibration) sample: graded response model per scale, bifactor model for scales made of components, composite route for totals |
| `transfer` | optional transfer module: Bayesian distributional regressions (brms), fitted in the donor samples, of the scales the reference sample did not answer on the bridge items |
| `norms` | draws of the scores and weights of the reference sample per scale (direct route: Bayesian bootstrap; transfer route: projection of the transfer model into the reference sample), from them percentile ranks and standard scores per reporting cell and overall, with norm bands |
| `continuous_norms` | norms as a smooth function of a continuous variable such as age (Epanechnikov kernel with adaptive half-width) |
| `intervals_conditional` | score-conditional intervals (summed-score EAP, Lord–Wingersky recursion, population prior of the cell) |
| `intervals_integrated` | integrated intervals (mixture over item-parameter draws and over the norm table of every norm-value draw) |
| `export` | long norm table, continuous table, compact per-version tables, JSON norm file, marginal reliabilities |

`intervals` stands for both interval stages. Every stage runs `prepare`
first when it has not been run yet.

## Design

- **Any response format.** Categories, item values (including reverse keying,
  non-uniform values such as 0/1/2/1/0 and half-unit values) and the number
  of items come from the configuration. Attainable scores follow from the
  value sets by convolution; scales can be sums, means or means of component
  means.
- **Any norm groups and norm scale.** The `groups` block names the
  variables of the reference sample whose groups are tabulated: categorical
  variables with their codes, numeric variables cut into bands, in any
  combination, or none (overall norms only). `norm_scale` sets the metric of
  the standard scores, which also names the score columns of every table.
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
  cumulative frequency PR = P(X ≤ x). The standard score is
  mean + sd · Φ⁻¹(PR), on the T scale 50 + 10 Φ⁻¹(PR).
- **Optional transfer module** for scales whose items the reference sample did
  not answer: a bridge of items answered by the reference sample and by donor
  samples predicts the missing scales.

## Output files

| File | Content |
|---|---|
| `norm_tables.csv` | one row per version, scale, reporting cell and raw score: percentile rank with band, standard score with band, score-conditional and integrated interval, extrapolation flag, effective N, the raw score as a whole-number sum (`sum_score`) |
| `norm_tables_continuous.csv` | norms along the continuous variable, e.g. per gender and year of age |
| `compact/norms_<version>.csv` | compact per-version tables (raw score, standard score, integrated interval, percentile rank, extrapolation flag, summed score) |
| `json/<instrument>.json` | JSON norm file per instrument (English or German keys, `export$language`) |
| `reliability.csv` | marginal reliability of every scale with an item model |
| `scale_registry.csv` | scales with their route |
| `run_info.json` | what the tables were computed with: package and format version, commit, R and package versions, checksum of the configuration file, engine, seed rule, settings, raking |
| `cache/` | cached intermediate results per scale and their manifest |

`export$formats` selects the formats. `migrate_output()` moves an output
folder written by normpipe < 0.4 to this layout without recomputing
anything.

Example output of the demo instrument is in `inst/examples/` (written by
`scripts/make_examples.R`).

## Stability

From version 1.0.0 on, normpipe changes only by addition:

- **Configuration.** Keys and their defaults are kept; new keys come with a
  default that reproduces the earlier results. An unknown key is an error
  (with `options(normpipe.unknown_keys = "warn")` a warning), and
  `requires: "1.1"` in a configuration file states the version it needs.
- **Results.** Stage names, files, columns and JSON keys are kept; new ones
  are added. Read tables by column name, not by position.
- **Numbers.** The same data, configuration and package versions give the
  same numbers. Every random step runs with a seed of its own and leaves
  the random-number state of the session alone; `run_info.json` records
  what a set of tables was computed with. A correction that changes results
  is released as a new version and listed in `NEWS.md` under "Changes to
  results"; the cache then recomputes what it affects and keeps everything
  else, also across versions.
- **Functions.** Every help page states in a section "Lifecycle" whether a
  function is stable, a building block (exported for advanced use, may
  change in a minor version), experimental (the measurement invariance
  functions) or superseded (a name of an earlier version, kept as an
  alias). The folder `cache/` and unexported functions are internal.

Version numbers follow this: 1.x.0 adds features, 1.0.x corrects errors;
nothing that is stable is removed.

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
monotonicity, percentile conventions, norm groups and norm scales) and a
reference test: the demo instrument is normed with fixed settings and
compared with stored tables (`tests/testthat/_ref/`), so that no change
of results goes unnoticed, together with checks that no exported function,
default or column disappears. `scripts/smoke.R` runs the whole
pipeline on simulated data for the demo instrument; `scripts/smoke_transfer.R`
runs the transfer module on simulated data and compares the projected norms
with the norms of the known full responses.

## Licence

GPL-3, see `LICENSE.md`.
