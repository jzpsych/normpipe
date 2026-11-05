# ==============================================================================
# weights.R -- calibration of survey weights by raking
#
# The reference sample may carry design weights; they can be raked to
# population margins of the reporting-cell variables (age band x gender)
# with the survey package. Margins are supplied as a data frame with the
# cell variables and a column Freq (population counts), e.g. aggregated from
# single-year population counts with make_margins(). Weights are trimmed at
# trim_quantile and renormalized to the net sample size.
# ==============================================================================

#' Population margins and raking of survey weights
#'
#' `make_margins()` aggregates single-year population counts to the
#' reporting cells of an instrument. `rake_to_margins()` rakes the design
#' weights of the reference sample to these margins with [survey::rake()],
#' trims them at `trim_quantile` and rescales them to the net sample size.
#' [run_norming()] calls it when `margins` is given.
#'
#' @param pop Tibble with `age` (single years), `gender` and `Freq`
#'   (population counts).
#' @param inst Instrument ([read_instrument()]).
#' @return `make_margins()`: data frame with the cell variables and `Freq`.
#' @export
make_margins <- function(pop, inst) {
  vars <- cell_vars(inst)
  d <- add_age_group(pop, inst)
  if ("gender" %in% vars) d <- filter(d, .data$gender %in% inst$cells$gender)
  if ("age_group" %in% vars) d <- filter(d, !is.na(.data$age_group))
  d %>% count(across(all_of(vars)), wt = .data$Freq, name = "Freq") %>%
    as.data.frame()
}

#' @rdname make_margins
#' @param meta Tibble with `id`, `age`, `gender`, `weight` (design weight; 1
#'   if none).
#' @param margins Output of `make_margins()`.
#' @param trim_quantile Raked weights are trimmed at this quantile.
#' @return `rake_to_margins()`: `meta` (rows with complete cell variables)
#'   with column `weight_raked` and attribute `"diagnostics"` (Kish effective
#'   N before and after raking, ratio of largest to smallest weight).
#' @export
rake_to_margins <- function(meta, margins, inst, trim_quantile = 0.99) {
  vars <- cell_vars(inst)
  d <- add_age_group(meta, inst) %>% filter(!is.na(.data$weight))
  for (v in vars) d <- d[!is.na(d[[v]]), ]
  if ("gender" %in% vars) d <- filter(d, .data$gender %in% inst$cells$gender)
  d$gender <- factor(d$gender)
  margins$gender <- if ("gender" %in% names(margins)) factor(margins$gender, levels = levels(d$gender))
  des <- survey::svydesign(ids = ~1, weights = ~weight, data = d)
  f <- stats::as.formula(paste("~", paste(vars, collapse = " + ")))
  raked <- survey::rake(des, sample.margins = list(f),
                        population.margins = list(margins),
                        control = list(maxit = 50, epsilon = 1e-7))
  w <- stats::weights(raked)
  w <- pmin(w, quantile(w, trim_quantile, na.rm = TRUE))
  d$weight_raked <- w * nrow(d) / sum(w)
  attr(d, "diagnostics") <- list(
    kish_neff_design = kish_neff(d$weight),
    kish_neff_raked = kish_neff(d$weight_raked),
    weight_ratio_max = max(d$weight_raked) / min(d$weight_raked))
  d
}
