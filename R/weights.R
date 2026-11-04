# ==============================================================================
# weights.R -- calibration of survey weights by raking (Supplement S5.1)
#
# The reference sample may carry design weights; they can be raked to
# population margins of the reporting-cell variables (age band x gender)
# with the survey package. Margins are supplied as a data frame with the
# cell variables and a column Freq (population counts), e.g. aggregated from
# single-year population counts with make_margins(). Weights are trimmed at
# trim_quantile and renormalized to the net sample size.
# ==============================================================================

#' Aggregate single-year population counts to the configured cells
#'
#' @param pop tibble with age (single years), gender, Freq
make_margins <- function(pop, inst) {
  vars <- cell_vars(inst)
  d <- add_age_group(pop, inst)
  if ("gender" %in% vars) d <- filter(d, .data$gender %in% inst$cells$gender)
  if ("age_group" %in% vars) d <- filter(d, !is.na(.data$age_group))
  d %>% count(across(all_of(vars)), wt = .data$Freq, name = "Freq") %>%
    as.data.frame()
}

#' Rake design weights to the population margins of the reporting cells
#'
#' @param meta tibble with id, age, gender, weight (design weight; 1 if none)
#' @return meta with column weight_raked and attribute "diagnostics"
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
