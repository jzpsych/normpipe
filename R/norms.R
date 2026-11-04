# ==============================================================================
# norms.R -- weighted percentile norms per reporting cell (Supplement S5.1)
#
# Input is a draw frame (frames.R). Per cell and draw, the weighted cumulative
# share at or below every attainable score is computed; across draws the
# median gives the norm value and the 2.5th/97.5th percentiles the norm band.
# T = 50 + 10 * qnorm(percentile), with the percentile clamped at pctl_clamp.
# ==============================================================================

TICK_TOL <- 1e-8

#' Weighted empirical CDF evaluated at the grid points
wecdf_at <- function(score, weight, grid) {
  o <- order(score)
  cw <- cumsum(weight[o])
  c(0, cw)[findInterval(grid + TICK_TOL, score[o]) + 1] / cw[length(cw)]
}

pctl_to_T <- function(p, clamp = 1e-4) 50 + 10 * qnorm(pmin(pmax(p, clamp), 1 - clamp))

#' Band norms of one scale: cells = gender x age band (as configured) and overall
#'
#' @param fr draw frame
#' @param s scale specification
#' @param inst instrument
#' @param cells character vector of cell variables; character(0) = overall
#' @return tibble with the cell variables, raw, pctl_med, pctl_lo, pctl_hi,
#'   T, T_lo, T_hi, n_eff (Kish effective N of the point frame)
weighted_norms <- function(fr, s, inst, cells = cell_vars(inst),
                           probs = inst$uncertainty$probs,
                           clamp = inst$uncertainty$pctl_clamp) {
  n_bad <- sum(is.na(fr$score) | is.na(fr$weight))
  if (n_bad > 0) stop("weighted_norms(): ", n_bad, " rows with NA score or weight.")
  grid <- s$grid
  fr <- add_age_group(fr, inst)
  if ("gender" %in% cells) fr <- filter(fr, .data$gender %in% inst$cells$gender)
  if ("age_group" %in% cells) fr <- filter(fr, !is.na(.data$age_group))
  if (!length(cells)) fr$.cell <- "overall"
  grp <- if (length(cells)) cells else ".cell"
  neff <- fr %>% filter(.data$.draw == min(.data$.draw), .data$shape >= 1) %>%
    group_by(across(all_of(grp))) %>%
    summarise(n_eff = sum(.data$weight)^2 / sum(.data$weight^2), n = n(), .groups = "drop")
  out <- fr %>%
    group_by(across(all_of(grp)), .data$.draw) %>%
    reframe(raw = grid, pctl = wecdf_at(.data$score, .data$weight, grid)) %>%
    group_by(across(all_of(grp)), .data$raw) %>%
    summarise(pctl_med = median(.data$pctl),
              pctl_lo = quantile(.data$pctl, probs[1]),
              pctl_hi = quantile(.data$pctl, probs[2]),
              .groups = "drop") %>%
    mutate(T = pctl_to_T(.data$pctl_med, clamp),
           T_lo = pctl_to_T(.data$pctl_lo, clamp),
           T_hi = pctl_to_T(.data$pctl_hi, clamp)) %>%
    left_join(neff, by = grp)
  if (!length(cells)) out$.cell <- NULL
  out
}

#' Band norms for all configured cells plus the overall table
norms_all_cells <- function(fr, s, inst) {
  cv <- cell_vars(inst)
  res <- list()
  if (length(cv)) res$cells <- weighted_norms(fr, s, inst, cells = cv)
  if (isTRUE(inst$cells$overall) || !length(cv))
    res$overall <- weighted_norms(fr, s, inst, cells = character(0))
  res
}

#' Continuous-age kernel norms (Supplement S5.1)
#'
#' Per gender (or pooled when the instrument has no gender cells) and target
#' age, respondents are weighted by an Epanechnikov kernel around the target
#' age whose half-width widens from half_width0 until the local effective
#' sample size reaches min_neff, up to max_half_width. Pseudo-observations of
#' the prior carry their own ages (prior_rows), so they enter the kernel like
#' everybody else.
weighted_norms_continuous <- function(fr, s, inst, ages = inst$cells$ages,
                                      probs = inst$uncertainty$probs,
                                      clamp = inst$uncertainty$pctl_clamp) {
  k <- inst$cells$kernel
  grid <- s$grid
  fr <- filter(fr, !is.na(.data$age))
  genders <- if (is.null(inst$cells$gender)) NA_character_ else inst$cells$gender
  map_dfr(genders, \(g) {
    pg <- if (is.na(g)) fr else filter(fr, .data$gender == g)
    draws <- sort(unique(pg$.draw))
    persons <- pg %>% group_by(.data$id) %>%
      summarise(age = first(.data$age), w_bar = mean(.data$weight), .groups = "drop")
    map_dfr(ages, \(a) {
      hw <- k$half_width0
      repeat {
        kern <- pmax(0, 1 - ((persons$age - a) / hw)^2)
        wk <- persons$w_bar * kern
        neff <- if (sum(wk) > 0) sum(wk)^2 / sum(wk^2) else 0
        if (neff >= k$min_neff || hw >= k$max_half_width) break
        hw <- hw + 1
      }
      kv <- setNames(kern, persons$id)
      pa <- pg %>% mutate(w = .data$weight * kv[as.character(.data$id)]) %>%
        filter(.data$w > 0)
      P <- sapply(draws, \(d) { x <- pa[pa$.draw == d, ]; wecdf_at(x$score, x$w, grid) })
      P <- matrix(P, nrow = length(grid))
      tibble(gender = g, age = a, half_width = hw, n_eff_local = neff, raw = grid,
             pctl_med = apply(P, 1, median),
             pctl_lo = apply(P, 1, quantile, probs[1]),
             pctl_hi = apply(P, 1, quantile, probs[2])) %>%
        mutate(T = pctl_to_T(.data$pctl_med, clamp),
               T_lo = pctl_to_T(.data$pctl_lo, clamp),
               T_hi = pctl_to_T(.data$pctl_hi, clamp))
    })
  })
}

#' Kish effective sample size
kish_neff <- function(w) sum(w)^2 / sum(w^2)
