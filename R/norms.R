# ==============================================================================
# norms.R -- weighted percentile norms per reporting cell
#
# Input is a draw frame (frames.R). Per cell and draw, the weighted percentile
# rank of every attainable score is computed; across draws the median gives
# the norm value and the quantiles `probs` (default 2.5 % and 97.5 %) the norm
# band. The standard score is mean + sd * qnorm(percentile rank) on the norm
# scale of the instrument (default T: 50 + 10 z), with the percentile rank
# clamped to [pctl_clamp, 1 - pctl_clamp]. The score columns carry the name
# of the norm scale (`T`, `T_lo`, `T_hi` by default).
#
# Percentile convention (uncertainty$pctl_convention):
#   "midpoint"    PR = P(X < x) + P(X = x) / 2, the classical percentile rank
#                 (Lienert & Raatz): the mean percentile rank of the norm
#                 sample is exactly 50 and floor scores are not shifted
#                 upwards (default)
#   "cumulative"  PR = P(X <= x), the relative cumulative frequency (option,
#                 e.g. to reproduce tables built with that convention)
# ==============================================================================

TICK_TOL <- 1e-8

#' Weighted percentile ranks at the grid points
#' @param convention "midpoint" (P(X < x) + P(X = x) / 2, default) or "cumulative" (P(X <= x))
#' @noRd
wecdf_at <- function(score, weight, grid, convention = "midpoint") {
  convention <- match.arg(convention, c("midpoint", "cumulative"))
  o <- order(score)
  cw <- cumsum(weight[o]); tot <- cw[length(cw)]
  at_or_below <- c(0, cw)[findInterval(grid + TICK_TOL, score[o]) + 1] / tot
  if (convention == "cumulative") return(at_or_below)
  below <- c(0, cw)[findInterval(grid - TICK_TOL, score[o]) + 1] / tot
  (at_or_below + below) / 2
}

#' Band norms of one scale with norm bands
#'
#' Percentile ranks and standard scores for every attainable score of a
#' scale, per reporting cell (the combinations of the configured norm groups,
#' e.g. gender x age band) or overall. In every draw of the draw frame the
#' weighted percentile rank is computed; the median over draws is the norm
#' value, the quantiles `probs` the norm band. The standard score is
#' `mean + sd * qnorm(PR)` on the norm scale of the instrument (`norm_scale`,
#' default T = 50 + 10 z), with PR clamped to `[clamp, 1 - clamp]`. Its
#' columns carry the name of the norm scale; the documentation writes them
#' for the default, `T`.
#'
#' Two percentile conventions are available (`uncertainty$pctl_convention`):
#' `"midpoint"` (default), PR = P(X < x) + P(X = x) / 2, the classical
#' percentile rank, whose mean over the norm sample is exactly 50; and
#' `"cumulative"`, PR = P(X <= x), the relative cumulative frequency.
#'
#' @param fr Draw frame ([direct_frame()], [frequency_frame()] or the
#'   projection of the transfer module).
#' @param s Scale specification ([get_scale()]).
#' @param inst Instrument ([read_instrument()]).
#' @param cells Character vector of cell variables; `character(0)` gives the
#'   overall table.
#' @param probs Lower and upper quantile of the norm band.
#' @param clamp Percentile ranks are clamped to `[clamp, 1 - clamp]` before
#'   the transformation to the norm scale.
#' @param convention `"midpoint"` or `"cumulative"` (see Details).
#' @return `weighted_norms()`: tibble with the cell variables, `raw`,
#'   `pctl_med`, `pctl_lo`, `pctl_hi`, `T`, `T_lo`, `T_hi`, `n_eff` (Kish
#'   effective N of the respondents' weights, averaged over the draws) and
#'   `n`. `norms_all_cells()`: list with
#'   `cells` (all configured cells) and `overall`.
#' @section Lifecycle:
#' Stable since normpipe 1.0.0: the name, the arguments and the documented
#' parts of the result are kept; later versions only add to them.
#' @export
weighted_norms <- function(fr, s, inst, cells = cell_vars(inst),
                           probs = inst$uncertainty$probs,
                           clamp = inst$uncertainty$pctl_clamp,
                           convention = inst$uncertainty$pctl_convention %||% "midpoint") {
  n_bad <- sum(is.na(fr$score) | is.na(fr$weight))
  if (n_bad > 0) stop("weighted_norms(): ", n_bad, " rows with NA score or weight.")
  grid <- s$grid
  sc <- inst$norm_scale
  fr <- add_groups(fr, inst)
  fr <- fr[in_groups(fr, inst, cells), , drop = FALSE]
  if (!length(cells)) fr$.cell <- "overall"
  grp <- if (length(cells)) cells else ".cell"
  # effective N of a cell: every draw of a draw frame carries bootstrap
  # weights, so the weights of a person are averaged over the draws first
  neff <- fr %>% filter(.data$shape >= 1) %>%
    group_by(across(all_of(grp)), .data$id) %>%
    summarise(w = mean(.data$weight), .groups = "drop_last") %>%
    summarise(n_eff = sum(.data$w)^2 / sum(.data$w^2), n = n(), .groups = "drop")
  out <- fr %>%
    group_by(across(all_of(grp)), .data$.draw) %>%
    reframe(raw = grid, pctl = wecdf_at(.data$score, .data$weight, grid, convention)) %>%
    group_by(across(all_of(grp)), .data$raw) %>%
    summarise(pctl_med = median(.data$pctl),
              pctl_lo = quantile(.data$pctl, probs[1]),
              pctl_hi = quantile(.data$pctl, probs[2]),
              .groups = "drop") %>%
    mutate(T = pctl_to_norm(.data$pctl_med, clamp, sc),
           T_lo = pctl_to_norm(.data$pctl_lo, clamp, sc),
           T_hi = pctl_to_norm(.data$pctl_hi, clamp, sc)) %>%
    left_join(neff, by = grp)
  if (!length(cells)) out$.cell <- NULL
  norm_public(out, inst)
}

#' @rdname weighted_norms
#' @export
norms_all_cells <- function(fr, s, inst) {
  cv <- cell_vars(inst)
  res <- list()
  if (length(cv)) res$cells <- weighted_norms(fr, s, inst, cells = cv)
  if (isTRUE(inst$groups$overall) || !length(cv))
    res$overall <- weighted_norms(fr, s, inst, cells = character(0))
  res
}

#' Norms as a smooth function of a continuous variable (kernel norms)
#'
#' For every value of the continuous variable of the instrument
#' (`groups$continuous`, typically age in years), respondents are weighted
#' by an Epanechnikov kernel around that value whose half-width widens from
#' `half_width0` until the local effective sample size reaches `min_neff`,
#' up to `max_half_width`. This is done separately for every combination of
#' the categorical variables in `groups$continuous$by` (default: all
#' categorical group variables, e.g. per gender), or pooled when there is
#' none. Pseudo-observations of the prior carry their own values of the
#' variable, so they enter the kernel like everybody else.
#'
#' @inheritParams weighted_norms
#' @param values Target values of the continuous variable (default: the
#'   configured range in steps of one).
#' @param ages Name of `values` in normpipe < 0.4; still accepted.
#' @return Tibble with the `by` variables, the continuous variable,
#'   `half_width`, `n_eff_local`, `raw`, `pctl_med`, `pctl_lo`, `pctl_hi`,
#'   `T`, `T_lo`, `T_hi`.
#' @section Lifecycle:
#' Stable since normpipe 1.0.0: the name, the arguments and the documented
#' parts of the result are kept; later versions only add to them. The
#' argument `ages` is the name of `values` in normpipe < 0.4 and is kept as
#' an alias.
#' @export
weighted_norms_continuous <- function(fr, s, inst, values = inst$groups$continuous$values,
                                      probs = inst$uncertainty$probs,
                                      clamp = inst$uncertainty$pctl_clamp,
                                      convention = inst$uncertainty$pctl_convention %||% "midpoint",
                                      ages = NULL) {
  cont <- inst$groups$continuous
  if (is.null(cont)) stop("The instrument has no `groups$continuous` block.")
  if (!is.null(ages)) values <- ages
  k <- cont$kernel; cv <- cont$variable; sc <- inst$norm_scale
  grid <- s$grid
  fr <- fr[!is.na(fr[[cv]]), , drop = FALSE]
  by <- cont$by
  # one row per combination of the `by` groups, in the configured order
  combos <- if (length(by)) do.call(expand.grid, c(map(inst$groups$vars[rev(by)], "levels"),
                                                   stringsAsFactors = FALSE))[by] else data.frame(.all = NA)
  map_dfr(seq_len(nrow(combos)), \(ci) {
    pg <- fr
    for (b in by) pg <- pg[!is.na(pg[[b]]) & pg[[b]] == combos[[b]][ci], , drop = FALSE]
    draws <- sort(unique(pg$.draw))
    persons <- pg %>% group_by(.data$id) %>%
      summarise(.x = first(.data[[cv]]), w_bar = mean(.data$weight), .groups = "drop")
    map_dfr(values, \(a) {
      hw <- k$half_width0
      repeat {
        kern <- pmax(0, 1 - ((persons$.x - a) / hw)^2)
        wk <- persons$w_bar * kern
        neff <- if (sum(wk) > 0) sum(wk)^2 / sum(wk^2) else 0
        if (neff >= k$min_neff || hw >= k$max_half_width) break
        hw <- hw + 1
      }
      kv <- setNames(kern, persons$id)
      pa <- pg %>% mutate(w = .data$weight * kv[as.character(.data$id)]) %>%
        filter(.data$w > 0)
      P <- sapply(draws, \(d) { x <- pa[pa$.draw == d, ]; wecdf_at(x$score, x$w, grid, convention) })
      P <- matrix(P, nrow = length(grid))
      out <- tibble(.x = a, half_width = hw, n_eff_local = neff, raw = grid,
                    pctl_med = apply(P, 1, median),
                    pctl_lo = apply(P, 1, quantile, probs[1]),
                    pctl_hi = apply(P, 1, quantile, probs[2])) %>%
        mutate(T = pctl_to_norm(.data$pctl_med, clamp, sc),
               T_lo = pctl_to_norm(.data$pctl_lo, clamp, sc),
               T_hi = pctl_to_norm(.data$pctl_hi, clamp, sc))
      names(out)[1] <- cv
      if (length(by)) out <- bind_cols(as_tibble(combos[rep(ci, nrow(out)), by, drop = FALSE]), out)
      out
    })
  }) %>% norm_public(inst)
}

#' Kish effective sample size
#' @noRd
kish_neff <- function(w) sum(w)^2 / sum(w^2)
