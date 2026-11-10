# ==============================================================================
# transfer-observed.R -- observed (bridge) items of a transfer outcome and the
# count structure of the beta-binomial
#
# The bridge items that belong to an outcome are answered in the reference
# sample. The transfer model therefore predicts only the part of the score that
# the other items contribute (the "rest"), and the projection adds the
# respondent's observed part. This holds for every family; it is not a
# configuration option.
# ==============================================================================

#' Observed items of a transfer outcome
#'
#' Items of the outcome that belong to a bridge block in use. Items of
#' excluded blocks (e.g. in the leave-own-domain-out cross-check) count as
#' unobserved.
#' @noRd
observed_items <- function(s, inst, exclude = NULL) {
  br <- inst$transfer$bridge
  used <- br[!map_chr(br, "id") %in% exclude]
  intersect(s$items, unique(unlist(map(used, "items"))))
}

#' Weight of each item value (integer units) in the raw score of a scale
#'
#' The raw score is linear in the item values: sum (weight = unit), mean
#' (unit / n) or mean of component means (unit / (n_component * n_components)).
#' @noRd
item_weights <- function(s, inst) {
  unit <- attr(inst$items, "unit")
  w <- switch(s$metric,
    sum = setNames(rep(unit, length(s$items)), s$items),
    mean = setNames(rep(unit / length(s$items), length(s$items)), s$items),
    mean_of_means = {
      C <- length(s$components)
      unlist(unname(map(s$components, \(it) setNames(rep(unit / (length(it) * C), length(it)), it))))
    },
    stop("Scale metric '", s$metric, "' is not linear in the item values; ",
         "observed bridge items cannot be separated from the rest."))
  w[as.character(s$items)]
}

#' Everything the projection needs about the observed items (stored with the fit)
#'
#' Weights, the range of the observed part (integer units, for the
#' beta-binomial) and of the rest (raw metric, for the censored families),
#' and, per observed item, the distribution of its value given the mean of
#' the other items of its bridge block in the donor data. The latter is a
#' fallback for a reference respondent whose item is still missing in the
#' prediction frame (without central imputation): a new value is drawn in
#' every posterior draw, so that the norm band carries this uncertainty.
#' @noRd
observed_info <- function(s, inst, known, dat = NULL) {
  w <- item_weights(s, inst)
  v <- setNames(inst$items$values[match(s$items, inst$items$item)], s$items)
  vmin <- map_dbl(v, min); vmax <- map_dbl(v, max)
  kc <- as.character(known); uc <- setdiff(names(w), kc)
  out <- list(known_items = known, w = w[kc],
              range_units = c(min = sum(vmin[kc]), max = sum(vmax[kc])),
              rest = c(min = sum(w[uc] * vmin[uc]), max = sum(w[uc] * vmax[uc])),
              cond = list())
  if (!length(known)) return(out)
  br <- inst$transfer$bridge
  for (it in known) {
    b <- br[[which(map_lgl(br, \(x) it %in% x$items))[1]]]
    partners <- setdiff(b$items, it)
    vals <- sort(unique(v[[as.character(it)]]))
    cnt <- NULL; keys <- numeric(0)
    if (!is.null(dat) && length(partners)) {
      pk <- round(rowMeans(as.matrix(dat[paste0("bi_", partners)]), na.rm = TRUE))
      x <- dat[[paste0("bi_", it)]]
      ok <- !is.na(pk) & !is.na(x)
      keys <- sort(unique(pk[ok]))
      cnt <- vapply(keys, \(k) tabulate(match(x[ok][pk[ok] == k], vals), length(vals)), numeric(length(vals)))
      cnt <- matrix(cnt, nrow = length(vals))          # values x keys
    }
    marg <- if (!is.null(cnt)) rowSums(cnt) else rep(1, length(vals))
    out$cond[[as.character(it)]] <- list(partners = partners, values = vals, keys = keys,
                                         prob = if (!is.null(cnt)) sweep(cnt + 0.5, 2, colSums(cnt + 0.5), `/`),
                                         marginal = (marg + 0.5) / sum(marg + 0.5))
  }
  out
}

#' Observed part of the score per draw (draws x persons)
#'
#' `wt`: weight per observed item (raw metric for the censored families,
#' 1 for the integer-unit sums of the beta-binomial). With `D = 1` and
#' complete data (fitting) the result is a vector.
#' @noRd
observed_part <- function(d, info, wt = info$w, D = 1) {
  known <- info$known_items
  if (!length(known)) return(if (D == 1) rep(0, nrow(d)) else matrix(0, D, nrow(d)))
  cols <- paste0("bi_", known)
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop("Bridge item column(s) ", paste(miss, collapse = ", "), " missing; ",
                         "build the frame with bridge_scores().")
  n <- nrow(d)
  tot <- matrix(0, D, n)
  for (j in seq_along(known)) {
    x <- d[[cols[j]]]
    vj <- matrix(x, D, n, byrow = TRUE)
    na <- which(is.na(x))
    if (length(na) && D > 1) {
      cj <- info$cond[[as.character(known[j])]]
      pk <- if (length(cj$partners))
        round(rowMeans(as.matrix(d[paste0("bi_", cj$partners)])[na, , drop = FALSE], na.rm = TRUE)) else
        rep(NA, length(na))
      for (k in seq_along(na)) {
        col <- match(pk[k], cj$keys)
        pr <- if (!is.na(col)) cj$prob[, col] else cj$marginal
        vj[, na[k]] <- cj$values[sample.int(length(cj$values), D, TRUE, prob = pr)]
      }
    }
    tot <- tot + wt[[j]] * vj
  }
  if (D == 1) drop(tot) else tot
}

#' Summed-score representation of a scale for the beta-binomial family
#'
#' The attainable raw scores must form an arithmetic sequence: the summed
#' score in integer units is `raw / unit_per_sum`, the trials of the
#' beta-binomial are the range of the attainable sums, and `offset` is the
#' smallest attainable sum (non-zero when item values do not start at 0).
#' With observed items, their value range is removed from the modelled sum.
#'
#' A dense scale (mean of facet means of unequal length, tabulated on a
#' reporting grid) has no count structure: modelled as a count out of the
#' grid steps, the beta-binomial loses the floor and widens the upper tail
#' (PID-5 domain pilot, 2026-11). Such scales keep the censored families.
#' @noRd
bb_scale <- function(s, known = NULL, values = NULL) {
  if (isTRUE(s$dense))
    stop("Scale '", s$version, "_", s$id, "': a dense scale (mean of means of unequal ",
         "length, tabulated on the reporting grid) has no count structure; use ",
         "family 'gaussian' or 'skew_normal' for it (transfer$outcomes).")
  if (!isTRUE(all.equal(diff(s$grid), rep(s$step, length(s$grid) - 1))))
    stop("Scale '", s$version, "_", s$id, "': the beta-binomial family needs an ",
         "arithmetic sequence of attainable scores (sum, mean, or mean_of_means ",
         "with equal component lengths).")
  lo <- min(s$sums); hi <- max(s$sums)
  if (!isTRUE(all.equal(sort(unique(s$sums)), lo:hi)))
    stop("Scale '", s$version, "_", s$id, "': the attainable sums have gaps; ",
         "the beta-binomial family is not available.")
  if (length(known)) {
    lo <- lo - values[["min"]]; hi <- hi - values[["max"]]
    if (hi <= lo) stop("Scale '", s$version, "_", s$id, "': all items are bridge items; ",
                       "the scale should be observed directly.")
  }
  list(unit_per_sum = s$step, offset = as.integer(round(lo)), trials = as.integer(round(hi - lo)))
}

#' Can a scale be modelled by the beta-binomial?
#' @noRd
bb_possible <- function(s) !inherits(try(bb_scale(s), silent = TRUE), "try-error")
