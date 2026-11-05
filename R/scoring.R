# ==============================================================================
# scoring.R -- from item responses to scale scores
#
# Item responses enter as the codes of the data file (e.g., 1 to 4). They are
# turned into (a) CATEGORY INDICES 0, ..., K - 1 for the item models and
# (b) VALUES in integer units for the raw scores, according to the item table
# of the instrument. Reverse keying is part of the value mapping, so the
# category indices always run in the order of the data codes and the item
# models see the same ordering in every sample.
# ==============================================================================

#' Recode item responses to category indices
#'
#' @param dat Data frame with the item columns, named as in the item table
#'   (`<prefix><item>`).
#' @param inst Instrument ([read_instrument()]).
#' @param items Item ids to recode; default all items of the instrument.
#' @return Tibble of integer category indices (0-based, in the order of the
#'   configured `categories`), `NA` preserved, columns named as in the item
#'   table. Codes outside the configured categories are an error.
#' @export
item_indices <- function(dat, inst, items = inst$items$item) {
  tab <- inst$items[match(items, inst$items$item), ]
  missing_cols <- setdiff(tab$col, names(dat))
  if (length(missing_cols))
    stop("Item columns missing in the data: ", paste(head(missing_cols, 5), collapse = ", "),
         if (length(missing_cols) > 5) " ...")
  out <- map2(tab$col, tab$categories, \(cn, cats) {
    x <- as.numeric(dat[[cn]])
    k <- match(x, cats) - 1L
    bad <- !is.na(x) & is.na(k)
    if (any(bad))
      stop("Item ", cn, ": ", sum(bad), " response(s) outside the configured ",
           "categories (", paste(cats, collapse = ","), "): e.g. ",
           paste(head(unique(x[bad]), 5), collapse = ","))
    k
  })
  names(out) <- tab$col
  as_tibble(out)
}

#' Item values from category indices
#'
#' @param idx Category indices as returned by [item_indices()].
#' @param inst Instrument ([read_instrument()]).
#' @param items Item ids; default all items of the instrument.
#' @return Tibble of item values in integer units (multiply by
#'   `attr(inst$items, "unit")` for the configured values), reverse keying
#'   applied.
#' @export
item_values <- function(idx, inst, items = inst$items$item) {
  tab <- inst$items[match(items, inst$items$item), ]
  out <- map2(tab$col, tab$values, \(cn, v) v[idx[[cn]] + 1L])
  names(out) <- tab$col
  as_tibble(out)
}

#' Score every scale of every version
#'
#' @param dat Data frame with the item columns (raw response codes).
#' @param inst Instrument ([read_instrument()]).
#' @param prorate If `TRUE`, a scale is scored from the available items when
#'   at least the proportion `min_answered` (scale key, default 1) of its
#'   items is answered (item means, or sums scaled up to the full length);
#'   otherwise any missing item gives `NA`.
#' @param partial If `TRUE`, item columns absent from `dat` are allowed (e.g.,
#'   a reference sample that answered only a subset of the items); scales that
#'   need an absent item are returned as `NA`.
#' @return Tibble with one column per outcome (`<version>_<scale>`) on the raw
#'   metric of the scale.
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' sim <- simulate_reference(inst, N = 50)
#' head(score_scales(sim$items, inst))
#' @export
score_scales <- function(dat, inst, prorate = TRUE, partial = FALSE) {
  present <- if (partial) inst$items$item[inst$items$col %in% names(dat)] else inst$items$item
  idx <- item_indices(dat, inst, items = present)
  vals <- item_values(idx, inst, items = present)
  out <- list()
  for (v in inst$versions) for (s in v$scales) {
    out[[paste0(v$id, "_", s$id)]] <- if (all(s$items %in% present))
      score_one_scale(vals, s, inst, prorate) else rep(NA_real_, nrow(dat))
  }
  as_tibble(out)
}

score_one_scale <- function(vals, s, inst, prorate = TRUE) {
  unit <- attr(inst$items, "unit")
  col <- function(ids) inst$items$col[match(ids, inst$items$item)]
  mean_block <- function(ids) {
    x <- as.matrix(vals[col(ids)])
    n_ok <- rowSums(!is.na(x))
    m <- rowSums(x, na.rm = TRUE) / n_ok * unit
    ok <- if (prorate) n_ok >= s$min_answered * length(ids) & n_ok > 0 else
      n_ok == length(ids)
    m[!ok] <- NA
    m
  }
  switch(s$metric,
    mean = mean_block(s$items),
    sum  = mean_block(s$items) * s$n_items,
    mean_of_means = {
      m <- map(s$components, mean_block)
      Reduce(`+`, m) / length(m)
    })
}

#' Snap raw scores to the attainable grid of a scale
#'
#' Prorated scores and rounded published values lie between grid points; the
#' nearest attainable score is used (ties go up). Values outside the range are
#' clamped.
#'
#' @param x Numeric vector of raw scores.
#' @param s Scale specification ([get_scale()]).
#' @return Numeric vector of grid values (`NA` preserved).
#' @export
snap_to_grid <- function(x, s) {
  g <- s$grid
  i <- findInterval(x, g)
  i <- pmin(pmax(i, 1), length(g) - 1)
  lo <- g[i]; hi <- g[i + 1]
  res <- ifelse(x - lo < hi - x, lo, hi)
  res[x <= g[1]] <- g[1]; res[x >= g[length(g)]] <- g[length(g)]
  res[is.na(x)] <- NA
  res
}

#' Summed score (integer units) of a scale from item indices, for the EAP
#' tables (complete cases only)
#' @noRd
summed_values <- function(idx, s, inst) {
  vals <- item_values(idx, inst, s$items)
  rowSums(as.matrix(vals))
}

#' Raw score of a scale from its summed value (integer units)
#' @noRd
sum_to_raw <- function(S, s) {
  if (s$metric == "sum") S * s$unit else S * s$unit / s$n_items
}
