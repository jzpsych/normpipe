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

#' Impute skipped items of the reference sample
#'
#' A respondent of the reference sample who skipped single items stays in the
#' sample (the inclusion rule is `min_answered` of each scale). With
#' `missing$reference: impute` in the instrument, [run_norming()] replaces each
#' skipped item by a draw from its distribution in the reference sample given
#' the answered items of its block, in `missing$n_imputations` completed data
#' sets. Every scale score, bridge predictor and observed part of a transfer
#' outcome is then computed from the completed items, and each draw of a draw
#' frame uses one of the completed data sets, so that the norm band carries
#' the imputation uncertainty. With `missing$reference: prorate` (default)
#' scores are computed from the answered items.
#'
#' The block of an item is its bridge predictor (`transfer$bridge`) when the
#' instrument has one that contains the item, otherwise the shortest scale
#' that contains it. The conditioning value is the rounded mean value of the
#' answered other items of the block. The conditional distributions are
#' estimated from the respondents who answered the item (counts + 0.5); when
#' no other item of the block is answered, or the conditioning value does not
#' occur among them, the marginal distribution is used.
#'
#' @param dat Data frame with item columns (raw response codes).
#' @param inst Instrument ([read_instrument()]).
#' @param M Number of completed data sets.
#' @param seed Random seed.
#' @return List of `M` data frames like `dat`, with the missing responses of
#'   the item columns filled in; attribute `n_imputed` gives the number of
#'   imputed responses per item column.
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' sim <- simulate_reference(inst, N = 100)
#' sim$items[1:5, 1] <- NA
#' imp <- impute_items(sim$items, inst, M = 3)
#' sapply(imp, function(d) d[1:5, 1])
#' @export
impute_items <- function(dat, inst, M = inst$missing$n_imputations %||% 20, seed = 1) {
  present <- inst$items$item[inst$items$col %in% names(dat)]
  idx <- item_indices(dat, inst, items = present)
  vals <- as.matrix(item_values(idx, inst, items = present))
  cols <- inst$items$col[match(present, inst$items$item)]
  colnames(vals) <- cols
  blocks <- imputation_blocks(inst, present)
  n_imp <- colSums(is.na(vals))
  out <- rep(list(dat), M)
  if (!any(n_imp > 0)) return(structure(out, n_imputed = n_imp))
  set.seed(seed)
  for (j in which(n_imp > 0)) {
    it <- present[j]
    cats <- inst$items$categories[[match(it, inst$items$item)]]
    partners <- setdiff(blocks[[as.character(it)]], it)
    pc <- cols[match(partners, present)]
    key <- if (length(pc)) round(rowMeans(vals[, pc, drop = FALSE], na.rm = TRUE)) else rep(NA_real_, nrow(vals))
    key[is.nan(key)] <- NA
    k <- idx[[cols[j]]] + 1L                         # category position, NA = skipped
    obs <- !is.na(k)
    marg <- tabulate(k[obs], length(cats)) + 0.5
    keys <- sort(unique(key[obs & !is.na(key)]))
    cond <- map(keys, \(q) tabulate(k[obs & !is.na(key) & key == q], length(cats)) + 0.5)
    na <- which(!obs)
    draws <- matrix(NA_integer_, length(na), M)
    for (r in seq_along(na)) {
      q <- match(key[na[r]], keys)
      pr <- if (is.na(q)) marg else cond[[q]]
      draws[r, ] <- sample.int(length(cats), M, replace = TRUE, prob = pr)
    }
    for (m in seq_len(M)) out[[m]][[cols[j]]][na] <- cats[draws[, m]]
  }
  structure(out, n_imputed = n_imp)
}

#' Block of items used to impute each item (bridge predictor or shortest scale)
#' @noRd
imputation_blocks <- function(inst, present) {
  br <- if (isTRUE(inst$transfer$enabled) || !is.null(inst$transfer$bridge)) map(inst$transfer$bridge, "items") else list()
  scales <- unlist(map(inst$versions, \(v) map(v$scales, "items")), recursive = FALSE)
  scales <- scales[map_lgl(scales, \(x) all(x %in% present))]
  scales <- scales[order(lengths(scales))]
  set_names(map(present, \(it) {
    b <- detect(br, \(x) it %in% x && all(x %in% present))
    if (is.null(b)) b <- detect(scales, \(x) it %in% x && length(x) > 1)
    as.integer(b %||% it)
  }), as.character(present))
}
