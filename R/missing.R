# ==============================================================================
# missing.R -- one handling of skipped items
#
# Who is in a sample is decided by the share of the presented items a
# respondent answered (missing$min_answered); people below it are unit
# non-respondents and are left to the weights. Single skipped items of the
# others are imputed in M completed data sets by a nearest-neighbour hot deck
# (impute_items()), and everything that needs scores uses the completed data:
# directly normed scales, bridge predictors, the observed part of transfer
# outcomes, and the donor samples of the transfer models. Item models are
# calibrated on the responses as given (full-information likelihood).
# ==============================================================================

#' Impute skipped items
#'
#' With `missing: {method: impute}` in the instrument, [run_norming()] keeps
#' every respondent who answered at least the share `missing$min_answered` of
#' the items presented to the sample and replaces each skipped item in
#' `missing$n_imputations` completed data sets. All scores are then computed
#' from the completed items, and each draw of a draw frame uses one of the
#' completed data sets, so that the norm bands carry the imputation
#' uncertainty. With `method: prorate` (default) scores are means of the
#' answered items, subject to `min_answered` of each scale.
#'
#' `impute_items()` is a nearest-neighbour hot deck. For a respondent with
#' skipped items in a block, the candidates are the respondents who answered
#' the whole block; their distance is the mean absolute difference of the
#' item values (scaled to the range of each item) over the block items the
#' respondent answered, plus optional covariates. One of the `k` nearest
#' candidates is drawn per completed data set, and all skipped items of the
#' block are copied from that candidate, which keeps the association between
#' the items of a block. The block of an item is the largest scale below the
#' total that contains it and whose items were all presented (for a
#' questionnaire with facets and domains: the domain). A respondent who
#' answered no other item of the block is matched on all items they answered,
#' and one who answered nothing receives draws from the marginal distribution.
#'
#' `answered_share()` gives the share of the presented items each respondent
#' answered.
#'
#' @param dat Data frame with item columns (raw response codes).
#' @param inst Instrument ([read_instrument()]).
#' @param M Number of completed data sets.
#' @param seed Random seed.
#' @param covariates Optional data frame with one row per row of `dat`
#'   (e.g. age and gender); numeric columns enter the distance scaled to their
#'   range, other columns as match or mismatch, each with the weight of one
#'   item. Respondents with a missing covariate are matched without it.
#' @param k Number of nearest candidates to draw from.
#' @return `impute_items()`: list of `M` data frames like `dat` with the
#'   missing responses of the item columns filled in; attribute `n_imputed`
#'   gives the number of imputed responses per item column.
#'   `answered_share()`: numeric vector.
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' sim <- simulate_reference(inst, N = 100)
#' sim$items[1:5, 1] <- NA
#' imp <- impute_items(sim$items, inst, M = 3)
#' sapply(imp, function(d) d[1:5, 1])
#' answered_share(sim$items, inst)[1:6]
#' @export
impute_items <- function(dat, inst, M = inst$missing$n_imputations %||% 20, seed = 1,
                         covariates = NULL, k = inst$missing$neighbours %||% 10) {
  present <- inst$items$item[inst$items$col %in% names(dat)]
  cols <- inst$items$col[match(present, inst$items$item)]
  idx <- item_indices(dat, inst, items = present)
  vals <- as.matrix(item_values(idx, inst, items = present)[cols])
  # item values on 0..1, so that every item has the same weight in the distance
  rng <- map(inst$items$values[match(present, inst$items$item)], range)
  Z <- sweep(sweep(vals, 2, map_dbl(rng, 1)), 2, pmax(map_dbl(rng, diff), 1e-12), `/`)
  colnames(Z) <- cols
  miss <- is.na(Z)
  n_imp <- colSums(miss); names(n_imp) <- cols
  out <- rep(list(dat), M)
  if (!any(miss)) return(structure(out, n_imputed = n_imp))
  C <- covariate_matrix(covariates, nrow(dat))
  block_of <- imputation_blocks(inst, present)
  blocks <- unique(block_of)
  codes <- as.matrix(as.data.frame(dat)[cols])            # raw codes to copy
  set.seed(seed)
  for (b in blocks) {
    bc <- cols[match(b, present)]
    own <- cols[map_lgl(block_of, \(x) identical(x, b))]   # items imputed within this block
    todo <- which(rowSums(miss[, own, drop = FALSE]) > 0)
    if (!length(todo)) next
    cand <- which(rowSums(miss[, bc, drop = FALSE]) == 0)
    for (i in todo) {
      mi <- own[miss[i, own]]
      use <- bc[!miss[i, bc]]                              # answered items of the block
      pool <- cand
      if (!length(use)) {                                  # nothing of the block answered
        use <- cols[!miss[i, ]]
        pool <- if (length(use)) which(rowSums(miss[, c(use, mi), drop = FALSE]) == 0) else integer(0)
      }
      if (length(pool) < 1) {
        # no candidate with the needed items: item by item, from everyone who answered it
        for (j in mi) {
          src <- which(!miss[, j])
          pick <- src[sample.int(length(src), M, replace = TRUE)]
          for (m in seq_len(M)) out[[m]][[j]][i] <- codes[pick[m], j]
        }
        next
      }
      d <- if (length(use)) rowMeans(abs(sweep(Z[pool, use, drop = FALSE], 2, Z[i, use]))) * length(use) else rep(0, length(pool))
      w <- length(use)
      if (!is.null(C)) {
        ok <- !is.na(C[i, ])
        if (any(ok)) {
          dc <- abs(sweep(C[pool, ok, drop = FALSE], 2, C[i, ok]))
          dc[is.na(dc)] <- 1
          d <- d + rowSums(dc); w <- w + sum(ok)
        }
      }
      d <- d / max(w, 1)
      near <- pool[order(d, stats::runif(length(pool)))[seq_len(min(k, length(pool)))]]
      pick <- near[sample.int(length(near), M, replace = TRUE)]
      for (j in mi) for (m in seq_len(M)) out[[m]][[j]][i] <- codes[pick[m], j]
    }
  }
  structure(out, n_imputed = n_imp)
}

#' @rdname impute_items
#' @export
answered_share <- function(dat, inst) {
  cols <- intersect(inst$items$col, names(dat))
  if (!length(cols)) return(rep(NA_real_, nrow(dat)))
  rowMeans(!is.na(as.data.frame(dat)[cols]))
}

#' Covariates of the hot deck on 0..1 (numeric: range; other: dummy columns)
#' @noRd
covariate_matrix <- function(covariates, n) {
  if (is.null(covariates) || !ncol(covariates)) return(NULL)
  stopifnot(nrow(covariates) == n)
  do.call(cbind, map(as.list(as.data.frame(covariates)), \(x) {
    if (is.numeric(x)) {
      r <- range(x, na.rm = TRUE)
      matrix((x - r[1]) / max(diff(r), 1e-12), ncol = 1)
    } else {
      lv <- sort(unique(stats::na.omit(as.character(x))))
      m <- outer(as.character(x), lv, `==`) / 2            # a mismatch counts 1 in total
      m[is.na(x), ] <- NA
      m
    }
  }))
}

#' Block of items used to impute each item: the largest scale below the total
#' whose items were all presented (all presented items when there is none)
#' @noRd
imputation_blocks <- function(inst, present) {
  sc <- unlist(map(inst$versions, \(v) v$scales), recursive = FALSE)
  sc <- sc[map_lgl(sc, \(s) !identical(s$level, "total") && length(s$items) > 1 && all(s$items %in% present))]
  items <- map(sc, \(s) as.integer(s$items))
  items <- items[order(-lengths(items))]
  set_names(map(present, \(it) {
    b <- detect(items, \(x) it %in% x)
    as.integer(b %||% present)
  }), as.character(present))
}

#' Settings of the `missing` block (accepts the 0.2.6 key `reference`)
#' @noRd
missing_cfg <- function(inst) {
  m <- inst$missing %||% list()
  method <- m$method %||% m$reference %||% "prorate"
  if (!method %in% c("prorate", "impute")) stop("missing$method must be 'prorate' or 'impute'.")
  list(method = method, n_imputations = as.integer(m$n_imputations %||% 20),
       min_answered = m$min_answered %||% 0.5, neighbours = as.integer(m$neighbours %||% 10),
       covariates = as.character(unlist(m$covariates)))
}
