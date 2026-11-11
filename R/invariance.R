# ==============================================================================
# invariance.R -- measurement invariance of ordered items across groups
#
# Multiple-group CFA for ordered items (lavaan, WLSMV, theta parameterization)
# with the model syntax written out by hand, identified as in Wu and Estabrook
# (2016) and with item intercepts fixed to zero in every group and model (the
# convention of Mplus, which has no item intercepts for categorical items).
#
#   configural    thresholds and loadings free; residual variances 1; all
#                 latent means and intercepts 0
#   scalar        thresholds and item loadings equal; residual variances 1 in
#                 group 1 and free in the other groups; means of the
#                 first-order factors 0 in group 1 and free in the other
#                 groups; second-order loadings free
#   strict        scalar, and residual variances 1 in every group
#   second_order  scalar, and second-order loadings equal, intercepts of the
#                 first-order factors 0 in every group, means of the
#                 second-order factors free in groups 2+ (group differences in
#                 the first-order factors run through the second-order
#                 factors only); supplementary, not part of the sequence
#
# "Scalar" is therefore invariance of the item parameters. The syntax is not
# generated with semTools::measEq.syntax(), which leaves latent means of a
# higher-order model free that are not identified.
#
# Noninvariance is flagged as in Chen (2007, p. 501, N > 300): a drop of the
# CFI by .010 or more, supplemented by a rise of the RMSEA by .015 or more or
# of the SRMR by .010 or more (Chen gives .030 for the SRMR when loadings are
# tested; the stricter value is used at every step).
# ==============================================================================

#' Factor structure for the invariance models
#'
#' Collects the first-order factors (blocks of items) and, where the
#' instrument defines scales from such blocks, the second-order factors above
#' them.
#'
#' @param inst Instrument ([read_instrument()]).
#' @param scales Outcome names (`<version>_<scale>`). A scale with components
#'   becomes a second-order factor over one first-order factor per component,
#'   any other scale a first-order factor of its items. Default (`NULL`): the
#'   bridge blocks of the transfer module (`transfer$bridge`) as first-order
#'   factors, and every scale of the instrument that is made of at least
#'   three of these blocks as a second-order factor; blocks that belong to no
#'   such scale stay first-order factors.
#' @return List with `first` (named list of item ids per first-order factor)
#'   and `second` (named list of first-order factor names per second-order
#'   factor, possibly empty).
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' invariance_structure(inst, c("full_a", "full_b"))
#' @export
invariance_structure <- function(inst, scales = NULL) {
  first <- list(); second <- list()
  comp_of <- function(s) map(s$components, as.integer)   # named list of item ids
  if (is.null(scales)) {
    br <- inst$transfer$bridge
    if (is.null(br))
      stop("The instrument has no transfer$bridge blocks: name the scales to test in `scales`.")
    first <- set_names(map(br, \(b) as.integer(b$items)), map_chr(br, "id"))
    for (v in inst$versions) for (s in v$scales) {
      if (is.null(s$components) || length(s$components) < 3) next
      cp <- comp_of(s)
      same <- all(names(cp) %in% names(first)) &&
        all(map_lgl(names(cp), \(n) setequal(cp[[n]], first[[n]])))
      taken <- any(map_lgl(second, \(x) setequal(x, names(cp))))
      if (same && !taken && !any(names(cp) %in% unlist(second))) second[[s$id]] <- names(cp)
    }
  } else {
    for (o in scales) {
      s <- get_scale(inst, o)
      if (!is.null(s$components) && !identical(s$model, "composite")) {
        cp <- comp_of(s)
        first[names(cp)] <- cp
        second[[s$id]] <- names(cp)
      } else first[[s$id]] <- as.integer(s$items)
    }
  }
  if (anyDuplicated(unlist(first)))
    stop("An item belongs to more than one first-order factor: ",
         paste(unique(unlist(first)[duplicated(unlist(first))]), collapse = ", "))
  list(first = first, second = second)
}

#' lavaan syntax of one invariance model
#'
#' @param structure Factor structure ([invariance_structure()]) with the
#'   items given as column names (as [invariance_sequence()] passes it).
#' @param nthr Named integer vector: number of thresholds of every item.
#' @param n_groups Number of groups.
#' @param level `"configural"`, `"scalar"`, `"strict"` or `"second_order"`
#'   (see Details).
#' @param form `"higher_order"` (default) or `"domain_residual"`: the
#'   second-order factors load on the items directly and the items of every
#'   first-order factor below them share a residual covariance. The two forms
#'   are equivalent for the configural model when these first-order factors
#'   have two items each; the second form has no boundary at a factor
#'   variance of zero and is available for the configural model only.
#' @details Thresholds and loadings are free in the configural model
#'   (residual variances 1, latent means and intercepts 0). The scalar model
#'   holds thresholds and item loadings equal, frees the residual variances
#'   and the means of the first-order factors in groups 2+ and leaves the
#'   second-order loadings free. The strict model adds residual variances of 1
#'   in all groups. `"second_order"` adds equal second-order loadings to the
#'   scalar model and moves the free means to the second-order factors. Item
#'   intercepts are 0 throughout. Items with two categories keep a residual
#'   variance of 1 in every group.
#' @return The model syntax (one string).
#' @export
invariance_syntax <- function(structure, nthr, n_groups,
                              level = c("configural", "scalar", "strict", "second_order"),
                              form = c("higher_order", "domain_residual")) {
  level <- match.arg(level); form <- match.arg(form)
  ng <- n_groups; inv <- level != "configural"; strict <- level == "strict"
  rp <- function(x) paste0("c(", paste(rep(x, ng), collapse = ", "), ")")
  g1 <- function(first, rest) paste0("c(", paste(c(first, rep(rest, ng - 1)), collapse = ", "), ")")
  first <- structure$first; second <- structure$second
  below <- unlist(second, use.names = FALSE)
  item_lines <- function(v) c(
    sprintf("%s | %s", v, paste(map_chr(seq_len(nthr[[v]]), \(k)
      sprintf("%s*t%d", if (inv) rp(sprintf("th_%s_%d", v, k)) else rp("NA"), k)),
      collapse = " + ")),
    sprintf("%s ~ 0*1", v),
    sprintf("%s ~~ %s*%s", v,
            if (inv && !strict && nthr[[v]] > 1) g1(1, "NA") else rp(1), v))
  out <- character(0)
  if (form == "domain_residual") {
    if (level != "configural")
      stop("The domain-residual form is available for the configural model only.")
    if (any(lengths(first[below]) != 2))
      stop("The domain-residual form needs two items per first-order factor below a second-order factor.")
    for (d in names(second)) {
      it <- unlist(first[second[[d]]], use.names = FALSE)
      out <- c(out, sprintf("D_%s =~ 1*%s + %s", d, it[1],
                            paste(sprintf("%s*%s", rp("NA"), it[-1]), collapse = " + ")),
               sprintf("D_%s ~ %s*1", d, rp(0)),
               map_chr(second[[d]], \(f) sprintf("%s ~~ %s", first[[f]][1], first[[f]][2])))
    }
    for (f in setdiff(names(first), below)) {
      it <- first[[f]]
      out <- c(out, sprintf("F_%s =~ 1*%s + %s", f, it[1],
                            paste(sprintf("%s*%s", rp("NA"), it[-1]), collapse = " + ")),
               sprintf("F_%s ~ %s*1", f, rp(0)))
    }
    return(paste(c(out, unlist(map(unlist(first, use.names = FALSE), item_lines))),
                 collapse = "\n"))
  }
  for (f in names(first)) {
    it <- first[[f]]
    out <- c(out, sprintf("F_%s =~ 1*%s%s", f, it[1], if (length(it) > 1) paste0(" + ",
      paste(sprintf("%s*%s", if (inv) map_chr(it[-1], \(v) rp(paste0("l_", v))) else rp("NA"),
                    it[-1]), collapse = " + ")) else ""))
    out <- c(out, unlist(map(it, item_lines)))
    # mean of the first-order factor: free in groups 2+ under item-level
    # invariance; factors below no second-order factor keep it in "second_order"
    free_mean <- level %in% c("scalar", "strict") ||
      (level == "second_order" && !f %in% below)
    out <- c(out, sprintf("F_%s ~ %s*1", f, if (free_mean) g1(0, "NA") else rp(0)))
  }
  for (d in names(second)) {
    fs <- paste0("F_", second[[d]])
    out <- c(out,
      sprintf("D_%s =~ 1*%s + %s", d, fs[1], paste(map_chr(fs[-1], \(x)
        sprintf("%s*%s", if (level == "second_order") rp(paste0("g_", x)) else rp("NA"), x)),
        collapse = " + ")),
      sprintf("D_%s ~ %s*1", d, if (level == "second_order") g1(0, "NA") else rp(0)))
  }
  paste(out, collapse = "\n")
}

INVARIANCE_FIT <- c("cfi.scaled", "rmsea.scaled", "srmr", "chisq.scaled", "df.scaled")

# One model: fit indices and diagnostics, no parameter estimates.
invariance_fit <- function(syntax, d, items, missing) {
  warn <- character(0)
  f <- withCallingHandlers(
    tryCatch(lavaan::cfa(syntax, data = d, group = ".grp", ordered = items,
                         estimator = "WLSMV", missing = missing,
                         parameterization = "theta"),
             error = \(e) { warn <<- c(warn, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = \(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  na <- set_names(rep(NA_real_, length(INVARIANCE_FIT)), INVARIANCE_FIT)
  v <- if (is.null(f)) na else tryCatch(lavaan::fitMeasures(f, INVARIANCE_FIT), error = \(e) na)
  neg <- c(NA_integer_, NA_integer_)
  if (!is.null(f)) {
    pe <- lavaan::parameterEstimates(f)
    var <- pe$op == "~~" & pe$lhs == pe$rhs & pe$est < 0
    neg <- c(sum(var & grepl("^[FD]_", pe$lhs)), sum(var & pe$lhs %in% items))
  }
  as_tibble_row(set_names(as.numeric(v), INVARIANCE_FIT)) %>%
    mutate(n = if (is.null(f)) NA_integer_ else as.integer(sum(lavaan::lavInspect(f, "nobs"))),
           converged = if (is.null(f)) NA else lavaan::lavInspect(f, "converged"),
           neg_lv_variances = neg[1], neg_item_variances = neg[2],
           warnings = paste(unique(gsub("\\s+", " ", warn)), collapse = " | "))
}

# Merge response categories for the invariance test: the top category goes
# into the one below when it has fewer than min_n responses in any group
# (once), and a top or bottom category that is empty in a group is merged
# with its neighbour until none is (an empty category stops the estimation).
collapse_sparse <- function(x, g, min_n) {
  top_down <- function(x) { lv <- sort(unique(stats::na.omit(x))); x[which(x == lv[length(lv)])] <- lv[length(lv) - 1]; x }
  cells <- function(x) { lv <- sort(unique(stats::na.omit(x))); table(g, factor(x, levels = lv)) }
  tab <- cells(x)
  if (min_n > 0 && ncol(tab) > 2 && any(tab[, ncol(tab)] < min_n)) x <- top_down(x)
  repeat {
    tab <- cells(x)
    if (ncol(tab) <= 2) break
    if (any(tab[, ncol(tab)] == 0)) x <- top_down(x)
    else if (any(tab[, 1] == 0)) { lv <- as.numeric(colnames(tab)); x[which(x == lv[1])] <- lv[2] }
    else break
  }
  x
}

#' Changes in fit between nested invariance models and the Chen (2007) rule
#'
#' @param tab Fit table with the columns `model`, `compared_with`,
#'   `cfi.scaled`, `rmsea.scaled`, `srmr`, one row per model of one sequence.
#' @param cfi,rmsea,srmr Cutoffs: a model is flagged when the CFI falls by at
#'   least `cfi` and the RMSEA rises by at least `rmsea` or the SRMR by at
#'   least `srmr`, relative to the model in `compared_with`.
#' @return `tab` with `d_cfi`, `d_rmsea`, `d_srmr` and `flag` (`NA` for the
#'   first model and when an index is not available).
#' @references Chen, F. F. (2007). Sensitivity of goodness of fit indexes to
#'   lack of measurement invariance. *Structural Equation Modeling, 14*(3),
#'   464-504.
#' @export
invariance_flags <- function(tab, cfi = 0.010, rmsea = 0.015, srmr = 0.010) {
  ref <- match(tab$compared_with, tab$model)
  tab$d_cfi <- tab$cfi.scaled - tab$cfi.scaled[ref]
  tab$d_rmsea <- tab$rmsea.scaled - tab$rmsea.scaled[ref]
  tab$d_srmr <- tab$srmr - tab$srmr[ref]
  eps <- 1e-9                                  # cutoffs are inclusive
  tab$flag <- tab$d_cfi <= -cfi + eps &
    (tab$d_rmsea >= rmsea - eps | tab$d_srmr >= srmr - eps)
  tab
}

#' Sequence of invariance models across groups
#'
#' Fits the configural, scalar and strict model (and, with second-order
#' factors, the supplementary model with invariant second-order loadings) to
#' the items of a factor structure and compares each with its predecessor.
#' Only fit indices are returned.
#'
#' @param items Data frame of item responses (raw codes, columns named as in
#'   the item table).
#' @param group Grouping: a vector with one entry per row of `items`. Rows
#'   with a missing group are left out.
#' @param inst Instrument ([read_instrument()]).
#' @param structure Factor structure; default [invariance_structure()] of the
#'   instrument.
#' @param levels Models to fit, in this order. `"second_order"` is dropped
#'   when the structure has no second-order factor.
#' @param missing `"pairwise"` (default) or `"listwise"`.
#' @param min_cat_n When the highest response category of an item has fewer
#'   than `min_cat_n` responses in any group, it is merged with the category
#'   below in all groups, for the test only (`0`: never). Independently of
#'   this, a highest or lowest category that is empty in a group is merged
#'   with its neighbour, because an empty category stops the estimation.
#' @param label Value of the column `grouping` of the result.
#' @details When the robust test statistic of the configural higher-order
#'   model cannot be computed (typically a first-order factor variance at the
#'   boundary in a small group) and every first-order factor below a
#'   second-order factor has two items, the configural model is refitted in
#'   the equivalent domain-residual form (see [invariance_syntax()]); the
#'   column `form` records this.
#' @return Tibble with one row per model: `grouping`, `model`, `form`,
#'   `compared_with`, the scaled fit indices, `n`, `converged`, the number of
#'   negative latent and item residual variances, lavaan's warnings, the
#'   changes in fit and `flag` ([invariance_flags()]). The attribute
#'   `collapsed` names the items whose categories were merged.
#' @examples
#' \dontrun{
#' inst <- read_instrument(system.file("instruments/pid5/pid5.yml", package = "normpipe"))
#' invariance_sequence(ref$items, ref$meta$gender, inst, label = "gender")
#' }
#' @export
invariance_sequence <- function(items, group, inst, structure = invariance_structure(inst),
                                levels = c("configural", "scalar", "strict", "second_order"),
                                missing = c("pairwise", "listwise"), min_cat_n = 10,
                                label = "group") {
  if (!requireNamespace("lavaan", quietly = TRUE))
    stop("Package 'lavaan' is needed for the invariance models.")
  missing <- match.arg(missing)
  levels <- match.arg(levels, several.ok = TRUE)
  if (!length(structure$second)) levels <- setdiff(levels, "second_order")
  stopifnot(length(group) == nrow(items), "configural" %in% levels)
  ids <- unlist(structure$first, use.names = FALSE)
  tab <- inst$items[match(ids, inst$items$item), ]
  keep <- !is.na(group)
  g <- droplevels(factor(group[keep]))
  d <- as.data.frame(item_indices(items[keep, , drop = FALSE], inst, ids))
  for (i in seq_along(ids))                      # all items keyed in one direction
    if (tab$reverse[i]) d[[i]] <- (tab$K[i] - 1L) - d[[i]]
  if (missing == "listwise") { ok <- stats::complete.cases(d); d <- d[ok, , drop = FALSE]; g <- droplevels(g[ok]) }
  before <- d
  d[] <- map(d, collapse_sparse, g = g, min_n = min_cat_n)
  collapsed <- names(d)[!map2_lgl(d, before, identical)]
  nthr <- map_int(d, \(x) length(unique(stats::na.omit(x))) - 1L)
  st <- list(first = map(structure$first, \(x) tab$col[match(x, tab$item)]),
             second = structure$second)
  d$.grp <- g; ng <- nlevels(g)
  out <- map_dfr(levels, \(lv)
    invariance_fit(invariance_syntax(st, nthr, ng, lv), d, tab$col, missing) %>%
      mutate(model = lv, form = "higher_order"))
  two <- all(lengths(st$first[unlist(st$second, use.names = FALSE)]) == 2)
  if (is.na(out$cfi.scaled[1]) && length(st$second) && two) {
    alt <- invariance_fit(invariance_syntax(st, nthr, ng, "configural", "domain_residual"),
                          d, tab$col, missing) %>%
      mutate(model = "configural", form = "domain_residual")
    if (!is.na(alt$cfi.scaled)) out[1, ] <- alt
  }
  ref <- c(configural = NA, scalar = "configural", strict = "scalar", second_order = "scalar")
  out <- out %>%
    mutate(grouping = label, compared_with = unname(ref[.data$model])) %>%
    relocate("grouping", "model", "form", "compared_with") %>%
    invariance_flags()
  attr(out, "collapsed") <- collapsed
  out
}

#' Scale-by-scale invariance screen
#'
#' Runs [invariance_sequence()] (configural, scalar, strict) for one
#' unidimensional model per scale, e.g. for every facet of the full
#' instrument in the donor samples.
#'
#' @inheritParams invariance_sequence
#' @param scales Outcome names; default all scales of the level `"facet"` of
#'   the first version.
#' @return Tibble as returned by [invariance_sequence()], with the additional
#'   column `outcome`.
#' @export
invariance_screen <- function(items, group, inst, scales = NULL,
                              missing = c("pairwise", "listwise"), min_cat_n = 10,
                              label = "group") {
  missing <- match.arg(missing)
  if (is.null(scales)) {
    reg <- scale_registry(inst)
    scales <- reg$outcome[reg$version == reg$version[1] & reg$level == "facet"]
  }
  map_dfr(scales, \(o) {
    s <- get_scale(inst, o)
    st <- list(first = set_names(list(as.integer(s$items)), s$id), second = list())
    invariance_sequence(items, group, inst, structure = st,
                        levels = c("configural", "scalar", "strict"),
                        missing = missing, min_cat_n = min_cat_n, label = label) %>%
      mutate(outcome = o, .before = 1)
  })
}

#' Invariance analyses of a norming project
#'
#' The standard set of analyses: the bridge items across the samples
#' (reference sample and donor samples, pairwise and listwise), the bridge
#' items across the reporting cells (gender, age band) within the reference
#' sample, and the scale-by-scale screen across the same cells in the donor
#' samples. With `missing: {method: impute}`, respondents below
#' `min_answered` are left out as in the norming run; skipped items of the
#' others are handled pairwise (not imputed). Writes `invariance_fit.csv`, `invariance_screen.csv` and
#' `invariance_decision.txt` to `out_dir`; none of them contains individual
#' data.
#'
#' @param reference Reference sample: `list(meta, items)` (as for
#'   [run_norming()]).
#' @param donors Donor samples: `list(meta, items)` with `meta$sample`;
#'   `NULL` skips the comparison of the samples and the screen.
#' @param inst Instrument ([read_instrument()]).
#' @param out_dir Output folder; `NULL` writes nothing.
#' @param structure Factor structure of the items answered by the reference
#'   sample; default [invariance_structure()].
#' @param screen_scales Scales of the screen (see [invariance_screen()]);
#'   `FALSE` skips it.
#' @param min_cat_n Passed to [invariance_sequence()] for the reporting cells
#'   (the comparison of the samples merges no categories).
#' @param reference_label Name of the reference sample in the sample grouping.
#' @return List with `fit`, `screen` and `decision` (see
#'   [invariance_decision()]), invisibly.
#' @export
run_invariance <- function(reference, donors = NULL, inst, out_dir = NULL,
                           structure = invariance_structure(inst), screen_scales = NULL,
                           min_cat_n = 10, reference_label = "reference") {
  cells <- cell_vars(inst)
  cell_group <- function(meta, v) {
    meta <- add_age_group(meta, inst)
    g <- meta[[v]]
    if (v == "gender") g[!g %in% inst$cells$gender] <- NA
    factor(g)
  }
  say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " invariance: ", ...)
  # same samples as the norming run: unit non-respondents leave (missing$min_answered)
  mc <- missing_cfg(inst)
  included <- function(x, what) {
    if (is.null(x) || mc$method != "impute") return(x)
    share <- answered_share(x$items, inst)
    out <- is.na(share) | share < mc$min_answered
    if (any(out)) say(sum(out), " ", what, " answered less than ", round(100 * mc$min_answered),
                      " % of the items and are left out")
    list(meta = x$meta[!out, , drop = FALSE], items = x$items[!out, , drop = FALSE])
  }
  reference <- included(reference, "respondent(s) of the reference sample")
  donors <- included(donors, "donor(s)")
  fit <- list()
  if (!is.null(donors)) {
    ids <- unlist(structure$first, use.names = FALSE)
    cols <- inst$items$col[match(ids, inst$items$item)]
    both <- bind_rows(as.data.frame(reference$items)[cols], as.data.frame(donors$items)[cols])
    samp <- factor(c(rep(reference_label, nrow(reference$items)), as.character(donors$meta$sample)),
                   levels = c(reference_label, unique(as.character(donors$meta$sample))))
    for (m in c("pairwise", "listwise")) {
      say("samples, ", m)
      fit[[paste0("samples_", m)]] <- invariance_sequence(
        both, samp, inst, structure, missing = m, min_cat_n = 0,
        label = if (m == "pairwise") "samples" else "samples_listwise")
    }
  }
  for (v in cells) {
    say(v, " (reference sample)")
    fit[[v]] <- invariance_sequence(reference$items, cell_group(reference$meta, v), inst,
                                    structure, min_cat_n = min_cat_n, label = v)
  }
  fit <- bind_rows(fit)
  screen <- NULL
  if (!is.null(donors) && !isFALSE(screen_scales)) {
    screen <- map_dfr(cells, \(v) {
      say("screen by ", v, " (donor samples)")
      invariance_screen(donors$items, cell_group(donors$meta, v), inst, scales = screen_scales,
                        min_cat_n = min_cat_n, label = v)
    })
  }
  decision <- invariance_decision(fit, screen)
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write_csv(fit, file.path(out_dir, "invariance_fit.csv"))
    if (!is.null(screen)) write_csv(screen, file.path(out_dir, "invariance_screen.csv"))
    writeLines(decision$text, file.path(out_dir, "invariance_decision.txt"))
  }
  invisible(list(fit = fit, screen = screen, decision = decision))
}

#' Decision on the invariance analyses
#'
#' @param fit Result of [invariance_sequence()] (one or several groupings
#'   bound together).
#' @param screen Result of [invariance_screen()], or `NULL`.
#' @return List with `ok` (`TRUE` when scalar invariance of the bridge items
#'   is not flagged in any grouping and every configural and scalar model
#'   could be evaluated), `flagged` (the flagged rows of `fit` and `screen`)
#'   and `text` (the summary written to `invariance_decision.txt`).
#' @export
invariance_decision <- function(fit, screen = NULL) {
  main <- fit[fit$model %in% c("scalar", "strict"), ]
  need <- fit[fit$model %in% c("configural", "scalar"), ]
  unavailable <- need[is.na(need$cfi.scaled), ]
  bad_scalar <- main[main$model == "scalar" & !is.na(main$flag) & main$flag, ]
  ok <- nrow(unavailable) == 0 && nrow(bad_scalar) == 0
  fl <- function(x) if (is.null(x)) NULL else x[!is.na(x$flag) & x$flag & x$model != "configural", ]
  num <- function(x) formatC(x, format = "f", digits = 3)
  line <- function(x, lead) sprintf(
    "  %s%s vs. %s: dCFI %s, dRMSEA %s, dSRMR %s%s", lead, x$model, x$compared_with,
    num(x$d_cfi), num(x$d_rmsea), num(x$d_srmr), ifelse(!is.na(x$flag) & x$flag, "  FLAGGED", ""))
  cmp <- fit[!is.na(fit$compared_with), ]
  text <- c(
    "Measurement invariance (WLSMV, theta parameterization, item intercepts 0)",
    "Rule (Chen, 2007): flagged when the CFI falls by >= .010 and the RMSEA rises by >= .015 or the SRMR by >= .010.",
    "", "Items of the factor structure:",
    line(cmp, paste0(format(cmp$grouping), "  ")),
    if (nrow(unavailable)) c("", paste0("  NOT AVAILABLE: ", unavailable$grouping, " / ", unavailable$model)),
    "", paste0("Scalar invariance of the items: ", if (ok) "not rejected in any grouping." else
      "REJECTED or not evaluable in at least one grouping (see above)."))
  if (!is.null(screen)) {
    sc <- screen[!is.na(screen$compared_with), ]
    n_by <- function(m) sprintf("%d of %d", sum(sc$flag[sc$model == m], na.rm = TRUE),
                                sum(sc$model == m))
    fs <- fl(screen)
    text <- c(text, "", "Scale-by-scale screen:",
              paste0("  scalar vs. configural flagged: ", n_by("scalar")),
              paste0("  strict vs. scalar flagged:     ", n_by("strict")),
              if (nrow(fs)) line(fs, paste0(fs$outcome, " by ", fs$grouping, "  ")))
  }
  list(ok = ok, flagged = bind_rows(fl(fit), fl(screen)), text = text)
}
