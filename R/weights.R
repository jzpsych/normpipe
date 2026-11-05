# ==============================================================================
# weights.R -- calibration of survey weights by raking
#
# The reference sample may carry design weights; they can be raked to
# population margins with the survey package. A margin is a data frame with
# one or more variables and a column Freq (population counts); several
# margins (e.g. age band x gender, education, region) are raked jointly by
# iterative proportional fitting. make_margins() builds the age band x gender
# margin from single-year population counts. Raked weights are trimmed at
# trim_quantile; by default the trimmed weights are raked again (and trimmed
# again) until the margins hold and no weight exceeds the cap, so that
# trimming does not undo the calibration. Weights are finally rescaled to the
# net sample size.
# ==============================================================================

#' Population margins and raking of survey weights
#'
#' `make_margins()` aggregates single-year population counts to the
#' reporting cells of an instrument. `rake_to_margins()` rakes the design
#' weights of the reference sample to one or more population margins with
#' [survey::rake()], trims them at `trim_quantile` and rescales them to the
#' net sample size. [run_norming()] calls it when `margins` is given.
#'
#' Each margin is a data frame whose columns other than `Freq` name the
#' raking variables of that margin (e.g. `age_group` and `gender`, or
#' `education`). `age_group` is derived from `age` and the age bands of the
#' instrument; every other variable must be a column of `meta`. A single data
#' frame is one margin (the joint distribution of its variables), a list of
#' data frames are several margins raked jointly.
#'
#' Trimming: the cap is the `trim_quantile` quantile of the raked weights.
#' With `rerake = TRUE` (default) the capped weights are raked again and
#' capped again, up to `max_rounds` times, until every margin deviates by
#' less than `tol` (relative) and no weight exceeds the cap by more than
#' `tol`. With `rerake = FALSE` the weights are capped once, as in normpipe
#' 0.1.0, and the margins hold only approximately.
#'
#' Respondents with a missing raking variable (or a gender outside the
#' configured cells, when `gender` is a raking variable) cannot be raked.
#' With `incomplete = "drop"` (default) they are dropped; with
#' `incomplete = "keep"` they keep their design weight, rescaled by the mean
#' ratio of raked to design weights, and enter the norms of every table that
#' does not condition on the missing variable (e.g. the overall norms).
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
#'   if none) and any further raking variables.
#' @param margins A margin data frame (e.g. the output of `make_margins()`)
#'   or a list of margin data frames (see Details).
#' @param trim_quantile Raked weights are capped at this quantile; `1`
#'   switches trimming off.
#' @param rerake Rake again after trimming until the margins hold (see
#'   Details).
#' @param max_rounds Maximum number of rake-and-trim rounds.
#' @param tol Relative tolerance for the margins and the cap.
#' @param incomplete `"drop"` or `"keep"` respondents who cannot be raked.
#' @return `rake_to_margins()`: `meta` with column `weight_raked` (rows that
#'   cannot be raked are dropped or kept, see `incomplete`; a logical column
#'   `raked` marks the raked rows) and attribute `"diagnostics"` (Kish
#'   effective N before and after raking, ratio of largest to smallest
#'   weight, number of rounds, share of capped weights, largest relative
#'   deviation from a margin, number of rows not raked).
#' @export
rake_to_margins <- function(meta, margins, inst, trim_quantile = 0.99,
                            rerake = TRUE, max_rounds = 50, tol = 1e-3,
                            incomplete = c("drop", "keep")) {
  incomplete <- match.arg(incomplete)
  if (is.data.frame(margins)) margins <- list(margins)
  margins <- map(margins, as.data.frame)
  vars_by_margin <- map(margins, \(m) setdiff(names(m), "Freq"))
  vars <- unique(unlist(vars_by_margin))
  d <- add_age_group(meta, inst) %>% filter(!is.na(.data$weight))
  miss <- setdiff(vars, names(d))
  if (length(miss)) stop("Raking variable(s) not in meta: ", paste(miss, collapse = ", "))
  ok <- rep(TRUE, nrow(d))
  for (v in vars) ok <- ok & !is.na(d[[v]])
  if ("gender" %in% vars && !is.null(inst$cells$gender))
    ok <- ok & as.character(d$gender) %in% inst$cells$gender
  dr <- d[ok, ]
  # common factor levels of sample and population for every raking variable
  for (v in vars) {
    lv <- unique(unlist(map(margins[map_lgl(vars_by_margin, \(x) v %in% x)], \(m) as.character(m[[v]]))))
    absent <- setdiff(unique(as.character(dr[[v]])), lv)
    if (length(absent)) stop("Raking variable '", v, "': sample categories without a ",
                             "population count: ", paste(absent, collapse = ", "))
    dr[[v]] <- factor(as.character(dr[[v]]), levels = lv)
    for (j in seq_along(margins)) if (v %in% vars_by_margin[[j]])
      margins[[j]][[v]] <- factor(as.character(margins[[j]][[v]]), levels = lv)
  }
  empty <- unlist(map(seq_along(margins), \(j) {
    m <- margins[[j]]; vs <- vars_by_margin[[j]]
    key_m <- do.call(paste, c(m[vs], sep = "|"))
    key_s <- unique(do.call(paste, c(dr[vs], sep = "|")))
    setdiff(key_m[m$Freq > 0], key_s)
  }))
  if (length(empty)) stop("Margin categories without respondents: ", paste(head(empty, 5), collapse = ", "))
  fs <- map(vars_by_margin, \(vs) stats::as.formula(paste("~", paste(vs, collapse = " + "))))
  rake_once <- function(w) {
    dr$.w <- w
    des <- survey::svydesign(ids = ~1, weights = ~.w, data = dr)
    stats::weights(survey::rake(des, sample.margins = fs, population.margins = margins,
                                control = list(maxit = 100, epsilon = 1e-9)))
  }
  w <- rake_once(dr$weight)
  cap <- if (trim_quantile < 1) stats::quantile(w, trim_quantile, names = FALSE) else Inf
  rounds <- 1L
  dev <- margin_deviation(w, dr, margins, vars_by_margin)
  if (is.finite(cap)) {
    repeat {
      w <- pmin(w, cap)
      if (!rerake) break
      dev <- margin_deviation(w, dr, margins, vars_by_margin)
      if (dev < tol || rounds >= max_rounds) break
      w <- rake_once(w)
      rounds <- rounds + 1L
      if (max(w) <= cap * (1 + tol)) break
    }
    if (rerake && rounds >= max_rounds)
      warning("rake_to_margins(): no convergence of rake-and-trim after ", max_rounds,
              " rounds (largest margin deviation ", signif(dev, 3), ").")
  }
  dev <- margin_deviation(w, dr, margins, vars_by_margin)
  n_cap <- sum(w >= cap * (1 - 1e-12))
  ratio <- sum(w) / sum(dr$weight)
  dr$weight_raked <- w
  dr$raked <- TRUE
  out <- dr
  if (incomplete == "keep" && any(!ok)) {
    di <- d[!ok, ]
    di$weight_raked <- di$weight * ratio
    di$raked <- FALSE
    out <- bind_rows(mutate(out, across(all_of(vars), as.character)),
                     mutate(di, across(all_of(vars), as.character)))
  }
  out$weight_raked <- out$weight_raked * nrow(out) / sum(out$weight_raked)
  out$.w <- NULL
  attr(out, "diagnostics") <- list(
    kish_neff_design = kish_neff(out$weight),
    kish_neff_raked = kish_neff(out$weight_raked),
    weight_ratio_max = max(out$weight_raked) / min(out$weight_raked),
    rounds = rounds,
    share_capped = n_cap / nrow(dr),
    max_margin_deviation = dev,
    n_not_raked = sum(!ok))
  out
}

#' Largest relative deviation of the weighted sample from the margins
#' @noRd
margin_deviation <- function(w, d, margins, vars_by_margin) {
  max(map_dbl(seq_along(margins), \(j) {
    m <- margins[[j]]; vs <- vars_by_margin[[j]]
    key_m <- do.call(paste, c(lapply(m[vs], as.character), sep = "|"))
    key_s <- do.call(paste, c(lapply(d[vs], as.character), sep = "|"))
    tot_s <- tapply(w, key_s, sum)[key_m]
    tot_s[is.na(tot_s)] <- 0
    target <- m$Freq / sum(m$Freq) * sum(w)
    keep <- target > 0
    max(abs(tot_s[keep] - target[keep]) / target[keep])
  }))
}
