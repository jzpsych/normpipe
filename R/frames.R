# ==============================================================================
# frames.R -- draw frames: the common input of every norm function
#
# A DRAW FRAME is a long tibble with one row per (respondent or pseudo
# respondent) x draw, columns
#   id, age, gender, weight, score, shape, .draw
# plus age_group when the instrument has age bands. score is on the raw metric
# of the scale. `shape` is the gamma shape of the row in the Bayesian
# bootstrap (1 for an observed respondent, the prior mass of a
# pseudo-observation, the count of a frequency row). Every norm function
# (weighted_norms, weighted_norms_continuous, integrated intervals) consumes
# this format, whichever route produced it:
#
#   direct_frame()      observed scores of the reference sample, augmented by
#                       the Dirichlet prior, Bayesian bootstrap over the
#                       weights (the default route)
#   frequency_frame()   the same from a table of counts per score and cell,
#                       for published or aggregated data (no item data)
#   project_outcome()   posterior predictive draws of a transfer model
#                       (transfer.R, optional)
#
# BAYESIAN BOOTSTRAP WITH A DIRICHLET PRIOR. In every draw
# the weight of a row is replaced by weight / shape * Gamma(shape):
#   - an observed respondent (shape 1) gets weight x Exp(1), the weighted
#     Bayesian bootstrap of Rubin (1981);
#   - a pseudo-observation of the prior (shape = prior mass of the row, weight
#     = shape x the weight of one effective observation of the cell) gets a
#     gamma draw, so that the pseudo part of a draw is a Dirichlet draw with
#     concentration prior_mass / K per attainable score;
#   - a frequency row with count n (shape n, weight = weighted count) gets
#     Gamma(n) scaled, which is the sum of n exponential draws: frequency and
#     individual input give the same posterior.
# Summed over all rows of a cell, the result is a draw from the Dirichlet
# posterior over the attainable scores with parameters (weighted counts +
# prior_mass / K). The percentile rank of a score is computed in every draw
# (norms.R, convention uncertainty$pctl_convention); the median over draws is
# the norm value, the quantiles `probs` of the draws the norm band.
# ==============================================================================

#' Draw frames: observed scores or frequency tables with a Dirichlet prior
#'
#' A draw frame is the common input of all norm functions: a long tibble with
#' one row per respondent (or pseudo-respondent) and draw, with columns `id`,
#' `age`, `gender`, `weight`, `score`, `shape`, `.draw` (plus `age_group` when
#' the instrument has age bands). Every reporting cell receives `prior_mass`
#' pseudo-observations spread uniformly over the attainable scores (a
#' Dirichlet prior with uniform base measure); `n_boot` Bayesian-bootstrap
#' draws then replace each weight by `weight / shape * Gamma(shape)`. Summed
#' over a cell, every draw is a draw from the Dirichlet posterior over the
#' attainable scores, so frequency and individual input give the same
#' posterior.
#'
#' `direct_frame()` takes the scale scores of the reference sample.
#' `frequency_frame()` takes a table of counts per score and cell, so that
#' published or aggregated frequency tables can be normed without raw data.
#'
#' @param ref Tibble with `id`, `weight`, `score` (raw metric; prorated
#'   scores are snapped to the grid) and optionally `age` and `gender`.
#' @param s Scale specification ([get_scale()]).
#' @param inst Instrument ([read_instrument()]).
#' @param n_boot Number of Bayesian-bootstrap draws; 0 gives the point frame
#'   (one draw with the original weights).
#' @param prior_mass Effective number of pseudo-observations per reporting
#'   cell; 0 switches the prior off.
#' @param seed Random seed; derived from the scale name when `NULL`.
#' @param scores_imp Optional matrix of scores from completed data sets
#'   ([impute_items()]), one row per row of `ref` and one column per
#'   imputation. Draw `r` then uses column `((r - 1) mod M) + 1` instead of
#'   `ref$score`, which only decides who is in the sample.
#' @return A draw frame (tibble, see Details).
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' s <- get_scale(inst, "full_b")
#' freq <- data.frame(score = 0:24, count = round(100 * dbinom(0:24, 24, 0.3)))
#' fr <- frequency_frame(freq, s, inst, n_boot = 20)
#' head(weighted_norms(fr, s, inst, cells = character(0)))
#' @export
direct_frame <- function(ref, s, inst, n_boot = 0L,
                         prior_mass = inst$uncertainty$prior_mass,
                         seed = NULL, scores_imp = NULL) {
  stopifnot(all(c("id", "weight", "score") %in% names(ref)))
  keep <- !is.na(ref$score) & !is.na(ref$weight)
  if (!is.null(scores_imp)) {
    scores_imp <- as.matrix(scores_imp)
    stopifnot(nrow(scores_imp) == nrow(ref))
    scores_imp <- scores_imp[keep, , drop = FALSE]
    if (anyNA(scores_imp)) stop("direct_frame(): imputed scores are missing for respondents with a score.")
  }
  obs <- ref %>%
    filter(keep) %>%
    transmute(id = as.character(.data$id),   # prior rows carry character ids
              age = if ("age" %in% names(ref)) .data$age else NA_real_,
              gender = if ("gender" %in% names(ref)) as.character(.data$gender) else NA_character_,
              weight = .data$weight,
              score = snap_to_grid(.data$score, s),
              shape = 1)
  obs <- add_age_group(obs, inst)
  n_person <- nrow(obs)
  if (prior_mass > 0) obs <- bind_rows(obs, prior_rows(obs, s, inst, prior_mass))
  fr <- bootstrap_frame(obs, n_boot, seed %||% sum(utf8ToInt(paste0("direct_", s$version, "_", s$id))))
  if (!is.null(scores_imp)) {
    # draw r takes the scores of completed data set ((r - 1) mod M) + 1
    pos <- (seq_len(nrow(fr)) - 1L) %% nrow(obs) + 1L
    person <- pos <= n_person
    m <- (fr$.draw[person] - 1L) %% ncol(scores_imp) + 1L
    fr$score[person] <- snap_to_grid(scores_imp[cbind(pos[person], m)], s)
  }
  fr
}

#' @rdname direct_frame
#' @param freq Tibble with columns `score` (raw metric), `count` (weighted
#'   counts allowed) and the cell variables of the instrument (`gender`,
#'   `age_group`) or `age` (single years). An optional column `n_obs` holds
#'   the unweighted counts (default `count`); it is the gamma shape of the
#'   row, so that weighted counts keep the right amount of uncertainty. Rows
#'   of cells not covered by the configuration enter the overall table only.
#' @export
frequency_frame <- function(freq, s, inst, n_boot = 0L,
                            prior_mass = inst$uncertainty$prior_mass,
                            seed = NULL) {
  stopifnot(all(c("score", "count") %in% names(freq)))
  obs <- tibble(
    id = paste0("freq_", seq_len(nrow(freq))),
    age = if ("age" %in% names(freq)) as.numeric(freq$age) else NA_real_,
    gender = if ("gender" %in% names(freq)) as.character(freq$gender) else NA_character_,
    weight = as.numeric(freq$count),
    score = snap_to_grid(as.numeric(freq$score), s),
    shape = if ("n_obs" %in% names(freq)) as.numeric(freq$n_obs) else as.numeric(freq$count))
  if ("age_group" %in% names(freq)) {
    obs$age_group <- factor(freq$age_group, levels = age_group_levels(inst))
  } else obs <- add_age_group(obs, inst)
  obs <- obs[obs$weight > 0, ]
  if (prior_mass > 0) obs <- bind_rows(obs, prior_rows(obs, s, inst, prior_mass))
  bootstrap_frame(obs, n_boot, seed %||% sum(utf8ToInt(paste0("freq_", s$version, "_", s$id))))
}

age_group_levels <- function(inst) {
  if (is.null(inst$cells$age_breaks)) return(NULL)
  levels(cut(inst$cells$age_breaks[1], inst$cells$age_breaks, right = FALSE))
}

#' Pseudo-observations of the Dirichlet prior
#'
#' For every reporting cell one row per attainable score (and, when ages are
#' known, per integer age of the band, so that the kernel norms see the same
#' prior). The rows of a cell add up to prior_mass effective observations:
#' shape = prior_mass / K / n_ages, weight = shape x (sum of weights / Kish N)
#' of the cell. Cells are the configured gender x age-band combinations; a
#' frame without gender or age uses the one overall cell.
#' @noRd
prior_rows <- function(obs, s, inst, prior_mass) {
  K <- length(s$grid)
  vars <- intersect(cell_vars(inst), names(obs))
  cells <- obs
  if ("gender" %in% vars) cells <- filter(cells, .data$gender %in% inst$cells$gender)
  if ("age_group" %in% vars) cells <- filter(cells, !is.na(.data$age_group))
  cells <- cells %>%
    group_by(across(all_of(vars))) %>%
    summarise(prior_scale = sum(.data$weight) / (sum(.data$weight)^2 / sum(.data$weight^2)),
              .groups = "drop")
  if (nrow(cells) == 0)       # no cell structure: one overall cell
    cells <- tibble(prior_scale = sum(obs$weight) / (sum(obs$weight)^2 / sum(obs$weight^2)))
  with_ages <- "age_group" %in% vars && !is.null(inst$cells$ages) && any(!is.na(obs$age))
  if (with_ages) {
    age_tbl <- tibble(age = inst$cells$ages) %>%
      mutate(age_group = cut(.data$age, inst$cells$age_breaks, right = FALSE)) %>%
      filter(!is.na(.data$age_group)) %>%
      add_count(.data$age_group, name = "n_ages")
    rows <- inner_join(cells, age_tbl, by = "age_group", relationship = "many-to-many")
  } else {
    rows <- mutate(cells, age = NA_real_, n_ages = 1L)
  }
  rows <- rows %>%
    crossing(score = s$grid) %>%
    mutate(shape = prior_mass / K / .data$n_ages,
           weight = .data$prior_scale * .data$shape,
           id = paste0("prior_", row_number()))
  if (!"gender" %in% names(rows)) rows$gender <- NA_character_
  if (!"age_group" %in% names(rows) && "age_group" %in% names(obs)) rows$age_group <- NA
  rows[, names(obs)]
}

#' Bayesian bootstrap: replicate a point frame n_boot times with gamma weights
#' @noRd
bootstrap_frame <- function(obs, n_boot, seed) {
  if (n_boot < 1L) return(mutate(obs, .draw = 1L))
  set.seed(seed)
  map_dfr(seq_len(n_boot), \(r) {
    g <- rgamma(nrow(obs), shape = obs$shape)
    mutate(obs, weight = .data$weight / .data$shape * g, .draw = r)
  })
}

#' Reduce a draw frame to the point frame (draw 1 with the original weights)
#' @noRd
point_frame <- function(fr) filter(fr, .data$.draw == min(.data$.draw))

#' Compact on-disk form of a frame (person rows once + score/weight matrices)
#' @noRd
frame_compact <- function(fr) {
  fr <- arrange(fr, .data$.draw, .data$id)
  draws <- sort(unique(fr$.draw))
  meta <- fr %>% filter(.data$.draw == draws[1]) %>%
    select(-"score", -"weight", -".draw")
  n <- nrow(meta)
  stopifnot(nrow(fr) == n * length(draws))
  list(meta = meta, draws = draws,
       score = matrix(fr$score, n, length(draws)),
       weight = matrix(fr$weight, n, length(draws)))
}

frame_long <- function(x) {
  n <- nrow(x$meta); D <- length(x$draws)
  bind_cols(x$meta[rep(seq_len(n), times = D), ],
            tibble(.draw = rep(x$draws, each = n),
                   score = as.vector(x$score), weight = as.vector(x$weight)))
}
