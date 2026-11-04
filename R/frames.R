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
#                       weights (the default route; Supplement S5.1)
#   frequency_frame()   the same from a table of counts per score and cell,
#                       for published or aggregated data (no item data)
#   project_outcome()   posterior predictive draws of a transfer model
#                       (transfer.R, optional)
#
# BAYESIAN BOOTSTRAP WITH A DIRICHLET PRIOR (Supplement S5.1). In every draw
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
# prior_mass / K). The percentile of a score is the cumulative share at or
# below it; the median over draws is the norm value, the 2.5th and 97.5th
# percentiles the norm band.
# ==============================================================================

#' Observed-score frame of the reference sample, with Dirichlet prior
#'
#' @param ref reference data: tibble with id, age, gender, weight and the
#'   scale score in column `score` (raw metric, prorated allowed)
#' @param s scale specification (get_scale)
#' @param inst instrument
#' @param n_boot number of Bayesian-bootstrap draws (0 = point frame)
#' @param prior_mass effective pseudo-observations per cell (0 = no prior)
#' @param seed seed; derived from the scale name when NULL
direct_frame <- function(ref, s, inst, n_boot = 0L,
                         prior_mass = inst$uncertainty$prior_mass,
                         seed = NULL) {
  stopifnot(all(c("id", "weight", "score") %in% names(ref)))
  obs <- ref %>%
    filter(!is.na(.data$score), !is.na(.data$weight)) %>%
    transmute(id = .data$id,
              age = if ("age" %in% names(ref)) .data$age else NA_real_,
              gender = if ("gender" %in% names(ref)) as.character(.data$gender) else NA_character_,
              weight = .data$weight,
              score = snap_to_grid(.data$score, s),
              shape = 1)
  obs <- add_age_group(obs, inst)
  if (prior_mass > 0) obs <- bind_rows(obs, prior_rows(obs, s, inst, prior_mass))
  bootstrap_frame(obs, n_boot, seed %||% sum(utf8ToInt(paste0("direct_", s$version, "_", s$id))))
}

#' Frame from a frequency table (no item data)
#'
#' @param freq tibble with columns score (raw metric), count (may be weighted
#'   counts), and the cell variables of the instrument (gender, age_group) or
#'   age (single years). Rows of cells not covered by the configuration enter
#'   the overall table only.
#' @param n_obs optional unweighted counts (defaults to count); the gamma
#'   shape of a row, so that weighted counts keep the right amount of
#'   uncertainty
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
    rows <- inner_join(cells, age_tbl, by = "age_group")
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
bootstrap_frame <- function(obs, n_boot, seed) {
  if (n_boot < 1L) return(mutate(obs, .draw = 1L))
  set.seed(seed)
  map_dfr(seq_len(n_boot), \(r) {
    g <- rgamma(nrow(obs), shape = obs$shape)
    mutate(obs, weight = .data$weight / .data$shape * g, .draw = r)
  })
}

#' Reduce a draw frame to the point frame (draw 1 with the original weights)
point_frame <- function(fr) filter(fr, .data$.draw == min(.data$.draw))

#' Compact on-disk form of a frame (person rows once + score/weight matrices)
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
