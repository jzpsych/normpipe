#' normpipe: population-referenced norms with two layers of uncertainty
#'
#' normpipe derives norm tables for self-report questionnaires from a weighted
#' reference sample (or from published frequency tables) and reports two kinds
#' of uncertainty: the uncertainty of the norm values themselves (norm bands
#' from a Bayesian bootstrap under a Dirichlet prior over the attainable
#' scores) and the measurement error of the individual score (score-conditional
#' and integrated intervals from calibrated item response models).
#'
#' An instrument is described by a YAML configuration (see
#' `vignette("instrument-config", package = "normpipe")`); [run_norming()]
#' then runs the pipeline stages and writes all tables to a folder. The
#' getting-started guide is `vignette("normpipe", package = "normpipe")`.
#'
#' @section Main functions:
#' * Configuration: [read_instrument()], [scale_registry()], [get_scale()],
#'   [scale_routes()]
#' * Scoring: [score_scales()], [item_indices()], [item_values()]
#' * Weights: [make_margins()], [rake_to_margins()]
#' * Norms: [direct_frame()], [frequency_frame()], [weighted_norms()],
#'   [norms_all_cells()], [weighted_norms_continuous()]
#' * Item models and intervals: [fit_scale_model()], [draw_model_params()],
#'   [make_blocks()], [conditional_tables()], [integrated_intervals()],
#'   [marginal_reliability()]
#' * Transfer module: [bridge_scores()], [fit_transfer()], [project_outcome()]
#' * Export: [long_table_scale()], [export_platform_tables()],
#'   [export_portal_json()]
#' * Driver and test data: [run_norming()], [simulate_reference()]
#'
#' @keywords internal
#' @import dplyr
#' @import purrr
#' @import tibble
#' @importFrom tidyr pivot_longer pivot_wider crossing
#' @importFrom readr read_csv write_csv
#' @importFrom rlang .data %||%
#' @importFrom stats approx qnorm pnorm plogis quantile median rgamma rnorm runif sd setNames optim cor cov2cor rexp weighted.mean var
#' @importFrom utils head tail modifyList
"_PACKAGE"

# symbols used inside brms::prior() expressions
utils::globalVariables(c("normal", "b"))
