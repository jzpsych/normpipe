# ==============================================================================
# groups.R -- norm groups (reporting cells) and the norm scale
#
# NORM GROUPS. The `groups` block of the instrument names the variables of
# the reference sample that define the groups for which norms are tabulated.
# A variable is either CATEGORICAL (its listed codes form the groups; column
# of the same name) or BANDED (a numeric variable cut into bands, closed on
# the left; column `<variable>_group`). A reporting cell is one combination
# of the groups of all variables; the overall table pools everybody.
# Optionally, norms are also computed as a smooth function of one CONTINUOUS
# variable (kernel norms, see norms.R).
#
#   groups:
#     by:
#       gender: [m, f]                        # categorical
#       age: {breaks: [18, 40, 65, Inf]}      # banded -> column age_group
#     overall: true
#     continuous: {variable: age, range: [18, 80],
#                  kernel: {half_width0: 4, min_neff: 150, max_half_width: 12}}
#
# Internally: inst$groups = list(vars = <named list of variable specs>,
# overall, continuous). A variable spec is list(name, type, levels | breaks,
# col).
#
# NORM SCALE. Percentile ranks are turned into standard scores
# mean + sd * qnorm(PR). The scale is T (50, 10) unless `norm_scale` says
# otherwise; its name is the name of the score columns in every table.
# ==============================================================================

# --- Configuration ------------------------------------------------------------

#' Build the groups of an instrument from the `groups` block (or from the
#' `cells` block of normpipe < 0.4)
#' @noRd
build_groups <- function(groups = NULL, cells = NULL) {
  if (!is.null(groups) && !is.null(cells))
    stop("The instrument has both `groups` and `cells`; keep `groups` only.")
  num <- function(x) as.numeric(unlist(lapply(x, \(b) if (identical(b, "Inf")) Inf else if (identical(b, "-Inf")) -Inf else b)))
  if (is.null(groups) && !is.null(cells)) {
    # cells: {gender, age_breaks, ages, overall, kernel}
    check_keys(cells, c("gender", "age_breaks", "ages", "overall", "kernel"), "`cells`")
    by <- list()
    if (!is.null(cells$gender)) by$gender <- cells$gender
    if (!is.null(cells$age_breaks)) by$age <- list(breaks = cells$age_breaks)
    groups <- list(by = by, overall = cells$overall,
                   continuous = if (!is.null(cells$ages)) list(variable = "age", range = cells$ages,
                                                               kernel = cells$kernel))
  }
  groups <- groups %||% list()
  check_keys(groups, c("by", "overall", "continuous"), "`groups`")
  by <- groups$by %||% list()
  if (length(by) && (is.null(names(by)) || any(!nzchar(names(by)))))
    stop("`groups$by` must be a mapping of variable names to their groups.")
  vars <- imap(by, \(spec, name) {
    if (is.list(spec) && !is.null(spec$breaks)) {
      check_keys(spec, "breaks", paste0("`groups$by$", name, "`"))
      br <- num(spec$breaks)
      if (length(br) < 2 || is.unsorted(br, strictly = TRUE))
        stop("groups$by$", name, ": `breaks` must be at least two increasing numbers.")
      list(name = name, type = "banded", breaks = br, col = paste0(name, "_group"))
    } else {
      lv <- as.character(unlist(spec))
      if (!length(lv)) stop("groups$by$", name, ": list the codes that form the groups, or give `breaks`.")
      list(name = name, type = "categorical", levels = lv, col = name)
    }
  })
  cont <- NULL
  if (!is.null(groups$continuous)) {
    cc <- groups$continuous
    check_keys(cc, c("variable", "range", "by", "kernel"), "`groups$continuous`")
    if (is.null(cc$variable) || length(cc$range) != 2)
      stop("`groups$continuous` needs `variable` and `range: [from, to]`.")
    cat_vars <- names(vars)[map_lgl(vars, \(v) v$type == "categorical")]
    cby <- if (is.null(cc$by)) cat_vars else as.character(unlist(cc$by))
    if (length(setdiff(cby, cat_vars)))
      stop("groups$continuous$by: ", paste(setdiff(cby, cat_vars), collapse = ", "),
           " is not a categorical variable of `groups$by`.")
    cont <- list(variable = as.character(cc$variable),
                 values = as.integer(cc$range[[1]]):as.integer(cc$range[[2]]),
                 by = cby,
                 kernel = modifyList(list(half_width0 = 4, min_neff = 150, max_half_width = 12),
                                     cc$kernel %||% list()))
  }
  list(vars = vars, overall = groups$overall %||% TRUE, continuous = cont)
}

# --- Accessors ------------------------------------------------------------------

#' Variables of the norm groups
#'
#' `cell_vars()` returns the columns that identify a reporting cell: the
#' name of every categorical group variable and `<variable>_group` for every
#' banded one. `group_sources()` returns the columns of the reference data
#' that the groups (and the continuous norms) are built from.
#'
#' @param inst Instrument ([read_instrument()]).
#' @return Character vector (possibly empty).
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' cell_vars(inst)
#' group_sources(inst)
#' @section Lifecycle:
#' Stable since normpipe 1.0.0: the name, the arguments and the documented
#' parts of the result are kept; later versions only add to them.
#' @export
cell_vars <- function(inst) unname(map_chr(inst$groups$vars, "col"))

#' @rdname cell_vars
#' @export
group_sources <- function(inst) {
  g <- inst$groups
  num <- c(g$continuous$variable, names(g$vars)[map_lgl(g$vars, \(v) v$type == "banded")])
  unique(c(num, names(g$vars)[map_lgl(g$vars, \(v) v$type == "categorical")]))
}

#' Group variable spec by its cell column
#' @noRd
group_of_col <- function(inst, col) {
  hit <- Filter(\(v) v$col == col, inst$groups$vars)
  if (length(hit)) hit[[1]] else NULL
}

#' Numeric source variables (banded or continuous) of an instrument
#' @noRd
numeric_sources <- function(inst) {
  g <- inst$groups
  unique(c(g$continuous$variable, names(g$vars)[map_lgl(g$vars, \(v) v$type == "banded")]))
}

#' Categorical source variables of an instrument
#' @noRd
categorical_sources <- function(inst)
  names(inst$groups$vars)[map_lgl(inst$groups$vars, \(v) v$type == "categorical")]

#' Levels of a band column
#' @noRd
band_levels <- function(v) levels(cut(v$breaks[1], v$breaks, right = FALSE))

#' Add the band columns of the norm groups to a data frame
#'
#' @param d Data frame with the numeric variables that are cut into bands.
#' @param inst Instrument ([read_instrument()]).
#' @return `d` with a factor column `<variable>_group` (bands closed on the
#'   left, as configured in `groups$by`) for every banded variable that is
#'   in `d`; existing band columns are kept.
#' @section Lifecycle:
#' `add_groups()` is stable since normpipe 1.0.0: the name, the arguments
#' and the documented parts of the result are kept; later versions only add
#' to them. `add_age_group()` is superseded: it is the name of this function
#' in normpipe < 0.4 and is kept as an alias.
#' @export
add_groups <- function(d, inst) {
  for (v in inst$groups$vars) {
    if (v$type != "banded") next
    if (!v$col %in% names(d) && v$name %in% names(d))
      d[[v$col]] <- cut(d[[v$name]], v$breaks, right = FALSE)
  }
  d
}

#' Rows of a data frame that belong to a configured group of every variable
#' in `cols` (cell columns)
#' @noRd
in_groups <- function(d, inst, cols = cell_vars(inst)) {
  ok <- rep(TRUE, nrow(d))
  for (cn in cols) {
    v <- group_of_col(inst, cn)
    if (is.null(v)) next
    ok <- ok & if (v$type == "categorical") as.character(d[[cn]]) %in% v$levels else !is.na(d[[cn]])
  }
  ok
}

#' Columns that every draw frame carries besides id, weight, score, shape:
#' the sources of the groups and of the continuous norms
#' @noRd
frame_vars <- function(inst) group_sources(inst)

#' The frame variables of a data frame as columns of the right type (NA when
#' the variable is absent)
#' @noRd
frame_columns <- function(d, inst) {
  cat_v <- categorical_sources(inst)
  out <- lapply(frame_vars(inst), \(v) {
    if (!v %in% names(d)) return(if (v %in% cat_v) rep(NA_character_, nrow(d)) else rep(NA_real_, nrow(d)))
    if (v %in% cat_v) as.character(d[[v]]) else as.numeric(d[[v]])
  })
  if (!length(out)) return(tibble(.rows = nrow(d)))
  as_tibble(stats::setNames(out, frame_vars(inst)))
}

# --- Norm scale ----------------------------------------------------------------

NORM_SCALES <- list(
  T = list(name = "T", mean = 50, sd = 10, digits = 1),
  z = list(name = "z", mean = 0, sd = 1, digits = 2),
  IQ = list(name = "IQ", mean = 100, sd = 15, digits = 0),
  stanine = list(name = "stanine", mean = 5, sd = 2, limits = c(1, 9), digits = 0),
  sten = list(name = "sten", mean = 5.5, sd = 2, limits = c(1, 10), digits = 0)
)

#' Build the norm scale from the `norm_scale` key: a preset name or a list
#' @noRd
build_norm_scale <- function(x = NULL) {
  if (is.null(x)) x <- "T"
  if (is.character(x)) {
    if (!x %in% names(NORM_SCALES))
      stop("norm_scale '", x, "' is not a preset (", paste(names(NORM_SCALES), collapse = ", "),
           "); give `{name, mean, sd}` instead.")
    return(NORM_SCALES[[x]])
  }
  check_keys(x, c("name", "mean", "sd", "limits", "digits"), "`norm_scale`")
  base <- if (!is.null(x$name) && x$name %in% names(NORM_SCALES)) NORM_SCALES[[x$name]] else list(digits = 1)
  sc <- modifyList(base, x)
  if (is.null(sc$name) || is.null(sc$mean) || is.null(sc$sd))
    stop("`norm_scale` needs `name`, `mean` and `sd`.")
  if (!grepl("^[A-Za-z][A-Za-z0-9]*$", sc$name) ||
      sc$name %in% c("raw", "n", "cell", "scale", "version", "level", "route", "pctl", "score", "weight"))
    stop("norm_scale$name '", sc$name, "' cannot be used as a column name.")
  if (!is.numeric(sc$sd) || sc$sd <= 0) stop("norm_scale$sd must be positive.")
  if (!is.null(sc$limits)) sc$limits <- as.numeric(unlist(sc$limits))
  sc
}

#' Percentile rank -> standard score of the norm scale
#' @noRd
pctl_to_norm <- function(p, clamp = 1e-4, scale = NORM_SCALES$T) {
  x <- scale$mean + scale$sd * qnorm(pmin(pmax(p, clamp), 1 - clamp))
  if (!is.null(scale$limits)) x <- pmin(pmax(x, scale$limits[1]), scale$limits[2])
  x
}

#' Grid of standard scores on which posteriors are mixed (mean +/- 4 sd)
#' @noRd
norm_grid <- function(scale) seq(scale$mean - 4 * scale$sd, scale$mean + 4 * scale$sd, by = scale$sd / 100)

#' Name of the standard-score columns
#' @noRd
norm_name <- function(inst) inst$norm_scale$name %||% "T"

# The package computes with the column prefix "T"; tables that leave a
# function carry the configured name, tables that enter one are turned back.
#' @noRd
norm_public <- function(x, inst) rename_norm(x, "T", norm_name(inst))
#' @noRd
norm_internal <- function(x, inst) rename_norm(x, norm_name(inst), "T")

rename_norm <- function(x, from, to) {
  if (identical(from, to) || is.null(x)) return(x)
  if (inherits(x, "data.frame")) {
    hit <- grepl(paste0("^", from, "($|_)"), names(x))
    names(x)[hit] <- paste0(to, substring(names(x)[hit], nchar(from) + 1L))
    return(x)
  }
  if (is.list(x)) return(lapply(x, rename_norm, from = from, to = to))
  x
}

#' @rdname add_groups
#' @details `add_age_group()` is the name of this function in normpipe < 0.4
#'   and is kept as an alias.
#' @export
add_age_group <- function(d, inst) add_groups(d, inst)
