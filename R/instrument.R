# ==============================================================================
# instrument.R -- the instrument configuration
#
# Everything that is specific to an instrument (items and response format,
# versions and scales, norm groups, norm scale, uncertainty settings,
# measurement and transfer options) is read from one YAML file per instrument and turned into
# an object of class "normpipe_instrument". Every other part of the package
# works from that object alone. The full list of keys is in
# vignette("instrument-config").
#
# ITEMS. Every item has an ordered vector of response CATEGORIES (the codes in
# the data file, low to high) and a vector of VALUES, the contribution of each
# category to the raw score. Values default to 0, 1, ..., K - 1, so that
# "categories: [1, 2, 3, 4]" gives the values 0 to 3. Non-uniform values
# (e.g., "values: [2, 1, 0, 1, 2]" for a bipolar item scored by distance from
# the midpoint) and item weights (values multiplied by a weight) are allowed. Values are
# stored in integer UNITS (items$unit): with unit 0.5, the values 0, 0.5, 1
# become 0, 1, 2 internally, so that the summed-score recursion of the item
# models runs on integers. Reverse-keyed items have their values reversed.
#
# SCALES. A scale names its items and its metric: "sum" (raw score = sum of
# the item values), "mean" (sum divided by the number of items), or
# "mean_of_means" (mean of the means of its components, e.g., a domain
# score as the mean of its facet scores). The attainable raw scores follow from the item values; a
# mean_of_means scale with unequal component lengths has a dense attainable
# set and is tabulated on a reporting grid (report_step) instead. The model
# of a scale ("grm", "bifactor", "composite", "none") determines how the
# individual intervals are computed; see irt.R.
# ==============================================================================

#' Read an instrument configuration
#'
#' Parses an instrument YAML file (or a list with the same structure), builds
#' the item table, the norm groups and every scale of every version
#' (attainable scores, reporting grid, model), and fills in the defaults of
#' the `uncertainty`, `measurement`, `transfer` and `weighting` blocks. The keys are
#' documented in `vignette("instrument-config", package = "normpipe")`.
#'
#' @param path Path to the YAML file, or a list with the same structure.
#' @return An object of class `normpipe_instrument`: a list with `id`, `name`,
#'   `items` (item table: `item`, `col`, `categories`, `values` in integer
#'   units, `reverse`, `K`; attribute `unit`), `versions` (each with its
#'   `scales`), `groups` (the norm groups, see [cell_vars()]), `norm_scale`,
#'   `export`, `measurement`, `transfer`, `weighting`, `uncertainty` and `source` (normalised path of the YAML file).
#' @seealso [scale_registry()], [get_scale()]
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' inst
#' scale_registry(inst)
#' @section Lifecycle:
#' Stable since normpipe 1.0.0: the keys of the configuration file and their
#' defaults are kept, later versions only add keys. The result is a list
#' whose documented elements are kept.
#' @export
read_instrument <- function(path) {
  cfg <- if (is.character(path)) yaml::yaml.load(paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")) else path
  stopifnot(!is.null(cfg$id), !is.null(cfg$items), !is.null(cfg$versions))
  check_requires(cfg$requires)
  check_keys(cfg, c("id", "name", "report_step", "items", "versions", "groups", "cells",
                    "norm_scale", "export", "reference", "donors", "measurement", "transfer",
                    "weighting", "missing", "uncertainty", "notes", "requires"),
             "the instrument")
  check_keys(cfg$export, c("formats", "language"), "`export`")
  check_keys(cfg$missing, c("method", "reference", "n_imputations", "min_answered", "covariates",
                            "neighbours"), "`missing`")
  check_keys(cfg$weighting, c("trim_quantile", "rerake", "incomplete"), "`weighting`")
  check_keys(cfg$uncertainty, c("prior_mass", "n_boot", "n_param_draws", "n_draws", "pctl_clamp",
                                "probs", "pctl_convention", "n_mc", "seed"), "`uncertainty`")
  check_keys(cfg$measurement, c("source", "engine", "calibration"), "`measurement`")
  if (!is.null(cfg$measurement$engine) && !cfg$measurement$engine %in% c("auto", "mirt", "internal"))
    stop("measurement$engine must be auto, mirt or internal.", call. = FALSE)
  if (!is.null(cfg$uncertainty$seed) &&
      !(is.numeric(cfg$uncertainty$seed) && length(cfg$uncertainty$seed) == 1 && !is.na(cfg$uncertainty$seed)))
    stop("uncertainty$seed must be a single number.", call. = FALSE)
  items <- build_item_table(cfg$items)
  groups <- build_groups(cfg$groups, cfg$cells)
  versions <- build_versions(cfg$versions, items, cfg$report_step %||% 0.01)
  inst <- structure(list(
    id = cfg$id, name = cfg$name %||% cfg$id,
    items = items, versions = versions, groups = groups,
    norm_scale = build_norm_scale(cfg$norm_scale),
    export = build_export(cfg$export),
    reference = cfg$reference %||% list(),
    donors = cfg$donors %||% list(),
    measurement = modifyList(list(source = "fit", engine = "auto",
                                  calibration = "reference"),
                             cfg$measurement %||% list()),
    transfer = modifyList(list(enabled = FALSE), cfg$transfer %||% list()),
    weighting = modifyList(list(trim_quantile = 0.99, rerake = TRUE, incomplete = "drop"),
                           cfg$weighting %||% list()),
    missing = cfg$missing %||% list(),
    uncertainty = modifyList(
      list(prior_mass = 1, n_boot = 200, n_param_draws = 50, n_draws = 500,
           pctl_clamp = 1e-4, probs = c(.025, .975),
           pctl_convention = "midpoint", n_mc = 1500L),
      cfg$uncertainty %||% list()),
    notes = cfg$notes,
    source = if (is.character(path)) normalizePath(path) else NA_character_
  ), class = "normpipe_instrument")
  if (!is.null(cfg$transfer)) check_transfer_config(inst)
  inst
}

#' Export settings: formats written by the export stage and the language of
#' the labels in the JSON file
#' @noRd
build_export <- function(x = NULL) {
  x <- x %||% list()
  fm <- as.character(unlist(x$formats %||% c("long", "compact", "json")))
  bad <- setdiff(fm, c("long", "compact", "json"))
  if (length(bad)) stop("export$formats: unknown format(s) ", paste(bad, collapse = ", "),
                        " (known: long, compact, json).")
  lang <- x$language %||% "en"
  if (!lang %in% c("de", "en")) stop("export$language must be \"de\" or \"en\".")
  list(formats = fm, language = lang)
}

#' Report configuration keys the package does not know
#'
#' An unknown key is an error, so that a typo, or a key of a later version
#' of the package, cannot silently fall back to a default.
#' `options(normpipe.unknown_keys = "warn")` turns the error into a warning.
#' @noRd
check_keys <- function(block, known, where) {
  bad <- setdiff(names(block), known)
  if (!length(bad)) return(invisible(NULL))
  msg <- paste0("Unknown key(s) in ", where, ": ", paste(bad, collapse = ", "), ".",
                near_miss(bad, known), " Known keys: ", paste(known, collapse = ", "), ".")
  if (identical(getOption("normpipe.unknown_keys", "error"), "warn")) warning(msg, call. = FALSE) else
    stop(msg, " (A key of a later version of normpipe? options(normpipe.unknown_keys = \"warn\") ",
         "ignores unknown keys with a warning.)", call. = FALSE)
  invisible(NULL)
}

#' " Did you mean ...?" for keys that differ from a known key by a typo
#' @noRd
near_miss <- function(bad, known) {
  if (!length(known)) return("")
  d <- utils::adist(tolower(bad), tolower(known))
  hint <- map_chr(seq_along(bad), \(i) {
    j <- which.min(d[i, ])
    if (d[i, j] <= max(1, min(2, nchar(known[j]) %/% 4))) paste0("`", known[j], "` instead of `", bad[i], "`") else NA_character_
  })
  hint <- hint[!is.na(hint)]
  if (length(hint)) paste0(" Did you mean ", paste(hint, collapse = "; "), "?") else ""
}

#' `requires: "1.1"` (or ">= 1.1"): the version of normpipe an instrument file
#' needs at least
#' @noRd
check_requires <- function(req) {
  if (is.null(req)) return(invisible(NULL))
  v <- trimws(sub("^\\s*>=\\s*", "", as.character(req)))
  need <- tryCatch(package_version(v), error = function(e)
    stop("`requires` must be a version number such as \"1.1\", not '", req, "'.", call. = FALSE))
  have <- utils::packageVersion("normpipe")
  if (have < need)
    stop("This instrument file needs normpipe >= ", v, "; installed is ", have,
         ". Update the package (remotes::install_github(\"jzpsych/normpipe\")).", call. = FALSE)
  invisible(NULL)
}

#' @rdname read_instrument
#' @param x An object of class `normpipe_instrument`.
#' @param ... Ignored.
#' @export
print.normpipe_instrument <- function(x, ...) {
  cat("normpipe instrument '", x$id, "': ", x$name, "\n", sep = "")
  cat("  items: ", nrow(x$items), " (", paste(range(x$items$K), collapse = "-"),
      " categories)\n", sep = "")
  for (v in names(x$versions)) {
    sc <- x$versions[[v]]$scales
    cat("  version '", v, "': ", length(sc), " scales (",
        paste(unique(map_chr(sc, "level")), collapse = ", "), ")\n", sep = "")
  }
  if (!length(x$groups$vars)) cat("  groups: none (overall norms only)\n")
  for (v in x$groups$vars)
    cat("  groups by ", v$name, ": ",
        if (v$type == "categorical") paste(v$levels, collapse = "/") else
          paste0("bands ", paste(v$breaks, collapse = "/")), "\n", sep = "")
  if (!is.null(x$groups$continuous))
    cat("  continuous norms over ", x$groups$continuous$variable, " (",
        paste(range(x$groups$continuous$values), collapse = "-"), ")\n", sep = "")
  cat("  norm scale: ", x$norm_scale$name, " (mean ", x$norm_scale$mean, ", sd ", x$norm_scale$sd, ")\n", sep = "")
  invisible(x)
}

# --- Items ---------------------------------------------------------------------
build_item_table <- function(it) {
  check_keys(it, c("n_items", "ids", "prefix", "unit", "categories", "values", "reverse", "overrides"),
             "`items`")
  for (ov in it$overrides %||% list())
    check_keys(ov, c("items", "categories", "values", "reverse"), "`items$overrides`")
  n <- it$n_items %||% length(it$ids)
  ids <- it$ids %||% seq_len(n)
  prefix <- it$prefix %||% "item_"
  unit <- it$unit %||% 1
  default_cat <- as.numeric(it$categories %||% c(0, 1, 2, 3))
  default_val <- it$values
  reverse <- as.integer(it$reverse %||% integer(0))
  tab <- tibble(item = as.integer(ids), col = paste0(prefix, ids),
                categories = rep(list(default_cat), length(ids)),
                values = rep(list(default_val), length(ids)),
                reverse = ids %in% reverse)
  for (ov in it$overrides %||% list()) {
    rows <- match(as.integer(ov$items), tab$item)
    stopifnot(!anyNA(rows))
    if (!is.null(ov$categories)) tab$categories[rows] <- list(as.numeric(ov$categories))
    if (!is.null(ov$values)) tab$values[rows] <- list(as.numeric(ov$values))
    if (!is.null(ov$reverse)) tab$reverse[rows] <- isTRUE(ov$reverse)
  }
  tab$K <- lengths(tab$categories)
  tab$values <- pmap(list(tab$values, tab$K, tab$reverse), \(v, K, rev) {
    v <- if (is.null(v)) seq_len(K) - 1 else v
    stopifnot(length(v) == K)
    v <- if (rev) rev(v) else v
    u <- v / unit
    if (any(abs(u - round(u)) > 1e-8))
      stop("Item values ", paste(v, collapse = ","), " are not multiples of ",
           "the unit ", unit, ". Set items$unit accordingly.")
    as.integer(round(u))
  })
  if (any(map_lgl(tab$values, \(v) min(v) < 0)))
    stop("Item values must be non-negative (shift them and record the shift ",
         "in the scale's offset if needed).")
  attr(tab, "unit") <- unit
  tab
}

# --- Versions and scales -------------------------------------------------------
build_versions <- function(vs, items, report_step) {
  out <- list()
  for (v in vs) {
    stopifnot(!is.null(v$id), !is.null(v$scales))
    check_keys(v, c("id", "label", "scales"), paste0("version '", v$id, "'"))
    scales <- list()
    for (s in v$scales) {
      scales[[s$id]] <- build_scale(s, v$id, scales, items, report_step)
    }
    out[[v$id]] <- list(id = v$id, label = v$label %||% v$id, scales = scales)
  }
  out
}

# Keys of a scale entry the package uses; further keys are kept as metadata
# (`extra`), with a warning when one looks like a typo of a known key
SCALE_KEYS <- c("id", "label", "level", "metric", "model", "items", "components", "min_answered",
                "observed", "calibrate_from", "equal_block_slopes")

build_scale <- function(s, version, scales_so_far, items, report_step) {
  stopifnot(!is.null(s$id))
  other <- setdiff(names(s), SCALE_KEYS)
  hint <- near_miss(other, SCALE_KEYS)
  if (nzchar(hint))
    warning("Scale '", s$id, "' of version '", version, "': key(s) ", paste(other, collapse = ", "),
            " are not used by normpipe and are kept as metadata.", hint, call. = FALSE)
  for (cp in s$components) if (is.list(cp)) check_keys(cp, c("id", "items"), paste0("a component of scale '", s$id, "'"))
  metric <- s$metric %||% "mean"
  stopifnot(metric %in% c("sum", "mean", "mean_of_means"))
  comps <- NULL
  if (!is.null(s$components)) {
    # components are scales of the same version defined earlier (their id) or
    # inline item lists; stored as named list of item id vectors in order
    comps <- map(s$components, \(cp) {
      if (is.character(cp)) {
        ref <- scales_so_far[[cp]]
        if (is.null(ref)) stop("Scale '", s$id, "': component '", cp,
                               "' is not defined before it.")
        ref$items
      } else as.integer(cp$items)
    })
    names(comps) <- map_chr(s$components, \(cp) if (is.character(cp)) cp else cp$id)
    item_ids <- unlist(comps, use.names = FALSE)
  } else {
    item_ids <- as.integer(s$items)
    if (metric == "mean_of_means")
      stop("Scale '", s$id, "': mean_of_means needs components.")
  }
  rows <- match(item_ids, items$item)
  if (anyNA(rows)) stop("Scale '", s$id, "': unknown items ",
                        paste(item_ids[is.na(rows)], collapse = ","))
  vals <- items$values[rows]
  unit <- attr(items, "unit")
  n <- length(item_ids)
  # attainable summed values (integer units) by convolution of the value sets
  sums <- attainable_sums(vals)
  smax <- sum(map_int(vals, max))
  model <- s$model %||% if (is.null(comps)) "grm" else "bifactor"
  stopifnot(model %in% c("grm", "bifactor", "composite", "none"))
  if (metric == "mean_of_means") {
    lens <- lengths(comps)
    if (length(unique(lens)) == 1 && all(map_lgl(vals, \(v) identical(sort(v), sort(vals[[1]]))))) {
      # equal component lengths and identical value sets (up to reverse
      # keying): the mean of means equals the item mean and the attainable set
      # is finite and short
      grid <- sums * unit / n
      step <- unit / n
      dense <- FALSE
    } else {
      # bounds: mean over components of the component means
      lo <- mean(map_dbl(comps, \(cp) sum(map_dbl(items$values[match(cp, items$item)], min)) * unit / length(cp)))
      hi <- mean(map_dbl(comps, \(cp) sum(map_dbl(items$values[match(cp, items$item)], max)) * unit / length(cp)))
      grid <- round(seq(lo, hi, by = report_step), 10)
      if (abs(tail(grid, 1) - hi) > 1e-9) grid <- c(grid, hi)
      step <- report_step
      dense <- TRUE
    }
  } else {
    grid <- if (metric == "sum") sums * unit else sums * unit / n
    step <- if (metric == "sum") unit else unit / n
    dense <- FALSE
    # with non-uniform value sets the attainable set may have gaps; the step
    # is the smallest spacing (used for censoring and prior intervals)
    if (length(grid) > 1) step <- min(diff(grid))
  }
  list(
    id = s$id, version = version, label = s$label %||% s$id,
    level = s$level %||% "scale", metric = metric, model = model,
    items = item_ids, n_items = n, components = comps,
    values = vals, unit = unit, sum_max = smax, sums = sums,
    grid = grid, step = step, dense = dense,
    range = c(min(grid), max(grid)),
    min_answered = s$min_answered %||% 1,
    observed = s$observed,       # NULL = decided from the reference data
    extra = s[setdiff(names(s), c("id", "label", "level", "metric", "model",
                                  "items", "components", "min_answered",
                                  "observed"))]
  )
}

#' Attainable summed values (integer units) of a set of items
#' @param vals list of integer value vectors, one per item
#' @noRd
attainable_sums <- function(vals) {
  reach <- 0L
  for (v in vals) reach <- sort(unique(as.vector(outer(reach, as.integer(v), `+`))))
  reach
}

#' Registry of all scales of an instrument
#'
#' @param inst Instrument ([read_instrument()]).
#' @return Tibble with one row per scale: `outcome` (`<version>_<scale>`, the
#'   name used throughout the package), `version`, `scale`, `level`, `metric`,
#'   `model`, `n_items`, `dense` (tabulated on the reporting grid), `n_rows`
#'   (rows of the norm table), `raw_min`, `raw_max`, `step`.
#' @section Lifecycle:
#' Stable since normpipe 1.0.0: the name, the arguments and the documented
#' parts of the result are kept; later versions only add to them.
#' @export
scale_registry <- function(inst) {
  map_dfr(inst$versions, \(v) map_dfr(v$scales, \(s)
    tibble(outcome = paste0(v$id, "_", s$id), version = v$id, scale = s$id,
           level = s$level, metric = s$metric, model = s$model,
           n_items = s$n_items, dense = s$dense, n_rows = length(s$grid),
           raw_min = s$range[1], raw_max = s$range[2], step = s$step)))
}

#' Specification of one scale
#'
#' @param inst Instrument ([read_instrument()]).
#' @param outcome Outcome name `<version>_<scale>` (see [scale_registry()]).
#' @return List with `id`, `version`, `label`, `level`, `metric`, `model`,
#'   `items`, `n_items`, `components`, `values`, `unit`, `sum_max`, `sums`
#'   (attainable summed values in integer units), `grid` (raw scores of the
#'   norm table), `step`, `dense`, `range`, `min_answered`, `observed` and
#'   `extra` (further keys of the scale, e.g. `calibrate_from`).
#' @section Lifecycle:
#' Stable since normpipe 1.0.0: the name, the arguments and the documented
#' parts of the result are kept; later versions only add to them.
#' @export
get_scale <- function(inst, outcome) {
  for (v in inst$versions) for (s in v$scales)
    if (paste0(v$id, "_", s$id) == outcome) return(s)
  stop("Unknown outcome '", outcome, "'.")
}

#' Route of every scale: direct norming or transfer
#'
#' A scale is normed directly when the reference sample answered all of its
#' items, otherwise it needs the transfer module. The scale key `observed`
#' (`true`/`false`) overrides this decision.
#'
#' @param inst Instrument ([read_instrument()]).
#' @param observed_items Item ids answered by the reference sample.
#' @return [scale_registry()] with an additional column `route`
#'   (`"direct"` or `"transfer"`).
#' @section Lifecycle:
#' Stable since normpipe 1.0.0: the name, the arguments and the documented
#' parts of the result are kept; later versions only add to them.
#' @export
scale_routes <- function(inst, observed_items) {
  reg <- scale_registry(inst)
  reg$route <- map_chr(reg$outcome, \(o) {
    s <- get_scale(inst, o)
    if (isTRUE(s$observed)) "direct" else if (isFALSE(s$observed)) "transfer"
    else if (all(s$items %in% observed_items)) "direct" else "transfer"
  })
  reg
}
