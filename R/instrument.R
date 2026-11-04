# ==============================================================================
# instrument.R -- the instrument configuration
#
# Everything that was PID-5-specific in the original pipeline (scoring keys,
# versions, response format, reference cells, sample loaders) is read from one
# YAML file per instrument and turned into an object of class
# "normpipe_instrument". Every other part of the package works from that
# object alone. The PID-5 configuration in inst/instruments/pid5 reproduces
# the scale registry (scale_specs) of the original driver.
#
# ITEMS. Every item has an ordered vector of response CATEGORIES (the codes in
# the data file, low to high) and a vector of VALUES, the contribution of each
# category to the raw score. Values default to 0, 1, ..., K - 1, so that
# "categories: [1, 2, 3, 4]" reproduces the 0-3 coding of the PID-5 after
# recoding. Non-uniform values ("values: [2, 1, 0, 1, 2]" for the PDS-ICD-11)
# and item weights (values multiplied by a weight) are allowed. Values are
# stored in integer UNITS (items$unit): with unit 0.5, the values 0, 0.5, 1
# become 0, 1, 2 internally, so that the summed-score recursion of the item
# models runs on integers. Reverse-keyed items have their values reversed.
#
# SCALES. A scale names its items and its metric: "sum" (raw score = sum of
# the item values), "mean" (sum divided by the number of items), or
# "mean_of_means" (mean of the means of its components, the domain rule of
# the PID-5). The attainable raw scores follow from the item values; a
# mean_of_means scale with unequal component lengths has a dense attainable
# set and is tabulated on a reporting grid (report_step) instead. The model
# of a scale ("grm", "bifactor", "composite", "none") determines how the
# individual intervals are computed; see irt.R.
# ==============================================================================

#' Read an instrument configuration
#'
#' @param path Path to the YAML file, or a list with the same structure.
#' @return An object of class normpipe_instrument.
read_instrument <- function(path) {
  cfg <- if (is.character(path)) yaml::yaml.load(paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")) else path
  stopifnot(!is.null(cfg$id), !is.null(cfg$items), !is.null(cfg$versions))
  items <- build_item_table(cfg$items)
  cells <- build_cells(cfg$cells)
  versions <- build_versions(cfg$versions, items, cfg$report_step %||% 0.01)
  structure(list(
    id = cfg$id, name = cfg$name %||% cfg$id,
    items = items, versions = versions, cells = cells,
    reference = cfg$reference %||% list(),
    donors = cfg$donors %||% list(),
    measurement = modifyList(list(source = "fit", engine = "auto",
                                  calibration = "reference"),
                             cfg$measurement %||% list()),
    transfer = modifyList(list(enabled = FALSE), cfg$transfer %||% list()),
    uncertainty = modifyList(
      list(prior_mass = 1, n_boot = 200, n_param_draws = 50, n_draws = 500,
           pctl_clamp = 1e-4, probs = c(.025, .975)),
      cfg$uncertainty %||% list()),
    source = if (is.character(path)) normalizePath(path) else NA_character_
  ), class = "normpipe_instrument")
}

print.normpipe_instrument <- function(x, ...) {
  cat("normpipe instrument '", x$id, "': ", x$name, "\n", sep = "")
  cat("  items: ", nrow(x$items), " (", paste(range(x$items$K), collapse = "-"),
      " categories)\n", sep = "")
  for (v in names(x$versions)) {
    sc <- x$versions[[v]]$scales
    cat("  version '", v, "': ", length(sc), " scales (",
        paste(unique(map_chr(sc, "level")), collapse = ", "), ")\n", sep = "")
  }
  cat("  cells: gender ", if (is.null(x$cells$gender)) "no" else
      paste(x$cells$gender, collapse = "/"), ", age bands ",
      if (is.null(x$cells$age_breaks)) "no" else
        paste(x$cells$age_breaks, collapse = "/"), "\n", sep = "")
  invisible(x)
}

# --- Items ---------------------------------------------------------------------
build_item_table <- function(it) {
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

# --- Reporting cells -----------------------------------------------------------
build_cells <- function(cl) {
  cl <- cl %||% list()
  list(
    gender = if (is.null(cl$gender)) NULL else as.character(cl$gender),
    age_breaks = if (is.null(cl$age_breaks)) NULL else
      as.numeric(unlist(lapply(cl$age_breaks, \(b) if (identical(b, "Inf")) Inf else b))),
    ages = if (is.null(cl$ages)) NULL else as.integer(cl$ages[1]):as.integer(cl$ages[2]),
    overall = cl$overall %||% TRUE,
    kernel = modifyList(list(half_width0 = 4, min_neff = 150, max_half_width = 12),
                        cl$kernel %||% list())
  )
}

#' Cell variables used by the norm tables of an instrument
cell_vars <- function(inst) {
  c(if (!is.null(inst$cells$gender)) "gender",
    if (!is.null(inst$cells$age_breaks)) "age_group")
}

#' Add the age_group column of the reporting cells to a frame
add_age_group <- function(d, inst) {
  if (is.null(inst$cells$age_breaks)) return(d)
  if (!"age_group" %in% names(d) && "age" %in% names(d))
    d$age_group <- cut(d$age, inst$cells$age_breaks, right = FALSE)
  d
}

# --- Versions and scales -------------------------------------------------------
build_versions <- function(vs, items, report_step) {
  out <- list()
  for (v in vs) {
    stopifnot(!is.null(v$id), !is.null(v$scales))
    scales <- list()
    for (s in v$scales) {
      scales[[s$id]] <- build_scale(s, v$id, scales, items, report_step)
    }
    out[[v$id]] <- list(id = v$id, label = v$label %||% v$id, scales = scales)
  }
  out
}

build_scale <- function(s, version, scales_so_far, items, report_step) {
  stopifnot(!is.null(s$id))
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
attainable_sums <- function(vals) {
  reach <- 0L
  for (v in vals) reach <- sort(unique(as.vector(outer(reach, as.integer(v), `+`))))
  reach
}

#' All scales of an instrument as a registry tibble (one row per scale)
scale_registry <- function(inst) {
  map_dfr(inst$versions, \(v) map_dfr(v$scales, \(s)
    tibble(outcome = paste0(v$id, "_", s$id), version = v$id, scale = s$id,
           level = s$level, metric = s$metric, model = s$model,
           n_items = s$n_items, dense = s$dense, n_rows = length(s$grid),
           raw_min = s$range[1], raw_max = s$range[2], step = s$step)))
}

#' Fetch one scale specification by outcome name (version_scale)
get_scale <- function(inst, outcome) {
  for (v in inst$versions) for (s in v$scales)
    if (paste0(v$id, "_", s$id) == outcome) return(s)
  stop("Unknown outcome '", outcome, "'.")
}

#' Which scales can be normed directly, i.e., all their items were answered by
#' the reference sample (observed_items), and which need the transfer module
scale_routes <- function(inst, observed_items) {
  reg <- scale_registry(inst)
  reg$route <- map_chr(reg$outcome, \(o) {
    s <- get_scale(inst, o)
    if (isTRUE(s$observed)) "direct" else if (isFALSE(s$observed)) "transfer"
    else if (all(s$items %in% observed_items)) "direct" else "transfer"
  })
  reg
}
