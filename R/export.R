# ==============================================================================
# export.R -- the published tables
#
# norm_tables.csv              one row per version, scale, reporting cell and
#                              raw score: percentile rank, standard score with
#                              norm band, score-conditional and integrated
#                              interval, extrapolation flag
# norm_tables_continuous.csv   kernel norms over the continuous variable
# compact/norms_<version>.csv  compact per-version tables for scoring software
# json/<instrument>.json       norm file in a JSON format for web applications
#
# The standard-score columns carry the name of the norm scale of the
# instrument; they are written here for the default, T.
# ==============================================================================

#' Export formats
#'
#' `long_table_scale()` assembles the long norm table of one scale: one row
#' per reporting cell and raw score with percentile rank, standard score,
#' norm band, score-conditional and integrated interval (interpolated to
#' every raw score) and an extrapolation flag (median percentile rank below
#' 0.5 % or above 99.5 %). [run_norming()] binds these tables to
#' `norm_tables.csv`.
#'
#' `export_compact()` writes compact per-version tables for scoring
#' software: `norms_<version>.csv` (raw score, standard score, integrated
#' interval, percentile rank, extrapolation flag; the cell variables are
#' `"all"` in the overall table) and, with continuous norms,
#' `norms_<version>_continuous.rds`. Scores are rounded to the `digits` of
#' the norm scale.
#'
#' `export_json()` writes one JSON norm file per instrument for web
#' applications: one source per version, the reporting cells as groups, and
#' per row raw score, cell, standard score and interval limits. Keys and
#' labels are English or German (`export$language` of the instrument).
#'
#' `export_platform_tables()` and `export_portal_json()` are the names of
#' these two functions in normpipe < 0.4 and are kept as aliases.
#'
#' The standard-score columns carry the name of the norm scale of the
#' instrument (`norm_scale`); they are documented for the default, `T`.
#'
#' @param nm Band norms of the scale ([norms_all_cells()]).
#' @param cond Score-conditional tables ([conditional_tables()]) or `NULL`.
#' @param intg Integrated intervals ([integrated_intervals()]) or `NULL`.
#' @param s Scale specification ([get_scale()]).
#' @param inst Instrument ([read_instrument()]).
#' @param route Route of the scale, written into the table.
#' @return `long_table_scale()`: tibble with `version`, `scale`, `level`,
#'   `route`, `cell`, the cell variables, `raw`, `pctl_med`, `pctl_lo`,
#'   `pctl_hi`, `T`, `T_lo`, `T_hi`, `T_cond_lo`, `T_cond_hi`, `T_int_lo`,
#'   `T_int_hi`, `extrapolation`, `n_eff`.
#' @export
long_table_scale <- function(nm, cond, intg, s, inst, route = "direct") {
  nm <- norm_internal(nm, inst); cond <- norm_internal(cond, inst); intg <- norm_internal(intg, inst)
  cv <- cell_vars(inst)
  rows <- list()
  if (!is.null(nm$cells))
    rows$cells <- nm$cells %>% mutate(cell = apply(as.data.frame(nm$cells[cv]), 1, paste, collapse = " "))
  if (!is.null(nm$overall)) rows$overall <- nm$overall %>% mutate(cell = "overall")
  out <- bind_rows(rows) %>%
    mutate(version = s$version, scale = s$id, level = s$level, route = route,
           extrapolation = .data$pctl_med < .005 | .data$pctl_med > .995,
           T_cond_lo = NA_real_, T_cond_hi = NA_real_,
           T_int_lo = NA_real_, T_int_hi = NA_real_)
  attach <- function(out, tabs, lo, hi) {
    for (key in names(tabs)) {
      tb <- tabs[[key]]; sel <- out$cell == key
      if (!any(sel) || is.null(tb) || nrow(tb) < 1) next
      ok <- !is.na(tb$T_lo)
      if (sum(ok) >= 2) {
        out[[lo]][sel] <- approx(tb$raw[ok], tb$T_lo[ok], xout = out$raw[sel], rule = 2)$y
        out[[hi]][sel] <- approx(tb$raw[ok], tb$T_hi[ok], xout = out$raw[sel], rule = 2)$y
      } else if (sum(ok) == 1) { out[[lo]][sel] <- tb$T_lo[ok]; out[[hi]][sel] <- tb$T_hi[ok] }
    }
    out
  }
  if (!is.null(cond)) out <- attach(out, cond, "T_cond_lo", "T_cond_hi")
  if (!is.null(intg)) out <- attach(out, intg, "T_int_lo", "T_int_hi")
  for (v in cv) if (!v %in% names(out)) out[[v]] <- NA
  out %>% select(all_of(c("version", "scale", "level", "route", "cell", cv, "raw",
                          "pctl_med", "pctl_lo", "pctl_hi", "T", "T_lo", "T_hi",
                          "T_cond_lo", "T_cond_hi", "T_int_lo", "T_int_hi",
                          "extrapolation", "n_eff"))) %>%
    norm_public(inst)
}

#' @rdname long_table_scale
#' @param long Long norm table (rows of `long_table_scale()`).
#' @param cont Continuous norms ([weighted_norms_continuous()] with
#'   `version`, `scale`, `level`) or `NULL`.
#' @param dir Output folder.
#' @return `export_compact()`: paths of the written files (invisibly).
#' @export
export_compact <- function(long, cont, dir, inst = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  nn <- if (is.null(inst)) "T" else norm_name(inst)
  digits <- if (is.null(inst)) 1 else inst$norm_scale$digits %||% 1
  long <- rename_norm(long, nn, "T")
  # cell variables: as configured, or the columns between `cell` and `raw`
  cv <- if (!is.null(inst)) cell_vars(inst) else {
    pos <- match(c("cell", "raw"), names(long))
    if (anyNA(pos) || pos[2] - pos[1] < 2) character(0) else names(long)[(pos[1] + 1):(pos[2] - 1)]
  }
  band <- long
  for (v in cv) band[[v]] <- if (v %in% names(long)) coalesce(as.character(long[[v]]), "all") else "all"
  band <- band %>%
    transmute(.data$version, .data$level, .data$scale, across(all_of(cv)),
              raw = round(.data$raw, 4), T = round(.data$T, digits),
              T_int_lo = round(.data$T_int_lo, digits), T_int_hi = round(.data$T_int_hi, digits),
              pctl_med = round(.data$pctl_med, 4),
              extrapolation = as.integer(.data$extrapolation)) %>%
    arrange(.data$version, .data$level, .data$scale, across(all_of(cv)), .data$raw)
  for (v in unique(band$version))
    write_csv(rename_norm(filter(band, .data$version == v), "T", nn),
              file.path(dir, paste0("norms_", v, ".csv")), na = "")
  if (!is.null(cont) && nrow(cont)) {
    cont <- rename_norm(cont, nn, "T")
    # the continuous variable and its `by` variables: the columns before
    # `half_width` other than the scale identifiers
    pos <- match("half_width", names(cont))
    lead <- setdiff(names(cont)[seq_len(pos - 1)], c("version", "level", "scale"))
    xv <- lead[length(lead)]; by <- setdiff(lead, xv)
    cc <- cont
    for (b in by) cc[[b]] <- coalesce(as.character(cc[[b]]), "all")
    cc[[xv]] <- as.integer(cc[[xv]])
    cc <- cc %>%
      transmute(.data$version, .data$level, .data$scale, across(all_of(c(by, xv))),
                raw = round(.data$raw, 4), T = round(.data$T, digits),
                pctl_med = round(.data$pctl_med, 4),
                extrapolation = as.integer(.data$pctl_med < .005 | .data$pctl_med > .995),
                .data$half_width) %>%
      arrange(.data$version, .data$level, .data$scale, across(all_of(c(by, xv))), .data$raw)
    for (v in unique(cc$version))
      saveRDS(as.data.frame(rename_norm(filter(cc, .data$version == v), "T", nn)),
              file.path(dir, paste0("norms_", v, "_continuous.rds")), compress = "xz")
  }
  invisible(list.files(dir, full.names = TRUE))
}

#' @rdname long_table_scale
#' @export
export_platform_tables <- function(long, cont, dir, inst = NULL) export_compact(long, cont, dir, inst)

# keys and labels of the JSON norm file
JSON_WORDS <- list(
  de = c(created = "erzeugt", grouping = "bezugsgruppen", sources = "quellen", type = "art",
         bayesian = "bayesianisch", title = "titel", citation = "zitation", score = "kennwert",
         score_digits = "dezimalen_kennwert", interval = "intervall", limits = "grenzen",
         level = "niveau", description = "erklaerung", scales = "skalen", groups = "gruppen",
         columns = "spalten", rows = "zeilen", raw = "roh", group = "gruppe", lo = "unten",
         hi = "oben", overall = "gesamt", Overall = "Gesamt", gender = "geschlecht", age = "alter",
         default_title = "Bayesianische Normen"),
  en = c(created = "created", grouping = "grouping", sources = "sources", type = "type",
         bayesian = "bayesian", title = "title", citation = "citation", score = "score",
         score_digits = "score_digits", interval = "interval", limits = "limits",
         level = "level", description = "description", scales = "scales", groups = "groups",
         columns = "columns", rows = "rows", raw = "raw", group = "group", lo = "lo",
         hi = "hi", overall = "overall", Overall = "Overall", gender = "gender", age = "age",
         default_title = "Bayesian norms"))

json_interval_text <- function(language, interval, level, nn) {
  if (language == "de") switch(interval,
    integrated = sprintf("%d-%%-Intervall f\u00fcr den wahren %s-Wert der Person, einschlie\u00dflich der Unsicherheit der Normwerte", level, nn),
    conditional = sprintf("%d-%%-Intervall f\u00fcr den wahren %s-Wert der Person bei festen Normwerten", level, nn),
    band = sprintf("%d-%%-Glaubw\u00fcrdigkeitsband des Normwerts", level))
  else switch(interval,
    integrated = sprintf("%d %% interval for the person's true %s score, including the uncertainty of the norm values", level, nn),
    conditional = sprintf("%d %% interval for the person's true %s score, given the norm values", level, nn),
    band = sprintf("%d %% credible band of the norm value", level))
}

#' @rdname long_table_scale
#' @param file Output file.
#' @param citation Citation string of the norms.
#' @param title Title prefix of every source (followed by the version
#'   label); default "Bayesian norms" in the language of the file.
#' @param interval Interval written per row: `"integrated"`, `"conditional"`
#'   or `"band"`.
#' @param language `"en"` or `"de"`: language of the keys and labels;
#'   default `export$language` of the instrument.
#' @return `export_json()`: `file` (invisibly).
#' @export
export_json <- function(long, inst, file, citation = "", title = NULL,
                        interval = c("integrated", "conditional", "band"),
                        language = inst$export$language %||% "en") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")
  interval <- match.arg(interval)
  w <- function(k) unname(JSON_WORDS[[language]][k])
  wv <- function(v) if (v %in% names(JSON_WORDS[[language]])) w(v) else v   # variable names
  nn <- norm_name(inst); digits <- inst$norm_scale$digits %||% 1
  long <- norm_internal(long, inst)
  title <- title %||% w("default_title")
  level <- as.integer(round(100 * diff(inst$uncertainty$probs)))
  lo <- switch(interval, integrated = "T_int_lo", conditional = "T_cond_lo", band = "T_lo")
  hi <- switch(interval, integrated = "T_int_hi", conditional = "T_cond_hi", band = "T_hi")
  cv <- cell_vars(inst)
  src_name <- function(cn) group_of_col(inst, cn)$name %||% cn
  sources <- map(split(long, long$version), \(lv) {
    v <- lv$version[1]
    scales <- map(split(lv, lv$scale), \(ls) {
      s <- get_scale(inst, paste0(v, "_", ls$scale[1]))
      cells <- split(ls, ls$cell)
      groups <- map(cells, \(cd) {
        cell <- cd$cell[1]
        g <- list(id = gsub("[^A-Za-z0-9]+", "_", cell),
                  label = if (cell == "overall") w("Overall") else cell,
                  n_eff = round(cd$n_eff[1]))
        if (cell != "overall") for (cn in cv) g[[wv(src_name(cn))]] <- as.character(cd[[cn]][1])
        g
      })
      rows <- map(seq_len(nrow(ls)), \(i) list(ls$raw[i], ls$cell[i], round(ls$T[i], digits),
                                                round(ls[[lo]][i], digits), round(ls[[hi]][i], digits),
                                                as.integer(ls$extrapolation[i])))
      stats::setNames(
        list(ls$scale[1], s$label, s$level, s$metric, s$items, unname(groups),
             c(w("raw"), w("group"), nn, paste0(nn, "_", w("lo")), paste0(nn, "_", w("hi")), "extrapolation"),
             rows),
        c("id", "label", "level", "metric", "items", w("groups"), w("columns"), w("rows")))
    })
    stats::setNames(
      list(v, w("bayesian"), paste(title, inst$versions[[v]]$label), citation, nn, digits,
           stats::setNames(list(w("limits"), level, json_interval_text(language, interval, level, nn)),
                           c(w("type"), w("level"), w("description"))),
           unname(scales)),
      c("id", w("type"), w("title"), w("citation"), w("score"), w("score_digits"), w("interval"), w("scales")))
  })
  out <- stats::setNames(
    list(inst$id, format(Sys.time(), "%Y-%m-%d"),
         c(map_chr(cv, \(cn) wv(src_name(cn))), w("overall")), unname(sources)),
    c("instrument", w("created"), w("grouping"), w("sources")))
  jsonlite::write_json(out, file, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  invisible(file)
}

#' @rdname long_table_scale
#' @export
export_portal_json <- function(long, inst, file, citation = "", title = NULL,
                               interval = c("integrated", "conditional", "band"),
                               language = "de")
  export_json(long, inst, file, citation = citation, title = title, interval = interval, language = language)
