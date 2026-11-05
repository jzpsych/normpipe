# ==============================================================================
# export.R -- the published tables
#
# norm_tables_long.csv   one row per version, scale, reporting cell and raw
#                        score: percentile, T, norm band (T_lo, T_hi), score-
#                        conditional interval (T_cond_lo/hi), integrated
#                        interval (T_int_lo/hi), extrapolation flag
# norm_tables_age_continuous.csv   continuous-age norms (kernel), per year
# platform/norms_<version>.csv     compact per-version tables for scoring software
# portal/<instrument>.json         norm file in the JSON format of a web portal
# ==============================================================================

#' Export formats
#'
#' `long_table_scale()` assembles the long norm table of one scale: one row
#' per reporting cell and raw score with percentile rank, T, norm band,
#' score-conditional and integrated interval (interpolated to every raw
#' score) and an extrapolation flag (median percentile rank below 0.5 % or
#' above 99.5 %). [run_norming()] binds these tables to
#' `norm_tables_long.csv`.
#'
#' `export_platform_tables()` writes compact per-version tables for scoring
#' software: `norms_<version>.csv` (raw score, T, integrated interval,
#' percentile rank, extrapolation flag; cells `gender`/`age_group`, `"all"`
#' for the overall table) and, with continuous-age norms,
#' `norms_<version>_age.rds`.
#'
#' `export_portal_json()` writes one JSON norm file per instrument for a web
#' portal: one source per version, the reporting cells as groups, and per row
#' raw score, cell, T and interval limits. The keys and the default labels
#' of this format are German (e.g. `quellen`, `gruppen`, `zeilen`).
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
                          "extrapolation", "n_eff")))
}

#' @rdname long_table_scale
#' @param long Long norm table (rows of `long_table_scale()`).
#' @param cont Continuous-age norms ([weighted_norms_continuous()] with
#'   `version`, `scale`, `level`) or `NULL`.
#' @param dir Output folder.
#' @return `export_platform_tables()`: paths of the written files
#'   (invisibly).
#' @export
export_platform_tables <- function(long, cont, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  band <- long %>%
    transmute(.data$version, .data$level, .data$scale,
              gender = if ("gender" %in% names(long)) coalesce(as.character(.data$gender), "all") else "all",
              age_group = if ("age_group" %in% names(long)) coalesce(as.character(.data$age_group), "all") else "all",
              raw = round(.data$raw, 4), T = round(.data$T, 1),
              T_int_lo = round(.data$T_int_lo, 1), T_int_hi = round(.data$T_int_hi, 1),
              pctl_med = round(.data$pctl_med, 4),
              extrapolation = as.integer(.data$extrapolation)) %>%
    arrange(.data$version, .data$level, .data$scale, .data$gender, .data$age_group, .data$raw)
  for (v in unique(band$version))
    write_csv(filter(band, .data$version == v), file.path(dir, paste0("norms_", v, ".csv")), na = "")
  if (!is.null(cont) && nrow(cont)) {
    cc <- cont %>%
      transmute(.data$version, .data$level, .data$scale, gender = coalesce(as.character(.data$gender), "all"),
                age = as.integer(.data$age), raw = round(.data$raw, 4), T = round(.data$T, 1),
                pctl_med = round(.data$pctl_med, 4),
                extrapolation = as.integer(.data$pctl_med < .005 | .data$pctl_med > .995),
                .data$half_width) %>%
      arrange(.data$version, .data$level, .data$scale, .data$gender, .data$age, .data$raw)
    for (v in unique(cc$version))
      saveRDS(as.data.frame(filter(cc, .data$version == v)),
              file.path(dir, paste0("norms_", v, "_age.rds")), compress = "xz")
  }
  invisible(list.files(dir, full.names = TRUE))
}

#' @rdname long_table_scale
#' @param file Output file.
#' @param citation Citation string of the norms.
#' @param title Title prefix of every source (followed by the version label).
#' @param interval Interval written per row: `"integrated"`, `"conditional"`
#'   or `"band"`.
#' @return `export_portal_json()`: `file` (invisibly).
#' @export
export_portal_json <- function(long, inst, file, citation = "", title = "Bayesianische Normen",
                               interval = c("integrated", "conditional", "band")) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is not installed.")
  interval <- match.arg(interval)
  lo <- switch(interval, integrated = "T_int_lo", conditional = "T_cond_lo", band = "T_lo")
  hi <- switch(interval, integrated = "T_int_hi", conditional = "T_cond_hi", band = "T_hi")
  cv <- cell_vars(inst)
  quellen <- map(split(long, long$version), \(lv) {
    v <- lv$version[1]
    scales <- map(split(lv, lv$scale), \(ls) {
      s <- get_scale(inst, paste0(v, "_", ls$scale[1]))
      cells <- split(ls, ls$cell)
      gruppen <- map(cells, \(cd) {
        cell <- cd$cell[1]
        g <- list(id = gsub("[^A-Za-z0-9]+", "_", cell),
                  label = if (cell == "overall") "Gesamt" else cell,
                  n_eff = round(cd$n_eff[1]))
        if ("gender" %in% cv && cell != "overall") g$geschlecht <- as.character(cd$gender[1])
        if ("age_group" %in% cv && cell != "overall") g$alter <- as.character(cd$age_group[1])
        g
      })
      zeilen <- map(seq_len(nrow(ls)), \(i) list(ls$raw[i], ls$cell[i], round(ls$T[i], 1),
                                                   round(ls[[lo]][i], 1), round(ls[[hi]][i], 1),
                                                   as.integer(ls$extrapolation[i])))
      list(id = ls$scale[1], label = s$label, level = s$level, metric = s$metric,
           items = s$items, gruppen = unname(gruppen),
           spalten = c("roh", "gruppe", "T", "T_unten", "T_oben", "extrapolation"),
           zeilen = zeilen)
    })
    list(id = v, art = "bayesianisch", titel = paste(title, inst$versions[[v]]$label),
         zitation = citation, kennwert = "T", dezimalen_kennwert = 1,
         intervall = list(art = "grenzen", niveau = 95,
                          erklaerung = switch(interval,
                            integrated = "95-%-Intervall f\u00fcr den wahren T-Wert der Person, einschlie\u00dflich der Unsicherheit der Normwerte",
                            conditional = "95-%-Intervall f\u00fcr den wahren T-Wert der Person bei festen Normwerten",
                            band = "95-%-Glaubw\u00fcrdigkeitsband des Normwerts")),
         skalen = unname(scales))
  })
  out <- list(instrument = inst$id, erzeugt = format(Sys.time(), "%Y-%m-%d"),
              bezugsgruppen = c(if ("gender" %in% cv) "geschlecht", if ("age_group" %in% cv) "alter", "gesamt"),
              quellen = unname(quellen))
  jsonlite::write_json(out, file, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  invisible(file)
}
