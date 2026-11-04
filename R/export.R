# ==============================================================================
# export.R -- the published tables
#
# norm_tables_long.csv   one row per version, scale, reporting cell and raw
#                        score: percentile, T, norm band (T_lo, T_hi), score-
#                        conditional interval (T_cond_lo/hi), integrated
#                        interval (T_int_lo/hi), extrapolation flag
# norm_tables_age_continuous.csv   continuous-age norms (kernel), per year
# platform/norms_<version>.csv     the reduced tables of the scoring platform
# portal/<instrument>.json         the norm file of the instrument portal
# ==============================================================================

#' Assemble the long table of one scale
#'
#' @param nm norms_all_cells() result, cond conditional_tables() result (or
#'   NULL), intg integrated_intervals() result (or NULL)
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

#' Platform tables: one csv per version with the columns the scoring
#' platform reads (T, integrated interval, percentile, extrapolation)
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

#' Norm file of the instrument portal (daten/normen/<id>.json): one source
#' per version with the cells as groups and the integrated interval per row
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
                            integrated = "95-%-Intervall für den wahren T-Wert der Person, einschließlich der Unsicherheit der Normwerte",
                            conditional = "95-%-Intervall für den wahren T-Wert der Person bei festen Normwerten",
                            band = "95-%-Glaubwürdigkeitsband des Normwerts")),
         skalen = unname(scales))
  })
  out <- list(instrument = inst$id, erzeugt = format(Sys.time(), "%Y-%m-%d"),
              bezugsgruppen = c(if ("gender" %in% cv) "geschlecht", if ("age_group" %in% cv) "alter", "gesamt"),
              quellen = unname(quellen))
  jsonlite::write_json(out, file, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  invisible(file)
}
