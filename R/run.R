# ==============================================================================
# run.R -- the driver
#
# run_norming() runs the stages of the pipeline for one instrument and writes
# every result to out_dir; stages are cached per scale and skipped when their
# file exists.
#
#   P0   scoring of the reference sample, raking, routes per scale
#   M1   calibration of the item models (reference sample by default)
#   M2   transfer models and projection (optional; transfer.R)
#   M3a  band norms per reporting cell and overall, with norm bands
#   M3b  continuous-age kernel norms
#   M3c  score-conditional EAP tables
#   M3d  integrated intervals
#   X1   exports (long table, per-version tables, portal JSON)
#
# INPUT. `ref` is a list(meta, items): meta has id, age, gender, weight (weight
# may be missing: unit weights), items the item columns named as in the
# instrument's item table. Alternatively `freq` is a named list of frequency
# tables (one per outcome, see frequency_frame()) for scales without item
# data; those scales get norms and norm bands but no individual intervals
# unless item parameters are supplied through `models`.
# ==============================================================================

#' Run the norming pipeline
#'
#' Runs the stages of the pipeline for one instrument and writes every result
#' to `out_dir`. Intermediate results are cached per scale as `.rds` files and
#' a stage is skipped for a scale when its cache file exists; delete the
#' files (or use a new `out_dir`) to recompute.
#'
#' Stages:
#' * `P0`: scoring of the reference sample, optional raking to `margins`,
#'   route per scale (direct, frequency or transfer).
#' * `M1`: calibration of the item models ([fit_scale_model()]) in the
#'   calibration sample (default: the reference sample).
#' * `M2`: transfer models and projection for scales whose items the
#'   reference sample did not answer (only with `transfer: enabled: true`;
#'   donor samples in `inst$donors_data`, see [fit_transfer()]).
#' * `M3a`: band norms per reporting cell and overall ([norms_all_cells()]).
#' * `M3b`: continuous-age kernel norms ([weighted_norms_continuous()]).
#' * `M3c`: score-conditional intervals ([conditional_tables()]).
#' * `M3d`: integrated intervals ([integrated_intervals()]).
#' * `X1`: exports: `norm_tables_long.csv`, `norm_tables_age_continuous.csv`,
#'   `platform/norms_<version>.csv` ([export_platform_tables()]),
#'   `portal/<instrument>.json` ([export_portal_json()], needs jsonlite) and
#'   `x1_reliability.csv`.
#'
#' Scales without an item model (`model: none`, or frequency input without
#' `models`) get norms and norm bands but no individual intervals. When
#' transfer is disabled, scales that would need it are dropped.
#'
#' @param inst Instrument ([read_instrument()]).
#' @param ref Reference sample: `list(meta, items)`. `meta` has `id`, `age`,
#'   `gender` and optionally `weight` (unit weights otherwise); `items` holds
#'   the item columns named `<prefix><item>`, in the same row order as `meta`.
#' @param out_dir Output folder (caches and tables).
#' @param freq Optional named list of frequency tables, one per outcome (see
#'   [frequency_frame()]); these outcomes take the frequency route.
#' @param calib Optional calibration sample `list(meta, items)` for stage
#'   `M1`; default `ref`.
#' @param margins Optional population margins for raking: a margin data
#'   frame ([make_margins()]) or a list of margin data frames
#'   ([rake_to_margins()]). Trimming, re-raking and the handling of
#'   respondents who cannot be raked follow the `weighting` block of the
#'   instrument (`trim_quantile` [0.99], `rerake` [true], `incomplete`
#'   [drop]).
#' @param models Optional named list of pre-fitted item models per outcome
#'   ([fit_scale_model()]); they take precedence over stage `M1`.
#' @param stages Stages to run (see Details).
#' @param outcomes Restrict to these outcomes (`<version>_<scale>`).
#' @param portal_citation Citation string written into the portal JSON.
#' @param verbose Print progress messages.
#' @return With stage `X1`, the long norm table (invisibly); otherwise
#'   `NULL` (invisibly). All results are written to `out_dir`.
#' @examples
#' \dontrun{
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' sim <- simulate_reference(inst, N = 1500)
#' res <- run_norming(inst, ref = list(meta = sim$meta, items = sim$items),
#'                    out_dir = "demo-out",
#'                    stages = c("P0", "M1", "M3a", "M3b", "M3c", "M3d", "X1"))
#' }
#' @export
run_norming <- function(inst, ref = NULL, out_dir, freq = NULL, calib = NULL,
                        margins = NULL, models = NULL,
                        stages = c("P0", "M1", "M2", "M3a", "M3b", "M3c", "M3d", "X1"),
                        outcomes = NULL, portal_citation = "", verbose = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  say <- function(...) if (verbose) message(...)
  cache <- function(name) file.path(out_dir, paste0(name, ".rds"))
  reg <- scale_registry(inst)
  if (!is.null(outcomes)) reg <- reg[reg$outcome %in% outcomes, ]

  # ---- P0: reference scores, raking, routes -------------------------------
  observed_items <- integer(0)
  ref_scores <- NULL; ref_imp <- NULL; scores_imp <- NULL
  if (!is.null(ref)) {
    meta <- as_tibble(ref$meta)
    if (!"weight" %in% names(meta)) meta$weight <- 1
    if (!"gender" %in% names(meta)) meta$gender <- NA_character_
    if (!"age" %in% names(meta)) meta$age <- NA_real_
    if (!is.null(margins)) {
      wcfg <- inst$weighting %||% list()
      raked <- rake_to_margins(meta, margins, inst,
                               trim_quantile = wcfg$trim_quantile %||% 0.99,
                               rerake = wcfg$rerake %||% TRUE,
                               incomplete = wcfg$incomplete %||% "drop")
      saveRDS(attr(raked, "diagnostics"), cache("p0_rake_diagnostics"))
      dg <- attr(raked, "diagnostics")
      say("Raking: Kish N_eff ", round(dg$kish_neff_design), " -> ", round(dg$kish_neff_raked))
      meta <- meta %>% select(-"weight") %>%
        inner_join(select(raked, "id", weight = "weight_raked"), by = "id")
    }
    item_cols <- intersect(inst$items$col, names(ref$items))
    observed_items <- inst$items$item[match(item_cols, inst$items$col)]
    sc <- score_scales(ref$items, inst, partial = TRUE)
    ref_scores <- bind_cols(meta, sc[match(meta$id, ref$meta$id), ])
    saveRDS(ref_scores, cache("p0_reference_scores"))
    # skipped items: completed data sets (missing$reference: impute)
    if (!inst$missing$reference %in% c("prorate", "impute"))
      stop("missing$reference must be 'prorate' or 'impute'.")
    if (identical(inst$missing$reference, "impute")) {
      f <- cache("p0_imputed_items")
      if (file.exists(f)) ref_imp <- readRDS(f) else {
        ref_imp <- impute_items(ref$items, inst, M = inst$missing$n_imputations,
                                seed = sum(utf8ToInt(inst$id)))
        saveRDS(ref_imp, f)
      }
      n_imp <- attr(ref_imp, "n_imputed")
      if (sum(n_imp) == 0) ref_imp <- NULL else {
        say("P0: ", sum(n_imp), " skipped item response(s) imputed in ", length(ref_imp),
            " completed data sets")
        scores_imp <- map(ref_imp, \(d) score_scales(d, inst, partial = TRUE)[match(meta$id, ref$meta$id), ])
      }
    }
  }
  routes <- scale_routes(inst, observed_items)
  routes$route[routes$outcome %in% names(freq %||% list())] <- "frequency"
  if (!isTRUE(inst$transfer$enabled))
    routes <- routes[routes$route != "transfer" | routes$outcome %in% names(models %||% list()), ]
  reg <- inner_join(reg, select(routes, "outcome", "route"), by = "outcome")
  write_csv(reg, file.path(out_dir, "p0_scale_registry.csv"))
  say("Scales: ", nrow(reg), " (", paste(names(table(reg$route)), table(reg$route), collapse = ", "), ")")

  # ---- M1: item models ------------------------------------------------------
  calib <- calib %||% ref
  fitted <- models %||% list()
  if ("M1" %in% stages && !is.null(calib)) {
    f <- cache("m1_models")
    if (file.exists(f)) fitted <- modifyList(readRDS(f), fitted)
    idx_all <- NULL
    for (o in reg$outcome) {
      s <- get_scale(inst, o)
      if (s$model %in% c("none", "composite") || !is.null(fitted[[o]])) next
      src <- s$extra$calibrate_from
      if (!is.null(src)) next           # scored from another scale's model
      cols <- inst$items$col[match(s$items, inst$items$item)]
      if (!all(cols %in% names(calib$items))) { say("M1 ", o, ": items not in the calibration data, skipped"); next }
      if (is.null(idx_all)) idx_all <- item_indices(calib$items, inst, items = inst$items$item[inst$items$col %in% names(calib$items)])
      idx <- idx_all[cols]
      idx <- idx[stats::complete.cases(idx), , drop = FALSE]
      say("M1 ", o, ": fitting ", s$model, " on ", nrow(idx), " complete cases")
      fitted[[o]] <- fit_scale_model(idx, s, inst)
    }
    saveRDS(fitted, f)
  }

  # ---- M2: transfer (optional) -----------------------------------------------
  if ("M2" %in% stages && isTRUE(inst$transfer$enabled)) {
    if (!exists("run_transfer_stage")) stop("transfer.R is not available")
    run_transfer_stage(inst, reg, ref, ref_scores, out_dir, say, ref_imp = ref_imp)
  }

  # ---- frames per outcome ---------------------------------------------------
  draws_for <- function(o, n_boot = inst$uncertainty$n_boot) {
    s <- get_scale(inst, o)
    route <- reg$route[reg$outcome == o]
    f <- cache(paste0("frame_", o))
    if (file.exists(f)) return(frame_long(readRDS(f)))
    fr <- switch(route,
      direct = direct_frame(ref_scores %>% transmute(.data$id, .data$age, .data$gender,
                                                     .data$weight, score = .data[[o]]),
                            s, inst, n_boot = n_boot,
                            scores_imp = if (!is.null(scores_imp)) sapply(scores_imp, \(x) x[[o]])),
      frequency = frequency_frame(freq[[o]], s, inst, n_boot = n_boot),
      transfer = project_frame(o, s, inst, out_dir))
    saveRDS(frame_compact(fr), f)
    fr
  }

  # ---- M3a: band norms ------------------------------------------------------
  if ("M3a" %in% stages) for (o in reg$outcome) {
    f <- cache(paste0("m3a_norms_", o))
    if (file.exists(f)) next
    s <- get_scale(inst, o)
    saveRDS(norms_all_cells(draws_for(o), s, inst), f)
    say("M3a ", o, " done")
  }

  # ---- C1: leave-own-domain-out cross-check (transfer$crosscheck) -----------
  cross <- intersect(unlist(inst$transfer$crosscheck), reg$outcome[reg$route != "transfer"])
  if ("M3a" %in% stages && isTRUE(inst$transfer$enabled) && length(cross)) {
    cc <- map_dfr(cross, \(o) {
      f <- cache(paste0("m3a_lodo_", o))
      if (!file.exists(f)) {
        if (!file.exists(file.path(out_dir, paste0("transfer_", o, "_lodo.rds")))) return(NULL)
        s <- get_scale(inst, o)
        saveRDS(norms_all_cells(project_frame(o, s, inst, out_dir, suffix = "_lodo"), s, inst), f)
        say("C1 ", o, ": leave-own-domain-out norms done")
      }
      direct <- cache(paste0("m3a_norms_", o))
      if (!file.exists(direct)) return(NULL)
      crosscheck_table(readRDS(direct), readRDS(f)) %>% mutate(outcome = o, .before = 1)
    })
    if (nrow(cc)) write_csv(cc, file.path(out_dir, "c1_crosscheck.csv"))
  }

  # ---- M3b: continuous-age norms -------------------------------------------
  if ("M3b" %in% stages && !is.null(inst$cells$ages)) for (o in reg$outcome) {
    f <- cache(paste0("m3b_continuous_", o))
    if (file.exists(f)) next
    s <- get_scale(inst, o)
    fr <- draws_for(o)
    if (all(is.na(fr$age))) next
    saveRDS(weighted_norms_continuous(fr, s, inst), f)
    say("M3b ", o, " done")
  }

  # ---- models and blocks per outcome -----------------------------------------
  models_for <- function(o) {
    s <- get_scale(inst, o)
    if (s$model == "none") return(NULL)
    if (s$model == "composite") {
      ms <- map(set_names(names(s$components)), \(cp) models_for(paste0(s$version, "_", cp)))
      if (any(map_lgl(ms, is.null))) return(NULL)
      return(ms)
    }
    src <- s$extra$calibrate_from %||% o
    fitted[[src]]
  }
  blocks_for <- function(o, mods) {
    s <- get_scale(inst, o)
    if (s$model == "composite") {
      lam_f <- cache(paste0("m1_loadings_", o))
      lam <- if (file.exists(lam_f)) readRDS(lam_f) else {
        cols <- inst$items$col[match(s$items, inst$items$item)]
        idx <- item_indices(calib$items, inst, items = s$items)
        l <- total_loadings(idx, s, inst, map(mods, \(m) m))
        saveRDS(l, lam_f); l
      }
      # component models may be calibrated on longer scales: subset by item id
      return(make_blocks(s, mods, lambda = lam))
    }
    subset <- match(s$items, mods$items)
    stopifnot(!anyNA(subset))
    make_blocks(s, mods, subset = subset)
  }

  # ---- M3c: conditional tables -----------------------------------------------
  if ("M3c" %in% stages) for (o in reg$outcome) {
    f <- cache(paste0("m3c_eap_", o))
    if (file.exists(f)) next
    mods <- models_for(o)
    if (is.null(mods)) { say("M3c ", o, ": no item model, skipped"); next }
    s <- get_scale(inst, o)
    nm <- readRDS(cache(paste0("m3a_norms_", o)))
    saveRDS(conditional_tables(blocks_for(o, mods), s, draws_for(o), nm, inst), f)
    say("M3c ", o, " done")
  }

  # ---- M3d: integrated intervals ---------------------------------------------
  if ("M3d" %in% stages) for (o in reg$outcome) {
    f <- cache(paste0("m3d_integrated_", o))
    if (file.exists(f)) next
    mods <- models_for(o)
    if (is.null(mods)) next
    s <- get_scale(inst, o)
    n_param <- inst$uncertainty$n_param_draws
    draws <- if (s$model == "composite") {
      comp <- map(mods, \(m) draw_model_params(m, n_param))
      map(seq_len(n_param), \(r) map(comp, \(d) d[[r]]))
    } else draw_model_params(mods, n_param)
    t0 <- Sys.time()
    saveRDS(integrated_intervals(\(m) blocks_for(o, m), draws, s, draws_for(o), inst), f)
    say("M3d ", o, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  }

  # ---- X1: exports --------------------------------------------------------------
  if ("X1" %in% stages) {
    long <- map_dfr(reg$outcome, \(o) {
      s <- get_scale(inst, o)
      nm <- readRDS(cache(paste0("m3a_norms_", o)))
      rd <- function(p) if (file.exists(cache(p))) readRDS(cache(p)) else NULL
      long_table_scale(nm, rd(paste0("m3c_eap_", o)), rd(paste0("m3d_integrated_", o)),
                       s, inst, route = reg$route[reg$outcome == o])
    })
    write_csv(long, file.path(out_dir, "norm_tables_long.csv"), na = "")
    cont <- map_dfr(reg$outcome, \(o) {
      f <- cache(paste0("m3b_continuous_", o))
      if (!file.exists(f)) return(NULL)
      s <- get_scale(inst, o)
      readRDS(f) %>% mutate(version = s$version, scale = s$id, level = s$level)
    })
    if (nrow(cont)) write_csv(cont, file.path(out_dir, "norm_tables_age_continuous.csv"), na = "")
    export_platform_tables(long, cont, file.path(out_dir, "platform"))
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      dir.create(file.path(out_dir, "portal"), showWarnings = FALSE)
      export_portal_json(long, inst, file.path(out_dir, "portal", paste0(inst$id, ".json")),
                         citation = portal_citation)
    }
    rel <- map_dfr(reg$outcome, \(o) {
      f <- cache(paste0("m3c_eap_", o))
      if (!file.exists(f)) return(NULL)
      e <- readRDS(f)
      if (is.null(e$overall)) return(NULL)
      tibble(outcome = o, reliability = marginal_reliability(e$overall))
    })
    if (nrow(rel)) write_csv(rel, file.path(out_dir, "x1_reliability.csv"))
    say("Export: ", nrow(long), " rows, ", n_distinct(paste(long$version, long$scale)), " scales")
    return(invisible(long))
  }
  invisible(NULL)
}
