# ==============================================================================
# run.R -- the driver
#
# norming_context() collects the inputs of a run, run_stage() runs one stage,
# run_norming() runs several. Every result is written to out_dir and cached
# with a dependency key (cache.R): a stage is skipped for a scale when its
# result exists and nothing it depends on has changed.
#
#   P0   scoring of the reference sample, raking, routes per scale
#   M1   calibration of the item models (reference sample by default)
#   M2   transfer models and projection (optional; transfer.R)
#   M3a  band norms per reporting cell and overall, with norm bands
#   M3b  continuous-age kernel norms
#   M3c  score-conditional EAP tables
#   M3d  integrated intervals
#   X1   exports (long table, per-version tables, portal JSON)
# ==============================================================================

# Bump an entry when the computation of that stage changes in a way that
# makes results of earlier package versions outdated.
CACHE_SCHEMA <- c(p0 = 1L, m1 = 1L, frame = 1L, m3a = 1L, m3b = 1L, m3c = 1L, m3d = 1L)

STAGES <- c("P0", "M1", "M2", "M3a", "M3b", "M3c", "M3d", "X1")

#' Run the norming pipeline
#'
#' `run_norming()` runs the stages of the pipeline for one instrument and
#' writes every result to `out_dir`. `norming_context()` only collects the
#' inputs of a run; `run_stage()` then runs single stages on it, which is
#' convenient for long runs and for inspecting intermediate results.
#'
#' Stages:
#' * `P0`: scoring of the reference sample, optional raking to `margins`,
#'   imputation of skipped items (`missing` block), route per scale (direct,
#'   frequency or transfer). `P0` is cheap and is run by every other stage
#'   when it has not been run on the context yet.
#' * `M1`: calibration of the item models ([fit_scale_model()]) in the
#'   calibration sample (default: the reference sample).
#' * `M2`: transfer models for scales whose items the reference sample did
#'   not answer (only with `transfer: enabled: true` and `donors`, see
#'   [fit_transfer()]).
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
#' Caching. Intermediate results are stored per scale as `.rds` files with a
#' dependency key ([cache_key()]). A stage is skipped for a scale when its
#' result exists and neither the configuration it uses, nor the data, nor the
#' results it builds on have changed; otherwise the result is recomputed and
#' the outdated file is moved to `out_dir/superseded/`. Changing, say, the
#' number of bootstrap draws therefore recomputes the frames and norm tables
#' but neither the item models nor the transfer fits.
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
#'   instrument.
#' @param models Optional named list of pre-fitted item models per outcome
#'   ([fit_scale_model()]); they take precedence over stage `M1`.
#' @param donors Donor samples for the transfer module: `list(meta, items)`
#'   with a column `meta$sample` naming the sample.
#' @param stages Stages to run (see Details).
#' @param outcomes Restrict to these outcomes (`<version>_<scale>`).
#' @param portal_citation Citation string written into the portal JSON.
#' @param verbose Print progress messages.
#' @return `run_norming()`: with stage `X1`, the long norm table (invisibly);
#'   otherwise `NULL` (invisibly). `norming_context()`: an environment of
#'   class `normpipe_run`. `run_stage()`: the context (invisibly); for `X1`
#'   the long table is in `ctx$long`. All results are written to `out_dir`.
#' @examples
#' \dontrun{
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' sim <- simulate_reference(inst, N = 1500)
#' res <- run_norming(inst, ref = list(meta = sim$meta, items = sim$items),
#'                    out_dir = "demo-out",
#'                    stages = c("P0", "M1", "M3a", "M3b", "M3c", "M3d", "X1"))
#' # the same, stage by stage
#' ctx <- norming_context(inst, ref = sim, out_dir = "demo-out")
#' run_stage(ctx, "M3a")
#' cache_status("demo-out")
#' }
#' @export
run_norming <- function(inst, ref = NULL, out_dir, freq = NULL, calib = NULL,
                        margins = NULL, models = NULL, donors = NULL,
                        stages = STAGES,
                        outcomes = NULL, portal_citation = "", verbose = TRUE) {
  ctx <- norming_context(inst, ref = ref, out_dir = out_dir, freq = freq, calib = calib,
                         margins = margins, models = models, donors = donors,
                         outcomes = outcomes, portal_citation = portal_citation,
                         verbose = verbose)
  for (st in intersect(STAGES, stages)) run_stage(ctx, st)
  invisible(ctx$long)
}

#' @rdname run_norming
#' @export
norming_context <- function(inst, ref = NULL, out_dir, freq = NULL, calib = NULL,
                            margins = NULL, models = NULL, donors = NULL,
                            outcomes = NULL, portal_citation = "", verbose = TRUE) {
  stopifnot(inherits(inst, "normpipe_instrument"))
  if (is.null(donors) && !is.null(inst$donors_data)) {
    donors <- inst$donors_data                 # accepted until 0.2.x scripts are updated
    inst$donors_data <- NULL
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ctx <- new.env(parent = emptyenv())
  ctx$inst <- inst; ctx$ref <- ref; ctx$out_dir <- out_dir; ctx$freq <- freq
  ctx$calib <- calib %||% ref; ctx$margins <- margins; ctx$models <- models %||% list()
  ctx$donors <- donors; ctx$outcomes <- outcomes; ctx$portal_citation <- portal_citation
  ctx$say <- if (verbose) function(...) message(...) else function(...) invisible(NULL)
  ctx$fitted <- ctx$models
  ctx$done <- character(0)
  structure(ctx, class = "normpipe_run")
}

#' @rdname run_norming
#' @param ctx Context of a run (`norming_context()`).
#' @param stage One of the stages (see Details).
#' @export
run_stage <- function(ctx, stage) {
  stopifnot(inherits(ctx, "normpipe_run"), length(stage) == 1)
  if (!stage %in% STAGES) stop("Unknown stage '", stage, "'; stages are ", paste(STAGES, collapse = ", "), ".")
  if (!"P0" %in% ctx$done) { stage_p0(ctx); ctx$done <- c(ctx$done, "P0") }
  if (stage != "P0") {
    switch(stage, M1 = stage_m1(ctx), M2 = stage_m2(ctx), M3a = stage_m3a(ctx),
           M3b = stage_m3b(ctx), M3c = stage_m3c(ctx), M3d = stage_m3d(ctx), X1 = stage_x1(ctx))
    ctx$done <- union(ctx$done, stage)
  }
  invisible(ctx)
}

#' @export
print.normpipe_run <- function(x, ...) {
  cat("normpipe run of '", x$inst$id, "' in ", x$out_dir, "\n", sep = "")
  cat("  stages run on this context: ", if (length(x$done)) paste(x$done, collapse = ", ") else "none", "\n", sep = "")
  if (!is.null(x$reg)) cat("  scales: ", nrow(x$reg), " (", paste(names(table(x$reg$route)), table(x$reg$route), collapse = ", "), ")\n", sep = "")
  invisible(x)
}

# ---- helpers -----------------------------------------------------------------

rds <- function(ctx, name) file.path(ctx$out_dir, paste0(name, ".rds"))
get_cached <- function(ctx, name, key, compute) cache_get(ctx$out_dir, name, key, compute, say = ctx$say)

#' Key of the definition of a scale (scale entry plus its items)
#' @noRd
scale_key <- function(inst, o) {
  s <- get_scale(inst, o)
  cache_key(s, as.data.frame(inst$items[inst$items$item %in% s$items, ]))
}

#' Key of a result some other part of the package wrote (transfer fits):
#' the manifest key when there is one, the file size and time otherwise.
#' @noRd
file_key <- function(ctx, name) {
  man <- read_manifest(ctx$out_dir)
  hit <- match(name, man$name)
  if (!is.na(hit)) return(man$key[hit])
  fi <- file.info(rds(ctx, name))
  cache_key(name, fi$size, as.numeric(fi$mtime))
}

# ---- P0: reference scores, raking, imputation, routes -------------------------

stage_p0 <- function(ctx) {
  inst <- ctx$inst; ref <- ctx$ref; say <- ctx$say
  reg <- scale_registry(inst)
  if (!is.null(ctx$outcomes)) reg <- reg[reg$outcome %in% ctx$outcomes, ]
  observed_items <- integer(0)
  ctx$ref_scores <- NULL; ctx$ref_imp <- NULL; ctx$scores_imp <- NULL
  mc <- missing_cfg(inst)
  if (!is.null(ref)) {
    if (mc$method == "impute") {
      # unit non-response: too few of the presented items answered
      share <- answered_share(ref$items, inst)
      out <- is.na(share) | share < mc$min_answered
      if (any(out)) {
        say("P0: ", sum(out), " respondent(s) answered less than ", round(100 * mc$min_answered),
            " % of the items and are not in the sample")
        ref <- list(meta = ref$meta[!out, , drop = FALSE], items = ref$items[!out, , drop = FALSE])
        ctx$ref <- ref
      }
    }
    ctx$key_p0 <- cache_key(CACHE_SCHEMA[["p0"]], ref$meta, ref$items, ctx$margins, inst$weighting,
                            mc, inst$cells, inst$versions, as.data.frame(inst$items))
    p0 <- get_cached(ctx, "p0_reference", ctx$key_p0, function() {
      meta <- as_tibble(ref$meta)
      if (!"weight" %in% names(meta)) meta$weight <- 1
      if (!"gender" %in% names(meta)) meta$gender <- NA_character_
      if (!"age" %in% names(meta)) meta$age <- NA_real_
      dg <- NULL
      if (!is.null(ctx$margins)) {
        wcfg <- inst$weighting %||% list()
        raked <- rake_to_margins(meta, ctx$margins, inst,
                                 trim_quantile = wcfg$trim_quantile %||% 0.99,
                                 rerake = wcfg$rerake %||% TRUE,
                                 incomplete = wcfg$incomplete %||% "drop")
        dg <- attr(raked, "diagnostics")
        say("Raking: Kish N_eff ", round(dg$kish_neff_design), " -> ", round(dg$kish_neff_raked))
        meta <- meta %>% select(-"weight") %>%
          inner_join(select(raked, "id", weight = "weight_raked"), by = "id")
      }
      rows <- match(meta$id, ref$meta$id)
      scores <- bind_cols(meta, score_scales(ref$items, inst, partial = TRUE)[rows, ])
      imp <- NULL; scores_imp <- NULL
      if (mc$method == "impute") {
        cv <- intersect(mc$covariates, names(ref$meta))
        imp <- impute_items(ref$items, inst, M = mc$n_imputations, seed = sum(utf8ToInt(inst$id)),
                            covariates = if (length(cv)) as.data.frame(ref$meta)[cv], k = mc$neighbours)
        n_imp <- attr(imp, "n_imputed")
        if (sum(n_imp) == 0) imp <- NULL else {
          say("P0: ", sum(n_imp), " skipped item response(s) imputed in ", length(imp),
              " completed data sets")
          scores_imp <- map(imp, \(d) score_scales(d, inst, partial = TRUE)[rows, ])
          # everyone in the sample has a score on every presented scale; the
          # plain score table shows the first completed data set
          scores <- bind_cols(meta, scores_imp[[1]])
        }
      }
      list(scores = scores, imputed = imp, scores_imp = scores_imp, rake_diagnostics = dg)
    })
    ctx$ref_scores <- p0$scores; ctx$ref_imp <- p0$imputed; ctx$scores_imp <- p0$scores_imp
    # plain copies for inspection and for scripts that read them
    saveRDS(p0$scores, rds(ctx, "p0_reference_scores"))
    if (!is.null(p0$rake_diagnostics)) saveRDS(p0$rake_diagnostics, rds(ctx, "p0_rake_diagnostics"))
    item_cols <- intersect(inst$items$col, names(ref$items))
    observed_items <- inst$items$item[match(item_cols, inst$items$col)]
  }
  routes <- scale_routes(inst, observed_items)
  routes$route[routes$outcome %in% names(ctx$freq %||% list())] <- "frequency"
  if (!isTRUE(inst$transfer$enabled))
    routes <- routes[routes$route != "transfer" | routes$outcome %in% names(ctx$models), ]
  reg <- inner_join(reg, select(routes, "outcome", "route"), by = "outcome")
  write_csv(reg, file.path(ctx$out_dir, "p0_scale_registry.csv"))
  ctx$reg <- reg
  say("Scales: ", nrow(reg), " (", paste(names(table(reg$route)), table(reg$route), collapse = ", "), ")")
  invisible(ctx)
}

# ---- M1: item models ------------------------------------------------------------

stage_m1 <- function(ctx) {
  inst <- ctx$inst; calib <- ctx$calib
  if (is.null(calib)) return(invisible(ctx))
  idx_all <- NULL
  for (o in ctx$reg$outcome) {
    s <- get_scale(inst, o)
    if (s$model %in% c("none", "composite") || !is.null(ctx$models[[o]])) next
    if (!is.null(s$extra$calibrate_from)) next           # scored from another scale's model
    cols <- inst$items$col[match(s$items, inst$items$item)]
    if (!all(cols %in% names(calib$items))) { ctx$say("M1 ", o, ": items not in the calibration data, skipped"); next }
    key <- cache_key(CACHE_SCHEMA[["m1"]], scale_key(inst, o), inst$measurement, as.data.frame(calib$items[cols]))
    ctx$fitted[[o]] <- get_cached(ctx, paste0("m1_model_", o), key, function() {
      if (is.null(idx_all)) idx_all <<- item_indices(calib$items, inst, items = inst$items$item[inst$items$col %in% names(calib$items)])
      idx <- idx_all[cols]
      # full-information likelihood: skipped items are missing, not excluded
      idx <- idx[rowSums(!is.na(idx)) > 0, , drop = FALSE]
      ctx$say("M1 ", o, ": fitting ", s$model, " on ", nrow(idx), " respondents (",
              sum(!stats::complete.cases(idx)), " with skipped items)")
      fit_scale_model(idx, s, inst)
    })
  }
  invisible(ctx)
}

#' Item model(s) of an outcome: fitted in this run, supplied, or cached by an
#' earlier run of stage M1 (NULL when there is none)
#' @noRd
models_for <- function(ctx, o) {
  s <- get_scale(ctx$inst, o)
  if (s$model == "none") return(NULL)
  if (s$model == "composite") {
    ms <- map(set_names(names(s$components)), \(cp) models_for(ctx, paste0(s$version, "_", cp)))
    if (any(map_lgl(ms, is.null))) return(NULL)
    return(ms)
  }
  src <- s$extra$calibrate_from %||% o
  if (is.null(ctx$fitted[[src]]) && file.exists(rds(ctx, paste0("m1_model_", src))))
    ctx$fitted[[src]] <- readRDS(rds(ctx, paste0("m1_model_", src)))
  ctx$fitted[[src]]
}

blocks_for <- function(ctx, o, mods) {
  inst <- ctx$inst; s <- get_scale(inst, o)
  if (s$model == "composite") {
    # loadings of the total come from the calibrated component models, also
    # when `mods` holds a parameter draw of them
    point <- models_for(ctx, o)
    key <- cache_key(CACHE_SCHEMA[["m1"]], scale_key(inst, o), point,
                     as.data.frame(ctx$calib$items[intersect(inst$items$col[match(s$items, inst$items$item)], names(ctx$calib$items))]))
    lam <- get_cached(ctx, paste0("m1_loadings_", o), key, function() {
      idx <- item_indices(ctx$calib$items, inst, items = s$items)
      total_loadings(idx, s, inst, map(point, \(m) m))
    })
    # component models may be calibrated on longer scales: subset by item id
    return(make_blocks(s, mods, lambda = lam))
  }
  subset <- match(s$items, mods$items)
  stopifnot(!anyNA(subset))
  make_blocks(s, mods, subset = subset)
}

# ---- M2: transfer (optional) -----------------------------------------------------

stage_m2 <- function(ctx) {
  inst <- ctx$inst
  if (!isTRUE(inst$transfer$enabled)) return(invisible(ctx))
  run_transfer_stage(inst, ctx$reg, ctx$ref, ctx$ref_scores, donors = ctx$donors,
                     out_dir = ctx$out_dir, say = ctx$say, ref_imp = ctx$ref_imp)
  invisible(ctx)
}

# ---- frames per outcome ----------------------------------------------------------

frame_key <- function(ctx, o, suffix = "") {
  inst <- ctx$inst
  route <- if (nzchar(suffix)) "transfer" else ctx$reg$route[ctx$reg$outcome == o]
  src <- switch(route,
    direct = ctx$key_p0,
    frequency = cache_key(ctx$freq[[o]]),
    transfer = cache_key(file_key(ctx, paste0("transfer_", o, suffix)), ctx$key_p0,
                         inst$transfer[c("n_draws", "projection_level", "bootstrap_weights", "bridge")]))
  cache_key(CACHE_SCHEMA[["frame"]], route, src, scale_key(inst, o), inst$cells,
            inst$uncertainty[c("n_boot", "prior_mass")])
}

#' Draw frame of a transfer outcome: the cached fit projected into the
#' reference sample (one prediction frame per completed data set)
#' @noRd
transfer_frame <- function(ctx, o, suffix = "") {
  inst <- ctx$inst; s <- get_scale(inst, o)
  fit <- readRDS(rds(ctx, paste0("transfer_", o, suffix)))
  frames <- prediction_frames(ctx$ref_scores, ctx$ref$meta$id, ctx$ref_imp %||% list(ctx$ref$items), fit, inst)
  set.seed(sum(utf8ToInt(o)))
  project_outcome(fit, frames, s, n_draws = inst$transfer$n_draws %||% 500,
                  level = inst$transfer$projection_level %||% NULL,
                  bootstrap_weights = inst$transfer$bootstrap_weights %||% TRUE)
}

draws_for <- function(ctx, o) {
  inst <- ctx$inst; s <- get_scale(inst, o)
  route <- ctx$reg$route[ctx$reg$outcome == o]
  if (route == "transfer" && !file.exists(rds(ctx, paste0("transfer_", o))))
    stop("No transfer model for '", o, "': run stage M2 first (needs `donors`).")
  frame_long(get_cached(ctx, paste0("frame_", o), frame_key(ctx, o), function() {
    fr <- switch(route,
      direct = direct_frame(ctx$ref_scores %>% transmute(.data$id, .data$age, .data$gender,
                                                         .data$weight, score = .data[[o]]),
                            s, inst, n_boot = inst$uncertainty$n_boot,
                            scores_imp = if (!is.null(ctx$scores_imp)) sapply(ctx$scores_imp, \(x) x[[o]])),
      frequency = frequency_frame(ctx$freq[[o]], s, inst, n_boot = inst$uncertainty$n_boot),
      transfer = transfer_frame(ctx, o))
    frame_compact(fr)
  }))
}

# ---- M3a: band norms, and the leave-own-domain-out cross-check --------------------

m3a_key <- function(ctx, o, suffix = "")
  cache_key(CACHE_SCHEMA[["m3a"]], frame_key(ctx, o, suffix), ctx$inst$cells, ctx$inst$uncertainty)

stage_m3a <- function(ctx) {
  inst <- ctx$inst
  for (o in ctx$reg$outcome) {
    s <- get_scale(inst, o)
    get_cached(ctx, paste0("m3a_norms_", o), m3a_key(ctx, o), function() {
      nm <- norms_all_cells(draws_for(ctx, o), s, inst)
      ctx$say("M3a ", o, " done"); nm
    })
  }
  cross <- intersect(as.character(unlist(inst$transfer$crosscheck)), ctx$reg$outcome[ctx$reg$route != "transfer"])
  if (isTRUE(inst$transfer$enabled) && length(cross)) {
    cc <- map_dfr(cross, \(o) {
      if (!file.exists(rds(ctx, paste0("transfer_", o, "_lodo")))) return(NULL)
      s <- get_scale(inst, o)
      lodo <- get_cached(ctx, paste0("m3a_lodo_", o), m3a_key(ctx, o, "_lodo"), function() {
        nm <- norms_all_cells(transfer_frame(ctx, o, suffix = "_lodo"), s, inst)
        ctx$say("C1 ", o, ": leave-own-domain-out norms done"); nm
      })
      crosscheck_table(readRDS(rds(ctx, paste0("m3a_norms_", o))), lodo) %>% mutate(outcome = o, .before = 1)
    })
    if (nrow(cc)) write_csv(cc, file.path(ctx$out_dir, "c1_crosscheck.csv"))
  }
  invisible(ctx)
}

# ---- M3b: continuous-age norms ----------------------------------------------------

stage_m3b <- function(ctx) {
  inst <- ctx$inst
  if (is.null(inst$cells$ages)) return(invisible(ctx))
  for (o in ctx$reg$outcome) {
    s <- get_scale(inst, o)
    key <- cache_key(CACHE_SCHEMA[["m3b"]], frame_key(ctx, o), inst$cells, inst$uncertainty)
    get_cached(ctx, paste0("m3b_continuous_", o), key, function() {
      fr <- draws_for(ctx, o)
      if (all(is.na(fr$age))) return(NULL)
      out <- weighted_norms_continuous(fr, s, inst)
      ctx$say("M3b ", o, " done"); out
    })
  }
  invisible(ctx)
}

# ---- M3c: conditional tables --------------------------------------------------------

stage_m3c <- function(ctx) {
  inst <- ctx$inst
  for (o in ctx$reg$outcome) {
    mods <- models_for(ctx, o)
    if (is.null(mods)) { ctx$say("M3c ", o, ": no item model, skipped"); next }
    s <- get_scale(inst, o)
    key <- cache_key(CACHE_SCHEMA[["m3c"]], m3a_key(ctx, o), mods, inst$uncertainty)
    get_cached(ctx, paste0("m3c_eap_", o), key, function() {
      nm <- get_cached(ctx, paste0("m3a_norms_", o), m3a_key(ctx, o),
                       function() norms_all_cells(draws_for(ctx, o), s, inst))
      out <- conditional_tables(blocks_for(ctx, o, mods), s, draws_for(ctx, o), nm, inst)
      ctx$say("M3c ", o, " done"); out
    })
  }
  invisible(ctx)
}

# ---- M3d: integrated intervals --------------------------------------------------------

stage_m3d <- function(ctx) {
  inst <- ctx$inst
  for (o in ctx$reg$outcome) {
    mods <- models_for(ctx, o)
    if (is.null(mods)) next
    s <- get_scale(inst, o)
    key <- cache_key(CACHE_SCHEMA[["m3d"]], frame_key(ctx, o), mods, inst$uncertainty, inst$cells)
    get_cached(ctx, paste0("m3d_integrated_", o), key, function() {
      n_param <- inst$uncertainty$n_param_draws
      draws <- if (s$model == "composite") {
        comp <- map(mods, \(m) draw_model_params(m, n_param))
        map(seq_len(n_param), \(r) map(comp, \(d) d[[r]]))
      } else draw_model_params(mods, n_param)
      t0 <- Sys.time()
      out <- integrated_intervals(\(m) blocks_for(ctx, o, m), draws, s, draws_for(ctx, o), inst)
      ctx$say("M3d ", o, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min"); out
    })
  }
  invisible(ctx)
}

# ---- X1: exports ------------------------------------------------------------------------

stage_x1 <- function(ctx) {
  inst <- ctx$inst; reg <- ctx$reg; out_dir <- ctx$out_dir
  rd <- function(p) if (file.exists(rds(ctx, p))) readRDS(rds(ctx, p)) else NULL
  missing_norms <- reg$outcome[!file.exists(rds(ctx, paste0("m3a_norms_", reg$outcome)))]
  if (length(missing_norms)) stop("No norms (stage M3a) for: ", paste(missing_norms, collapse = ", "))
  long <- map_dfr(reg$outcome, \(o) {
    s <- get_scale(inst, o)
    long_table_scale(rd(paste0("m3a_norms_", o)), rd(paste0("m3c_eap_", o)), rd(paste0("m3d_integrated_", o)),
                     s, inst, route = reg$route[reg$outcome == o])
  })
  write_csv(long, file.path(out_dir, "norm_tables_long.csv"), na = "")
  cont <- map_dfr(reg$outcome, \(o) {
    x <- rd(paste0("m3b_continuous_", o))
    if (is.null(x)) return(NULL)
    s <- get_scale(inst, o)
    x %>% mutate(version = s$version, scale = s$id, level = s$level)
  })
  if (nrow(cont)) write_csv(cont, file.path(out_dir, "norm_tables_age_continuous.csv"), na = "")
  export_platform_tables(long, cont, file.path(out_dir, "platform"))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    dir.create(file.path(out_dir, "portal"), showWarnings = FALSE)
    export_portal_json(long, inst, file.path(out_dir, "portal", paste0(inst$id, ".json")),
                       citation = ctx$portal_citation)
  }
  rel <- map_dfr(reg$outcome, \(o) {
    e <- rd(paste0("m3c_eap_", o))
    if (is.null(e$overall)) return(NULL)
    tibble(outcome = o, reliability = marginal_reliability(e$overall))
  })
  if (nrow(rel)) write_csv(rel, file.path(out_dir, "x1_reliability.csv"))
  ctx$say("Export: ", nrow(long), " rows, ", n_distinct(paste(long$version, long$scale)), " scales")
  ctx$long <- long
  invisible(ctx)
}
