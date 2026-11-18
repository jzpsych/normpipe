# ==============================================================================
# run.R -- the driver
#
# norming_context() collects the inputs of a run, run_stage() runs one stage,
# run_norming() runs several. Every result is written to out_dir and cached
# with a dependency key (cache.R): a stage is skipped for a scale when its
# result exists and nothing it depends on has changed.
#
#   prepare                scoring of the reference sample, raking, routes
#   calibrate              item models (reference sample by default)
#   transfer               transfer models and projection (optional)
#   norms                  norms per reporting cell and overall, with bands
#   continuous_norms       kernel norms along a continuous variable
#   intervals_conditional  score-conditional intervals
#   intervals_integrated   integrated intervals
#   export                 result tables
# ==============================================================================

# Bump an entry when the computation of that stage changes in a way that
# makes results of earlier package versions outdated.
CACHE_SCHEMA <- c(prepare = 1L, calibrate = 1L, frame = 1L, norms = 1L, continuous_norms = 1L,
                  intervals_conditional = 1L, intervals_integrated = 1L)

STAGES <- c("prepare", "calibrate", "transfer", "norms", "continuous_norms",
            "intervals_conditional", "intervals_integrated", "export")

# Stage names of normpipe < 0.4, and the shorthand for both interval stages
STAGE_ALIASES <- list(P0 = "prepare", M1 = "calibrate", M2 = "transfer", M3a = "norms",
                      M3b = "continuous_norms", M3c = "intervals_conditional",
                      M3d = "intervals_integrated", X1 = "export",
                      intervals = c("intervals_conditional", "intervals_integrated"))

#' Stage names in their current form, in the order of the pipeline
#' @noRd
resolve_stages <- function(stages) {
  stages <- unlist(lapply(as.character(stages), \(s) STAGE_ALIASES[[s]] %||% s))
  bad <- setdiff(stages, STAGES)
  if (length(bad))
    stop("Unknown stage '", paste(bad, collapse = "', '"), "'; stages are ",
         paste(STAGES, collapse = ", "), " (and `intervals` for both interval stages).")
  intersect(STAGES, stages)
}

#' Run the norming pipeline
#'
#' `run_norming()` runs the stages of the pipeline for one instrument and
#' writes every result to `out_dir`. `norming_context()` only collects the
#' inputs of a run; `run_stage()` then runs single stages on it, which is
#' convenient for long runs and for inspecting intermediate results.
#'
#' Stages:
#' * `prepare`: scoring of the reference sample, optional raking to
#'   `margins`, imputation of skipped items (`missing` block), route per
#'   scale (direct, frequency or transfer). This stage is cheap and is run by
#'   every other stage when it has not been run on the context yet.
#' * `calibrate`: calibration of the item models ([fit_scale_model()]) in the
#'   calibration sample (default: the reference sample).
#' * `transfer`: transfer models for scales whose items the reference sample
#'   did not answer (only with `transfer: enabled: true` and `donors`, see
#'   [fit_transfer()]).
#' * `norms`: norms per reporting cell and overall, with norm bands
#'   ([norms_all_cells()]).
#' * `continuous_norms`: kernel norms along the continuous variable of
#'   `groups$continuous`, e.g. age ([weighted_norms_continuous()]).
#' * `intervals_conditional`: score-conditional intervals
#'   ([conditional_tables()]).
#' * `intervals_integrated`: integrated intervals ([integrated_intervals()]).
#'   `intervals` stands for both interval stages.
#' * `export`: the result tables (see below).
#'
#' Results. The `export` stage writes to `out_dir`, depending on
#' `export$formats` of the instrument: `norm_tables.csv` (one row per
#' version, scale, reporting cell and raw score) and
#' `norm_tables_continuous.csv` (format `long`), `compact/norms_<version>.csv`
#' ([export_compact()], format `compact`) and `json/<instrument>.json`
#' ([export_json()], format `json`, needs jsonlite), and `reliability.csv`.
#' Other stages write `scale_registry.csv`, `reference_scores.rds`,
#' `rake_diagnostics.rds`, `transfer_family_selection.csv` and
#' `transfer_crosscheck.csv`.
#'
#' Scales without an item model (`model: none`, or frequency input without
#' `models`) get norms and norm bands but no individual intervals. When
#' transfer is disabled, scales that would need it are dropped.
#'
#' Caching. Intermediate results are stored per scale as `.rds` files with a
#' dependency key ([cache_key()]) in `out_dir/cache`. A stage is skipped for a scale when its
#' result exists and neither the configuration it uses, nor the data, nor the
#' results it builds on have changed; otherwise the result is recomputed and
#' the outdated file is moved to `out_dir/cache/superseded/`. Changing, say, the
#' number of bootstrap draws therefore recomputes the frames and norm tables
#' but neither the item models nor the transfer fits.
#'
#' @param inst Instrument ([read_instrument()]).
#' @param ref Reference sample: `list(meta, items)`. `meta` has `id`, the
#'   variables of the norm groups ([group_sources()], e.g. `age` and
#'   `gender`) and optionally `weight` (unit weights otherwise); `items` holds
#'   the item columns named `<prefix><item>`, in the same row order as `meta`.
#' @param out_dir Output folder: result tables, and the cache in `cache/`.
#' @param freq Optional named list of frequency tables, one per outcome (see
#'   [frequency_frame()]); these outcomes take the frequency route.
#' @param calib Optional calibration sample `list(meta, items)` for stage
#'   `calibrate`; default `ref`.
#' @param margins Optional population margins for raking: a margin data
#'   frame ([make_margins()]) or a list of margin data frames
#'   ([rake_to_margins()]). Trimming, re-raking and the handling of
#'   respondents who cannot be raked follow the `weighting` block of the
#'   instrument.
#' @param models Optional named list of pre-fitted item models per outcome
#'   ([fit_scale_model()]); they take precedence over stage `calibrate`.
#' @param donors Donor samples for the transfer module: `list(meta, items)`
#'   with a column `meta$sample` naming the sample.
#' @param stages Stages to run (see Details); they are run in the order of
#'   the pipeline. The stage names of normpipe < 0.4 (`P0`, `M1`, `M2`,
#'   `M3a`, `M3b`, `M3c`, `M3d`, `X1`) are still understood.
#' @param outcomes Restrict to these outcomes (`<version>_<scale>`).
#' @param citation Citation string written into the JSON export.
#' @param portal_citation Name of `citation` in normpipe < 0.4.
#' @param verbose Print progress messages.
#' @return `run_norming()`: with stage `export`, the long norm table
#'   (invisibly); otherwise `NULL` (invisibly). `norming_context()`: an
#'   environment of class `normpipe_run`. `run_stage()`: the context
#'   (invisibly); after `export` the long table is in `ctx$long`. All results are written to `out_dir`.
#' @examples
#' \dontrun{
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' sim <- simulate_reference(inst, N = 1500)
#' res <- run_norming(inst, ref = list(meta = sim$meta, items = sim$items),
#'                    out_dir = "demo-out",
#'                    stages = c("calibrate", "norms", "continuous_norms",
#'                               "intervals", "export"))
#' # the same, stage by stage
#' ctx <- norming_context(inst, ref = sim, out_dir = "demo-out")
#' run_stage(ctx, "norms")
#' cache_status("demo-out")
#' }
#' @export
run_norming <- function(inst, ref = NULL, out_dir, freq = NULL, calib = NULL,
                        margins = NULL, models = NULL, donors = NULL,
                        stages = STAGES,
                        outcomes = NULL, citation = "", verbose = TRUE, portal_citation = NULL) {
  ctx <- norming_context(inst, ref = ref, out_dir = out_dir, freq = freq, calib = calib,
                         margins = margins, models = models, donors = donors,
                         outcomes = outcomes, citation = portal_citation %||% citation,
                         verbose = verbose)
  for (st in resolve_stages(stages)) run_stage(ctx, st)
  invisible(ctx$long)
}

#' @rdname run_norming
#' @export
norming_context <- function(inst, ref = NULL, out_dir, freq = NULL, calib = NULL,
                            margins = NULL, models = NULL, donors = NULL,
                            outcomes = NULL, citation = "", verbose = TRUE, portal_citation = NULL) {
  stopifnot(inherits(inst, "normpipe_instrument"))
  if (is.null(donors) && !is.null(inst$donors_data)) {
    donors <- inst$donors_data                 # donors attached to the instrument
    inst$donors_data <- NULL
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ctx <- new.env(parent = emptyenv())
  ctx$inst <- inst; ctx$ref <- ref; ctx$out_dir <- out_dir; ctx$freq <- freq
  ctx$calib <- calib %||% ref; ctx$margins <- margins; ctx$models <- models %||% list()
  ctx$donors <- donors; ctx$outcomes <- outcomes; ctx$citation <- portal_citation %||% citation
  ctx$say <- if (verbose) function(...) message(...) else function(...) invisible(NULL)
  ctx$fitted <- ctx$models
  ctx$done <- character(0)
  structure(ctx, class = "normpipe_run")
}

#' @rdname run_norming
#' @param ctx Context of a run (`norming_context()`).
#' @param stage One stage (see Details), or `intervals`.
#' @export
run_stage <- function(ctx, stage) {
  stopifnot(inherits(ctx, "normpipe_run"), length(stage) == 1)
  if (!"prepare" %in% ctx$done) { stage_prepare(ctx); ctx$done <- c(ctx$done, "prepare") }
  for (st in setdiff(resolve_stages(stage), "prepare")) {
    switch(st, calibrate = stage_calibrate(ctx), transfer = stage_transfer(ctx), norms = stage_norms(ctx),
           continuous_norms = stage_continuous_norms(ctx),
           intervals_conditional = stage_intervals_conditional(ctx),
           intervals_integrated = stage_intervals_integrated(ctx), export = stage_export(ctx))
    ctx$done <- union(ctx$done, st)
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

rds <- function(ctx, name) cache_file(ctx$out_dir, name)
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

# ---- prepare: reference scores, raking, imputation, routes --------------------

stage_prepare <- function(ctx) {
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
        say("prepare: ", sum(out), " respondent(s) answered less than ", round(100 * mc$min_answered),
            " % of the items and are not in the sample")
        ref <- list(meta = ref$meta[!out, , drop = FALSE], items = ref$items[!out, , drop = FALSE])
        ctx$ref <- ref
      }
    }
    ctx$key_p0 <- cache_key(CACHE_SCHEMA[["prepare"]], ref$meta, ref$items, ctx$margins, inst$weighting,
                            mc, inst$groups, inst$versions, as.data.frame(inst$items))
    p0 <- get_cached(ctx, "prepare_reference", ctx$key_p0, function() {
      meta <- as_tibble(ref$meta)
      if (!"weight" %in% names(meta)) meta$weight <- 1
      for (v in setdiff(frame_vars(inst), names(meta)))
        meta[[v]] <- if (v %in% categorical_sources(inst)) NA_character_ else NA_real_
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
          say("prepare: ", sum(n_imp), " skipped item response(s) imputed in ", length(imp),
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
    saveRDS(p0$scores, file.path(ctx$out_dir, "reference_scores.rds"))
    if (!is.null(p0$rake_diagnostics)) saveRDS(p0$rake_diagnostics, file.path(ctx$out_dir, "rake_diagnostics.rds"))
    item_cols <- intersect(inst$items$col, names(ref$items))
    observed_items <- inst$items$item[match(item_cols, inst$items$col)]
  }
  routes <- scale_routes(inst, observed_items)
  routes$route[routes$outcome %in% names(ctx$freq %||% list())] <- "frequency"
  if (!isTRUE(inst$transfer$enabled))
    routes <- routes[routes$route != "transfer" | routes$outcome %in% names(ctx$models), ]
  reg <- inner_join(reg, select(routes, "outcome", "route"), by = "outcome")
  write_csv(reg, file.path(ctx$out_dir, "scale_registry.csv"))
  ctx$reg <- reg
  say("Scales: ", nrow(reg), " (", paste(names(table(reg$route)), table(reg$route), collapse = ", "), ")")
  invisible(ctx)
}

# ---- calibrate: item models -----------------------------------------------------

stage_calibrate <- function(ctx) {
  inst <- ctx$inst; calib <- ctx$calib
  if (is.null(calib)) return(invisible(ctx))
  idx_all <- NULL
  for (o in ctx$reg$outcome) {
    s <- get_scale(inst, o)
    if (s$model %in% c("none", "composite") || !is.null(ctx$models[[o]])) next
    if (!is.null(s$extra$calibrate_from)) next           # scored from another scale's model
    cols <- inst$items$col[match(s$items, inst$items$item)]
    if (!all(cols %in% names(calib$items))) { ctx$say("calibrate ", o, ": items not in the calibration data, skipped"); next }
    key <- cache_key(CACHE_SCHEMA[["calibrate"]], scale_key(inst, o), inst$measurement, as.data.frame(calib$items[cols]))
    ctx$fitted[[o]] <- get_cached(ctx, paste0("calibrate_model_", o), key, function() {
      if (is.null(idx_all)) idx_all <<- item_indices(calib$items, inst, items = inst$items$item[inst$items$col %in% names(calib$items)])
      idx <- idx_all[cols]
      # full-information likelihood: skipped items are missing, not excluded
      idx <- idx[rowSums(!is.na(idx)) > 0, , drop = FALSE]
      ctx$say("calibrate ", o, ": fitting ", s$model, " on ", nrow(idx), " respondents (",
              sum(!stats::complete.cases(idx)), " with skipped items)")
      fit_scale_model(idx, s, inst)
    })
  }
  invisible(ctx)
}

#' Item model(s) of an outcome: fitted in this run, supplied, or cached by an
#' earlier run of stage `calibrate` (NULL when there is none)
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
  if (is.null(ctx$fitted[[src]]) && file.exists(rds(ctx, paste0("calibrate_model_", src))))
    ctx$fitted[[src]] <- readRDS(rds(ctx, paste0("calibrate_model_", src)))
  ctx$fitted[[src]]
}

blocks_for <- function(ctx, o, mods) {
  inst <- ctx$inst; s <- get_scale(inst, o)
  if (s$model == "composite") {
    # loadings of the total come from the calibrated component models, also
    # when `mods` holds a parameter draw of them
    point <- models_for(ctx, o)
    key <- cache_key(CACHE_SCHEMA[["calibrate"]], scale_key(inst, o), point,
                     as.data.frame(ctx$calib$items[intersect(inst$items$col[match(s$items, inst$items$item)], names(ctx$calib$items))]))
    lam <- get_cached(ctx, paste0("calibrate_loadings_", o), key, function() {
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

# ---- transfer (optional) -------------------------------------------------------

stage_transfer <- function(ctx) {
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
  cache_key(CACHE_SCHEMA[["frame"]], route, src, scale_key(inst, o), inst$groups,
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
    stop("No transfer model for '", o, "': run stage `transfer` first (needs `donors`).")
  frame_long(get_cached(ctx, paste0("frame_", o), frame_key(ctx, o), function() {
    fr <- switch(route,
      direct = direct_frame(bind_cols(select(ctx$ref_scores, all_of(c("id", frame_vars(inst), "weight"))),
                                      score = ctx$ref_scores[[o]]),
                            s, inst, n_boot = inst$uncertainty$n_boot,
                            scores_imp = if (!is.null(ctx$scores_imp)) sapply(ctx$scores_imp, \(x) x[[o]])),
      frequency = frequency_frame(ctx$freq[[o]], s, inst, n_boot = inst$uncertainty$n_boot),
      transfer = transfer_frame(ctx, o))
    frame_compact(fr)
  }))
}

# ---- norms per reporting cell, and the cross-check of the projection -------------

norms_key <- function(ctx, o, suffix = "")
  cache_key(CACHE_SCHEMA[["norms"]], frame_key(ctx, o, suffix), ctx$inst$groups, ctx$inst$norm_scale,
            ctx$inst$uncertainty)

stage_norms <- function(ctx) {
  inst <- ctx$inst
  for (o in ctx$reg$outcome) {
    s <- get_scale(inst, o)
    get_cached(ctx, paste0("norms_", o), norms_key(ctx, o), function() {
      nm <- norms_all_cells(draws_for(ctx, o), s, inst)
      ctx$say("norms ", o, " done"); nm
    })
  }
  cross <- intersect(as.character(unlist(inst$transfer$crosscheck)), ctx$reg$outcome[ctx$reg$route != "transfer"])
  if (isTRUE(inst$transfer$enabled) && length(cross)) {
    cc <- map_dfr(cross, \(o) {
      if (!file.exists(rds(ctx, paste0("transfer_", o, "_crosscheck")))) return(NULL)
      s <- get_scale(inst, o)
      proj <- get_cached(ctx, paste0("norms_crosscheck_", o), norms_key(ctx, o, "_crosscheck"), function() {
        nm <- norms_all_cells(transfer_frame(ctx, o, suffix = "_crosscheck"), s, inst)
        ctx$say("crosscheck ", o, ": norms of the projection done"); nm
      })
      crosscheck_table(readRDS(rds(ctx, paste0("norms_", o))), proj, inst) %>% mutate(outcome = o, .before = 1)
    })
    if (nrow(cc)) write_csv(cc, file.path(ctx$out_dir, "transfer_crosscheck.csv"))
  }
  invisible(ctx)
}

# ---- continuous_norms: kernel norms along the continuous variable ----------------

stage_continuous_norms <- function(ctx) {
  inst <- ctx$inst
  cont <- inst$groups$continuous
  if (is.null(cont)) return(invisible(ctx))
  for (o in ctx$reg$outcome) {
    s <- get_scale(inst, o)
    key <- cache_key(CACHE_SCHEMA[["continuous_norms"]], frame_key(ctx, o), inst$groups, inst$norm_scale,
                     inst$uncertainty)
    get_cached(ctx, paste0("continuous_norms_", o), key, function() {
      fr <- draws_for(ctx, o)
      if (all(is.na(fr[[cont$variable]]))) return(NULL)
      out <- weighted_norms_continuous(fr, s, inst)
      ctx$say("continuous_norms ", o, " done"); out
    })
  }
  invisible(ctx)
}

# ---- intervals_conditional ---------------------------------------------------------

stage_intervals_conditional <- function(ctx) {
  inst <- ctx$inst
  for (o in ctx$reg$outcome) {
    mods <- models_for(ctx, o)
    if (is.null(mods)) { ctx$say("intervals_conditional ", o, ": no item model, skipped"); next }
    s <- get_scale(inst, o)
    key <- cache_key(CACHE_SCHEMA[["intervals_conditional"]], norms_key(ctx, o), mods, inst$uncertainty)
    get_cached(ctx, paste0("intervals_conditional_", o), key, function() {
      nm <- get_cached(ctx, paste0("norms_", o), norms_key(ctx, o),
                       function() norms_all_cells(draws_for(ctx, o), s, inst))
      out <- conditional_tables(blocks_for(ctx, o, mods), s, draws_for(ctx, o), nm, inst)
      ctx$say("intervals_conditional ", o, " done"); out
    })
  }
  invisible(ctx)
}

# ---- intervals_integrated -----------------------------------------------------------

stage_intervals_integrated <- function(ctx) {
  inst <- ctx$inst
  for (o in ctx$reg$outcome) {
    mods <- models_for(ctx, o)
    if (is.null(mods)) next
    s <- get_scale(inst, o)
    key <- cache_key(CACHE_SCHEMA[["intervals_integrated"]], frame_key(ctx, o), mods, inst$uncertainty,
                     inst$groups, inst$norm_scale)
    get_cached(ctx, paste0("intervals_integrated_", o), key, function() {
      n_param <- inst$uncertainty$n_param_draws
      draws <- if (s$model == "composite") {
        comp <- map(mods, \(m) draw_model_params(m, n_param))
        map(seq_len(n_param), \(r) map(comp, \(d) d[[r]]))
      } else draw_model_params(mods, n_param)
      t0 <- Sys.time()
      out <- integrated_intervals(\(m) blocks_for(ctx, o, m), draws, s, draws_for(ctx, o), inst)
      ctx$say("intervals_integrated ", o, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min"); out
    })
  }
  invisible(ctx)
}

# ---- export -------------------------------------------------------------------------

stage_export <- function(ctx) {
  inst <- ctx$inst; reg <- ctx$reg; out_dir <- ctx$out_dir
  formats <- inst$export$formats %||% c("long", "compact", "json")
  rd <- function(p) if (file.exists(rds(ctx, p))) readRDS(rds(ctx, p)) else NULL
  missing_norms <- reg$outcome[!file.exists(rds(ctx, paste0("norms_", reg$outcome)))]
  if (length(missing_norms)) stop("No norms (stage `norms`) for: ", paste(missing_norms, collapse = ", "))
  long <- map_dfr(reg$outcome, \(o) {
    s <- get_scale(inst, o)
    long_table_scale(rd(paste0("norms_", o)), rd(paste0("intervals_conditional_", o)),
                     rd(paste0("intervals_integrated_", o)),
                     s, inst, route = reg$route[reg$outcome == o])
  })
  cont <- map_dfr(reg$outcome, \(o) {
    x <- rd(paste0("continuous_norms_", o))
    if (is.null(x)) return(NULL)
    s <- get_scale(inst, o)
    x %>% mutate(version = s$version, scale = s$id, level = s$level)
  })
  if ("long" %in% formats) {
    write_csv(long, file.path(out_dir, "norm_tables.csv"), na = "")
    if (nrow(cont)) write_csv(cont, file.path(out_dir, "norm_tables_continuous.csv"), na = "")
  }
  if ("compact" %in% formats) export_compact(long, cont, file.path(out_dir, "compact"), inst)
  if ("json" %in% formats && requireNamespace("jsonlite", quietly = TRUE)) {
    dir.create(file.path(out_dir, "json"), showWarnings = FALSE)
    export_json(long, inst, file.path(out_dir, "json", paste0(inst$id, ".json")), citation = ctx$citation)
  }
  rel <- map_dfr(reg$outcome, \(o) {
    e <- rd(paste0("intervals_conditional_", o))
    if (is.null(e$overall)) return(NULL)
    tibble(outcome = o, reliability = marginal_reliability(e$overall))
  })
  if (nrow(rel)) write_csv(rel, file.path(out_dir, "reliability.csv"))
  ctx$say("Export: ", nrow(long), " rows, ", n_distinct(paste(long$version, long$scale)), " scales")
  ctx$long <- long
  invisible(ctx)
}
