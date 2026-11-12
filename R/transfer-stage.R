# ==============================================================================
# transfer-stage.R -- stage M2 of run_norming(): donor frame, prediction
# frames, transfer fits (with the choice among candidate families)
# ==============================================================================

#' Cache of stage M2
#'
#' Uses `cache_get()` (fingerprinted cache of normpipe >= 0.3.0) when
#' available, otherwise a plain file cache: an existing `<name>.rds` in
#' `out_dir` is reused.
#' @noRd
transfer_cache <- function(out_dir, say = message) {
  if (exists("cache_get", mode = "function", inherits = TRUE))
    return(function(name, key, compute) cache_get(out_dir, name, key, compute, say = say))
  function(name, key, compute) {
    f <- file.path(out_dir, paste0(name, ".rds"))
    if (file.exists(f)) return(readRDS(f))
    x <- compute()
    saveRDS(x, f)
    x
  }
}

#' Fingerprint of a transfer model: everything that changes the fit
#' @noRd
transfer_key <- function(td, backend) {
  vars <- intersect(union(all.vars(td$formula$formula), unlist(lapply(td$formula$pforms, all.vars))),
                    names(td$data))
  rlang::hash(list(outcome = td$outcome, family = td$family,
                   formula = format(td$formula$formula), pforms = lapply(td$formula$pforms, format),
                   prior = as.data.frame(td$prior), data = as.data.frame(td$data[vars]),
                   spec = td$spec[c("adapt_delta", "max_treedepth", "chains", "iter", "init")],
                   backend = backend))
}

#' Stage M2 of run_norming(): transfer models
#'
#' Builds the donor frame and the prediction frame(s) of the reference sample,
#' then fits a transfer model for every scale on the transfer route and for
#' every cross-check outcome (`<outcome>_lodo`). With several candidate
#' families (`transfer$family: [beta_binomial, gaussian]`), each candidate is
#' fitted (`transfer_<key>__<family>`) and the best converged one by PSIS-LOO
#' becomes `transfer_<key>`; the choice is written to
#' `m2_family_selection.csv`. Dense scales keep the censored candidates
#' (the beta-binomial and the cumulative model need a count structure).
#'
#' @param inst Instrument.
#' @param reg Scale registry (stage P0).
#' @param ref Reference data (`meta`, `items`).
#' @param ref_scores Scores of the reference sample (`id`, `age`, `gender`, `weight`).
#' @param donors Donor samples: `list(meta, items)`.
#' @param out_dir Output directory.
#' @param say Message function.
#' @param ref_imp Optional list of completed item data sets of the reference
#'   sample ([impute_items()]); one prediction frame each.
#' @param cache Cache function `(name, key, compute)`; default
#'   [transfer_cache()].
#' @noRd
run_transfer_stage <- function(inst, reg, ref, ref_scores, donors = inst$donors_data, out_dir,
                               say = message, ref_imp = NULL, cache = transfer_cache(out_dir, say)) {
  if (is.null(donors)) stop("Transfer is enabled but no donor samples were supplied ",
                            "(donors = list(meta, items)).")
  check_transfer_config(inst)
  need <- reg$outcome[reg$route == "transfer"]
  cross <- intersect(as.character(unlist(inst$transfer$crosscheck)), reg$outcome[reg$route != "transfer"])
  if (!length(need) && !length(cross)) return(invisible(NULL))
  # donors: one inclusion rule and completed data sets, as for the reference sample
  mc <- missing_cfg(inst)
  donor_frames <- NULL
  if (mc$method == "impute") {
    share <- answered_share(donors$items, inst)
    out <- is.na(share) | share < mc$min_answered
    if (any(out)) {
      say("M2: ", sum(out), " donor(s) answered less than ", round(100 * mc$min_answered),
          " % of the items and are not used")
      donors <- list(meta = donors$meta[!out, , drop = FALSE], items = donors$items[!out, , drop = FALSE])
    }
    cv <- intersect(mc$covariates, names(donors$meta))
    n_sets <- max(transfer_spec_defaults()$chains, inst$transfer$chains %||% 0)
    imp <- impute_items(donors$items, inst, M = n_sets, seed = sum(utf8ToInt(paste0(inst$id, "_donors"))),
                        covariates = if (length(cv)) as.data.frame(donors$meta)[cv], k = mc$neighbours)
    if (sum(attr(imp, "n_imputed")) > 0) {
      say("M2: ", sum(attr(imp, "n_imputed")), " skipped item response(s) of the donors imputed; ",
          "each chain of a transfer model uses one of ", n_sets, " completed data sets")
      donor_frames <- map(imp, \(it) transfer_donor_frame(list(meta = donors$meta, items = it), inst, ref_scores))
    }
  }
  donor_frame <- donor_frames[[1]] %||% transfer_donor_frame(donors, inst, ref_scores)
  saveRDS(donor_frame, file.path(out_dir, "m2_donor_frame.rds"))
  n_na <- sum(is.na(ref_scores$age) | is.na(ref_scores$gender) |
                !as.character(ref_scores$gender) %in% levels(donor_frame$gender))
  if (n_na) say("M2: ", n_na, " reference respondent(s) without age or with a gender ",
                "category absent from the donor samples are not projected")

  backend <- transfer_backend(inst)
  todo <- c(set_names(need, need), set_names(cross, paste0(cross, "_lodo")))
  sel_all <- list()
  for (i in seq_along(todo)) {
    o <- todo[[i]]; k <- names(todo)[i]
    s <- get_scale(inst, o)
    sp <- transfer_spec(inst, o)
    if (endsWith(k, "_lodo") && !length(sp$exclude)) sp$exclude <- own_bridge(inst, o)
    fams <- sp$family
    # count families need the count structure (not for dense scales); if no
    # other candidate is left, bb_scale() reports the problem
    if (any(count_families %in% fams) && !bb_possible(s) && length(setdiff(fams, count_families)))
      fams <- setdiff(fams, count_families)
    excl <- if (length(sp$exclude)) paste0(", without ", paste(sp$exclude, collapse = "/")) else ""
    tds <- set_names(map(fams, \(fam) transfer_data(o, donor_frame, inst, family = fam, spec = sp, s = s)), fams)
    # with imputed donors: one transfer_data per completed data set and chain
    tds_imp <- if (!is.null(donor_frames)) set_names(map(fams, \(fam)
      map(donor_frames[seq_len(min(sp$chains, length(donor_frames)))],
          \(d) transfer_data(o, d, inst, family = fam, spec = sp, s = s))), fams)
    key_of <- function(fam) if (is.null(tds_imp)) transfer_key(tds[[fam]], backend) else
      rlang::hash(map_chr(tds_imp[[fam]], transfer_key, backend = backend))
    fit_one <- function(fam, name) {
      cache(name, key_of(fam), function() {
        say("M2 ", k, ": fitting ", fam, excl)
        if (is.null(tds_imp)) transfer_model(tds[[fam]], backend = backend) else
          transfer_model_multiple(tds_imp[[fam]], backend = backend)
      })
    }
    if (length(fams) == 1) {
      fit_one(fams, paste0("transfer_", k))
      next
    }
    fits <- set_names(map(fams, \(fam) fit_one(fam, paste0("transfer_", k, "__", fam))), fams)
    keys <- map_chr(fams, key_of)
    # "phi_plateau": selection rule of normpipe 0.3.1 (recomputed for cached fits)
    sel <- cache(paste0("transfer_", k, "_selection"),
                 rlang::hash(list(keys, sp$select_rhat, sp$select_rhat_chain, "phi_plateau")), function() {
      x <- select_transfer_family(fits, rhat_max = sp$select_rhat, rhat_chain = sp$select_rhat_chain %||% 1.05)
      x$outcome <- k
      x
    })
    best <- sel$family[sel$chosen]
    cache(paste0("transfer_", k), rlang::hash(list(keys, sp$select_rhat, best)), function() {
      why <- ifelse(sel$eligible, "", sprintf(", not eligible: R-hat %.3f, %d at max. tree depth%s",
                                               ifelse(is.na(sel$rhat_within), sel$rhat_max, sel$rhat_within),
                                               sel$treedepth_hits,
                                               ifelse(sel$phi_plateau > 0, sprintf(", %.0f%% of draws at phi > 1e10",
                                                                                   100 * sel$phi_plateau), "")))
      say("M2 ", k, ": ", best, " chosen by LOO (",
          paste(sprintf("%s %.1f%s", sel$family, sel$elpd_loo, why), collapse = "; "), ")",
          if (!any(sel$eligible)) " -- no candidate converged, chosen among all; check the fits" else "")
      fits[[best]]
    })
    sel_all[[k]] <- sel
  }
  if (length(sel_all))
    readr::write_csv(bind_rows(sel_all) %>% select("outcome", everything()),
                     file.path(out_dir, "m2_family_selection.csv"))
  invisible(NULL)
}
