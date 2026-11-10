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
#' `m2_family_selection.csv`. Dense scales keep the censored candidates.
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
  donor_frame <- transfer_donor_frame(donors, inst, ref_scores)
  saveRDS(donor_frame, file.path(out_dir, "m2_donor_frame.rds"))
  pred_of <- function(items) transfer_prediction_frame(items, ref_scores, ref$meta$id, donor_frame, inst)
  pred <- pred_of(ref$items)
  f_imp <- file.path(out_dir, "m2_prediction_frames_imputed.rds")
  if (!is.null(ref_imp)) saveRDS(map(ref_imp, pred_of), f_imp) else if (file.exists(f_imp)) file.remove(f_imp)
  n_na <- sum(is.na(pred$age) | is.na(pred$gender))
  if (n_na) say("M2: ", n_na, " reference respondent(s) without age or with a gender ",
                "category absent from the donor samples are not projected")
  saveRDS(pred, file.path(out_dir, "m2_prediction_frame.rds"))

  backend <- transfer_backend(inst)
  todo <- c(set_names(need, need), set_names(cross, paste0(cross, "_lodo")))
  sel_all <- list()
  for (i in seq_along(todo)) {
    o <- todo[[i]]; k <- names(todo)[i]
    s <- get_scale(inst, o)
    sp <- transfer_spec(inst, o)
    if (endsWith(k, "_lodo") && !length(sp$exclude)) sp$exclude <- own_bridge(inst, o)
    fams <- sp$family
    if ("beta_binomial" %in% fams && length(fams) > 1 && !bb_possible(s))
      fams <- setdiff(fams, "beta_binomial")
    excl <- if (length(sp$exclude)) paste0(", without ", paste(sp$exclude, collapse = "/")) else ""
    tds <- set_names(map(fams, \(fam) transfer_data(o, donor_frame, inst, family = fam, spec = sp, s = s)), fams)
    fit_one <- function(fam, name) {
      key <- transfer_key(tds[[fam]], backend)
      cache(name, key, function() {
        say("M2 ", k, ": fitting ", fam, excl)
        transfer_model(tds[[fam]], backend = backend)
      })
    }
    if (length(fams) == 1) {
      fit_one(fams, paste0("transfer_", k))
      next
    }
    fits <- set_names(map(fams, \(fam) fit_one(fam, paste0("transfer_", k, "__", fam))), fams)
    keys <- map_chr(fams, \(fam) transfer_key(tds[[fam]], backend))
    sel <- cache(paste0("transfer_", k, "_selection"), rlang::hash(list(keys, sp$select_rhat)), function() {
      x <- select_transfer_family(fits, rhat_max = sp$select_rhat)
      x$outcome <- k
      x
    })
    best <- sel$family[sel$chosen]
    cache(paste0("transfer_", k), rlang::hash(list(keys, sp$select_rhat, best)), function() {
      why <- ifelse(sel$eligible, "", sprintf(", not converged: R-hat %.3f, %d at max. tree depth",
                                               sel$rhat_max, sel$treedepth_hits))
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
