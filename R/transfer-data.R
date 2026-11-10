# ==============================================================================
# transfer-data.R -- data of the transfer models: bridge scores, donor and
# prediction frames, and transfer_data() (response, censoring, formula,
# priors and initial values of one model; pure R, no Stan)
# ==============================================================================

#' Bridge scores
#'
#' One predictor per bridge block (mean item value of the block, raw metric),
#' their overall mean, and the scored single bridge items. Missing items are
#' ignored in the means; a block without any answered item is `NA`.
#'
#' @param items_df Data frame with the bridge item columns.
#' @param inst Instrument ([read_instrument()]) with a `transfer` block.
#' @return Tibble with one column `br_<id>` per bridge block, `br_mean`
#'   (mean over all bridge items) and `bi_<item>` (scored bridge item,
#'   integer units).
#' @export
bridge_scores <- function(items_df, inst) {
  br <- inst$transfer$bridge
  if (is.null(br) || !length(br)) stop("The instrument has no transfer$bridge blocks.")
  ids <- map_chr(br, "id")
  if (anyDuplicated(ids)) stop("Bridge block ids must be unique: ", paste(ids[duplicated(ids)], collapse = ", "))
  all_items <- unique(unlist(map(br, "items")))
  idx <- item_indices(items_df, inst, items = all_items)
  vals <- item_values(idx, inst, items = all_items)
  unit <- attr(inst$items, "unit")
  col_of <- \(it) inst$items$col[match(it, inst$items$item)]
  block_mean <- \(items) {
    m <- rowMeans(as.matrix(vals[col_of(items)]), na.rm = TRUE) * unit
    m[is.nan(m)] <- NA
    m
  }
  out <- as_tibble(setNames(map(br, \(b) block_mean(b$items)), paste0("br_", ids)))
  out$br_mean <- block_mean(all_items)
  for (it in all_items) out[[paste0("bi_", it)]] <- vals[[col_of(it)]]
  out
}

#' Donor frame of the transfer models
#'
#' Meta data (`id`, `age`, `gender`, `sample`), bridge scores and scale scores
#' (not prorated) of the donor samples. `sample` is a factor whose first level
#' is `transfer$reference_level`; gender levels follow `inst$cells$gender`.
#' With `ref_scores`, donors outside the age range of the reference sample are
#' dropped.
#'
#' @param donors List with `meta` and `items` (one completed data set).
#' @param inst Instrument.
#' @param ref_scores Optional scores of the reference sample (for the age range).
#' @return Tibble.
#' @export
transfer_donor_frame <- function(donors, inst, ref_scores = NULL) {
  dm <- as_tibble(donors$meta)
  need <- c("age", "gender", "sample")
  if (length(setdiff(need, names(dm))))
    stop("Donor meta data need the columns ", paste(need, collapse = ", "), ".")
  lv <- inst$transfer$reference_level %||% levels(factor(dm$sample))[1]
  if (!lv %in% dm$sample) stop("transfer$reference_level '", lv, "' is not a donor sample.")
  dm$sample <- factor(dm$sample, levels = c(lv, setdiff(unique(as.character(dm$sample)), lv)))
  g_don <- unique(stats::na.omit(as.character(dm$gender)))
  g_lv <- c(intersect(inst$cells$gender, g_don), sort(setdiff(g_don, inst$cells$gender)))
  dm$gender <- factor(as.character(dm$gender), levels = g_lv)
  out <- bind_cols(dm, bridge_scores(donors$items, inst),
                   score_scales(donors$items, inst, prorate = FALSE))
  if (!is.null(ref_scores))
    out <- out %>% filter(between(.data$age, min(ref_scores$age, na.rm = TRUE),
                                  max(ref_scores$age, na.rm = TRUE)))
  out
}

#' Prediction frame of the reference sample
#'
#' @param items Item data of the reference sample (one completed data set).
#' @param ref_scores Scores of the reference sample (`id`, `age`, `gender`, `weight`).
#' @param ref_ids Ids of the rows of `items`.
#' @param donor_frame Donor frame ([transfer_donor_frame()]) for the factor levels.
#' @param inst Instrument.
#' @return Tibble with `id`, `age`, `gender`, `weight`, the bridge scores and
#'   `sample` (reference level).
#' @noRd
transfer_prediction_frame <- function(items, ref_scores, ref_ids, donor_frame, inst) {
  lv <- levels(donor_frame$sample)
  bind_cols(ref_scores %>% select("id", "age", "gender", "weight"),
            bridge_scores(items, inst)[match(ref_scores$id, ref_ids), ]) %>%
    mutate(gender = factor(as.character(.data$gender), levels = levels(donor_frame$gender)),
           sample = factor(lv[1], levels = lv))
}

#' Data and model definition of one transfer model
#'
#' Builds everything brms needs for one outcome and one family, without
#' calling Stan: the response (only the part of the score that the reference
#' sample did not observe, see below), censoring or trials, the formula, the
#' brms family, the priors and the initial values.
#'
#' Bridge items that belong to the outcome are observed in the reference
#' sample and are never predicted: the censored families (`gaussian`,
#' `skew_normal`) model the outcome minus its weighted observed items,
#' interval-censored at the attainable range of that rest; the beta-binomial
#' models the sum of the unobserved items out of their range. The projection
#' ([project_outcome()]) adds the observed part again.
#'
#' @param outcome Outcome name (`<version>_<scale>`), a column of `dat`.
#' @param dat Donor frame ([transfer_donor_frame()]).
#' @param inst Instrument.
#' @param family One of `gaussian`, `skew_normal`, `beta_binomial`; default:
#'   the first family of the specification.
#' @param spec Transfer specification ([transfer_spec()]); arguments of
#'   [fit_transfer()] override its entries.
#' @param s Scale specification ([get_scale()]).
#' @return List of class `normpipe_transfer_data` with `outcome`, `family`,
#'   `data`, `formula` (brmsformula), `bfamily`, `prior`, `init` (function,
#'   0 or `NULL` for random), `info` (observed items), `preds`, `spec`.
#' @export
transfer_data <- function(outcome, dat, inst, family = NULL,
                          spec = transfer_spec(inst, outcome), s = get_scale(inst, outcome)) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("brms is not installed.")
  family <- family %||% spec$family[1]
  if (!family %in% transfer_families) stop("Unknown transfer family '", family, "'.")
  if (!outcome %in% names(dat)) stop("Outcome '", outcome, "' is not a column of the donor frame.")
  exclude <- spec$exclude
  preds <- setdiff(grep("^br_", names(dat), value = TRUE), c("br_mean", paste0("br_", exclude)))
  ok <- stats::complete.cases(dat[c(outcome, preds)])
  r2 <- summary(stats::lm(stats::reformulate(preds, outcome), data = dat[ok, ]))$r.squared
  if (r2 > 0.999)
    stop("Outcome '", outcome, "' is a (near-)exact linear combination of the ",
         "bridge predictors (R^2 = ", round(r2, 5), "); exclude its own bridge blocks.")
  known <- observed_items(s, inst, exclude)
  info <- observed_info(s, inst, known, dat)
  rhs <- expand_rhs(spec$formula, preds, dat)
  y_full <- dat[[outcome]]
  if (family == "beta_binomial") {
    bb <- bb_scale(s, known = known, values = info$range_units)
    kn <- observed_part(dat, info, wt = rep(1, length(known)))
    dat$y_sum <- as.integer(round(y_full / bb$unit_per_sum - kn)) - bb$offset
    dat$y_trials <- bb$trials
    f <- brms::bf(stats::as.formula(paste0("y_sum | trials(y_trials) ~ ", rhs)),
                  stats::as.formula(paste0("phi ~ ", expand_rhs(spec$phi, preds, dat))))
    bfam <- brms::beta_binomial()
    # All parameters start at zero on the unconstrained scale (phi = 1, mean at
    # the logit of 0.5 for centred predictors). From random values a chain can
    # start on the plateau of large phi, where the likelihood is flat in phi,
    # and never leave it; see also the proper prior on the phi coefficients.
    init <- 0
  } else {
    h <- s$step; lo <- s$range[1]; hi <- s$range[2]
    y <- y_full
    if (length(known)) {
      y <- y - observed_part(dat, info)
      lo <- info$rest[["min"]]; hi <- info$rest[["max"]]
    }
    eps <- 1e-9
    dat <- dat %>%
      mutate(cens_ind = case_when(is.na(y) ~ NA_character_,
                                  y <= lo + eps ~ "left", y >= hi - eps ~ "right", TRUE ~ "interval"),
             cens_lo = case_when(.data$cens_ind == "left" ~ lo + h / 2,
                                 .data$cens_ind == "right" ~ hi - h / 2,
                                 TRUE ~ y - h / 2),
             cens_hi = case_when(.data$cens_ind == "interval" ~ y + h / 2, TRUE ~ .data$cens_lo))
    f <- brms::bf(stats::as.formula(paste0("cens_lo | cens(cens_ind, cens_hi) ~ ", rhs)),
                  stats::as.formula(paste0("sigma ~ ", expand_rhs(spec$sigma, preds, dat))))
    bfam <- switch(family, gaussian = stats::gaussian(), skew_normal = brms::skew_normal())
    # Intercepts at the location and log scale of the outcome, every other
    # parameter at zero (lower-bounded ones at one). From random values the
    # censored likelihood of an outcome far from zero can be -Inf. The list is
    # complete: rstan fails on some models with partial initial values.
    y_mid <- (dat$cens_lo + dat$cens_hi) / 2
    init <- full_init(f, dat, bfam,
                      list(Intercept = mean(y_mid, na.rm = TRUE),
                           Intercept_sigma = log(stats::sd(y_mid, na.rm = TRUE))))
  }
  if (!is.null(spec$init)) init <- if (identical(spec$init, "random")) NULL else 0
  structure(list(outcome = outcome, family = family, data = dat, formula = f, bfamily = bfam,
                 prior = transfer_priors(spec, family), init = init, info = info,
                 preds = preds, spec = spec),
            class = "normpipe_transfer_data")
}

#' Complete initial values for a brms model
#'
#' Reads the parameters block of the Stan code (generated, not compiled) and
#' returns a function that gives every parameter a value: `values` where
#' named, otherwise zero (one for parameters with a lower bound of zero), each
#' with a small jitter.
#' @noRd
full_init <- function(f, dat, fam, values = list()) {
  code <- strsplit(brms::make_stancode(f, data = dat, family = fam), "\n")[[1]]
  sd <- brms::make_standata(f, data = dat, family = fam)
  i0 <- grep("^parameters \\{", code); i1 <- i0 + which(grepl("^\\}", code[-(1:i0)]))[1]
  decl <- trimws(sub("//.*$", "", code[(i0 + 1):(i1 - 1)]))
  decl <- decl[nzchar(decl)]
  pars <- lapply(decl, function(d) {
    m <- regmatches(d, regexec("^(real|vector)(<[^>]*>)?(\\[(.*)\\])?\\s+([A-Za-z0-9_]+);$", d))[[1]]
    if (!length(m)) stop("normpipe cannot set initial values for the declaration '", d, "'.")
    n <- if (nzchar(m[5])) as.integer(eval(parse(text = m[5]), envir = sd)) else 1L
    pos <- grepl("lower\\s*=\\s*0", m[3])
    list(name = m[6], vector = m[2] == "vector", n = n, base = if (pos) 1 else 0)
  })
  function() {
    out <- lapply(pars, function(p) {
      v <- values[[p$name]] %||% p$base
      v <- rep_len(v, p$n) + stats::rnorm(p$n, 0, 0.05)
      if (p$base == 1) v <- abs(v)
      if (p$vector) array(v, dim = p$n) else v
    })
    stats::setNames(out, map_chr(pars, "name"))
  }
}
