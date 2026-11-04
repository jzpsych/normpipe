# ==============================================================================
# irt.R -- item models, calibration, parameter draws, expected-score curves,
#          summed-score EAP (Supplement S3.3, S5.2, S5.3)
#
# The item models are stored in a plain representation that does not depend
# on mirt: a slope matrix a (items x factors), a list of intercept vectors d
# (one per item, decreasing; mirt's parameterization of the graded response
# model, P(X >= k | theta) = logistic(a' theta + d_k)), the category VALUES of
# every item, and a parameter table with the free parameters and their
# asymptotic covariance matrix for the parameter draws of the integrated
# intervals. mirt (or the internal estimator for unidimensional models) is
# used only in fit_scale_model(); everything downstream works from the plain
# object, so the summed-score recursion, the expected-score curves, and the
# parameter draws are generic in the number of categories and in the item
# values (sum of non-negative integers in units).
#
# MODEL TYPES.
#   grm       one factor, all items of the scale (facets, short domain scales)
#   bifactor  one general factor plus one specific factor per block
#             (components of a domain, or doublets); the EAP target is the
#             true score of the reported composite (omega-total logic)
#   composite several unidimensional component models coupled by a
#             second-order structure theta_b = lambda_b g + sqrt(1 - lambda_b^2)
#             u_b (total scores; Supplement S5.3)
# ==============================================================================

# --- Graded response model ----------------------------------------------------

#' Category probabilities of one graded item
#' @param a slope vector (length nfact), d intercepts (length K - 1, decreasing)
#' @param Theta G x nfact matrix
#' @return G x K matrix
grm_probs <- function(a, d, Theta) {
  z <- as.vector(Theta %*% a)
  Pstar <- cbind(1, sapply(d, \(dk) plogis(z + dk)), 0)
  if (!is.matrix(Pstar)) Pstar <- matrix(Pstar, nrow = 1)
  P <- Pstar[, -ncol(Pstar), drop = FALSE] - Pstar[, -1, drop = FALSE]
  pmax(P, 1e-300)
}

#' Category probabilities of item i of a model (G x K_i)
item_probs <- function(model, i, Theta) {
  if (ncol(Theta) < ncol(model$a)) {
    Th <- matrix(0, nrow(Theta), ncol(model$a)); Th[, seq_len(ncol(Theta))] <- Theta
    Theta <- Th
  }
  grm_probs(model$a[i, ], model$d[[i]], Theta)
}

#' Build a model object from its parts
new_scale_model <- function(type, items, K, values, a, d, spec = NULL,
                            par_table = NULL, vcov = NULL, constrain = list(),
                            engine = "manual", info = list()) {
  a <- as.matrix(a)
  stopifnot(nrow(a) == length(items), length(d) == length(items),
            all(lengths(d) == K - 1), length(values) == length(items))
  structure(list(type = type, items = as.integer(items), n_items = length(items),
                 K = as.integer(K), values = values, a = a, d = d,
                 spec = spec, blocks = if (!is.null(spec)) split(seq_along(spec), spec),
                 par_table = par_table, vcov = vcov, constrain = constrain,
                 engine = engine, info = info),
            class = "normpipe_model")
}

print.normpipe_model <- function(x, ...) {
  cat("normpipe item model (", x$type, ", ", x$n_items, " items, ",
      ncol(x$a), " factor(s), engine ", x$engine, ")\n", sep = "")
  invisible(x)
}

#' Parameter table of a model (one row per parameter, mirt layout)
#' columns: parnum, item (index), name (a1.., d1..), value, est
model_par_table <- function(a, d, est = NULL) {
  rows <- list(); pn <- 0L
  for (i in seq_len(nrow(a))) {
    for (f in seq_len(ncol(a))) {
      pn <- pn + 1L
      rows[[pn]] <- tibble(parnum = pn, item = i, name = paste0("a", f),
                           value = a[i, f], est = a[i, f] != 0)
    }
    for (k in seq_along(d[[i]])) {
      pn <- pn + 1L
      rows[[pn]] <- tibble(parnum = pn, item = i, name = paste0("d", k),
                           value = d[[i]][k], est = TRUE)
    }
  }
  tab <- bind_rows(rows)
  if (!is.null(est)) tab$est <- est
  tab
}

#' Rebuild a and d from a parameter table
apply_par_table <- function(model, tab) {
  a <- model$a; d <- model$d
  for (i in seq_len(model$n_items)) {
    ri <- tab[tab$item == i, ]
    ar <- ri[startsWith(ri$name, "a"), ]
    a[i, as.integer(sub("a", "", ar$name))] <- ar$value
    dr <- ri[startsWith(ri$name, "d"), ]
    d[[i]] <- dr$value[order(as.integer(sub("d", "", dr$name)))]
  }
  model$a <- a; model$d <- d; model$par_table <- tab
  model
}

# --- Calibration --------------------------------------------------------------

#' Calibrate the item model of a scale
#'
#' @param idx data frame of category indices (0-based) for the scale's items,
#'   columns in model order (item_indices on the calibration sample)
#' @param s scale specification
#' @param inst instrument
#' @param engine "mirt" (default when installed; required for bifactor),
#'   "internal" (unidimensional GRM by penalized marginal ML, no dependency)
#' @param subset for short forms scored from a longer calibration: the item
#'   positions actually calibrated (NULL = all items of s)
fit_scale_model <- function(idx, s, inst, engine = inst$measurement$engine %||% "auto",
                            ncycles = 2000, constrain_blocks = NULL) {
  stopifnot(s$model %in% c("grm", "bifactor"))
  items <- s$items
  tab <- inst$items[match(items, inst$items$item), ]
  dat <- as.data.frame(idx[, tab$col])
  # every configured category must occur, otherwise the category indices of
  # the fitted model would not match the configured values
  for (j in seq_along(items)) {
    seen <- sort(unique(dat[[j]][!is.na(dat[[j]])]))
    if (!identical(as.integer(seen), seq_len(tab$K[j]) - 1L))
      stop("Item ", tab$col[j], ": categories ", paste(setdiff(seq_len(tab$K[j]) - 1L, seen),
           collapse = ","), " do not occur in the calibration data. Collapse ",
           "them in the item configuration (categories/values) before fitting.")
  }
  if (engine == "auto") engine <- if (requireNamespace("mirt", quietly = TRUE)) "mirt" else "internal"
  spec <- if (s$model == "bifactor") rep(seq_along(s$components), lengths(s$components))
  # the item models run on VALUE-ORDERED categories: model category c of item
  # j is data index ord[[j]][c] - 1, and the model's values are increasing, so
  # reverse-keyed items keep a positive slope (as after reversing the item in
  # the original PID-5 pipeline); nothing downstream needs the data indices
  ord <- map(tab$values, order)
  for (j in seq_along(items)) dat[[j]] <- match(dat[[j]] + 1L, ord[[j]]) - 1L
  values <- map2(tab$values, ord, \(v, o) v[o])
  if (engine == "internal") {
    m <- fit_grm_internal(dat, tab$K, values, spec = spec,
                          constrain_blocks = constrain_blocks %||% isTRUE(s$extra$equal_block_slopes))
  } else {
    m <- fit_with_mirt(dat, tab$K, values, spec, ncycles,
                       constrain_blocks = constrain_blocks %||% isTRUE(s$extra$equal_block_slopes))
  }
  m$recode <- ord
  m$items <- as.integer(items)
  m$scale <- paste0(s$version, "_", s$id)
  m
}

fit_with_mirt <- function(dat, K, values, spec, ncycles, constrain_blocks = FALSE) {
  if (!requireNamespace("mirt", quietly = TRUE)) stop("mirt is not installed.")
  n <- ncol(dat)
  if (is.null(spec)) {
    mod <- mirt::mirt(dat, 1, itemtype = "graded", SE = TRUE, verbose = FALSE,
                      technical = list(NCYCLES = ncycles))
    nf <- 1L; cons <- list()
  } else {
    cons <- list()
    if (constrain_blocks) {
      # equal specific slopes within a block (doublets: just identified)
      sv <- mirt::bfactor(dat, spec, pars = "values")
      cons <- map(unique(spec), \(b) {
        pn <- sv$parnum[sv$name == paste0("a", 1 + b) & sv$est]
        pn
      })
    }
    mod <- mirt::bfactor(dat, spec, SE = TRUE, verbose = FALSE,
                         constrain = if (length(cons)) cons else NULL,
                         technical = list(NCYCLES = ncycles))
    nf <- 1L + length(unique(spec))
  }
  co <- mirt::coef(mod, simplify = TRUE)$items
  a <- as.matrix(co[, paste0("a", seq_len(nf)), drop = FALSE])
  d <- map(seq_len(n), \(i) as.numeric(co[i, paste0("d", seq_len(K[i] - 1))]))
  sv <- mirt::mod2values(mod)
  sv <- sv[sv$class == "graded", c("parnum", "item", "name", "value", "est")]
  sv$item <- match(sv$item, colnames(dat))
  tab <- as_tibble(sv)
  V <- tryCatch(mirt::extract.mirt(mod, "vcov"), error = \(e) NULL)
  parvec <- tryCatch(mirt::extract.mirt(mod, "parvec"), error = \(e) NULL)
  m <- new_scale_model(if (is.null(spec)) "grm" else "bifactor", seq_len(n), K, values,
                       a, d, spec = spec, par_table = tab, vcov = V,
                       constrain = cons, engine = "mirt",
                       info = list(converged = mirt::extract.mirt(mod, "converged"),
                                   parvec = parvec,
                                   logLik = mirt::extract.mirt(mod, "logLik")))
  m
}

#' Internal estimator: graded response model by penalized marginal maximum
#' likelihood, without mirt
#'
#' Unidimensional models and bifactor models (general factor plus one
#' specific factor per block, integrated by Gauss-Hermite quadrature over the
#' general trait and, per block, over the specific deviation, the same
#' dimension reduction as mirt::bfactor). Weak normal penalties on the log
#' slopes and the intercepts keep sparse categories finite. The covariance
#' matrix of the free parameters is the inverse observed information
#' (numerical Hessian) for small models and the outer product of the
#' per-person scores (OPG) for larger ones, transformed to the (a, d) scale
#' by the delta method. Slow for long scales: intended for tests and for
#' short instruments, mirt is the production engine.
fit_grm_internal <- function(dat, K, values, spec = NULL, gh_points = NULL,
                             gh_points_u = 9, penalty = c(a = 1.5, d = 5),
                             constrain_blocks = FALSE, maxit = 300) {
  n <- ncol(dat); N <- nrow(dat)
  bif <- !is.null(spec)
  gh_points <- gh_points %||% if (bif) 21 else 41
  gh <- statmod::gauss.quad.prob(gh_points, "normal")
  ghu <- statmod::gauss.quad.prob(gh_points_u, "normal")
  blocks <- if (bif) split(seq_len(n), spec) else list(seq_len(n))
  # parameter vector: per item log(a_g), [log(a_s)], then e_1, ..., e_{K-1}
  # with d_1 = e_1 and d_k = d_{k-1} - exp(e_k) (decreasing intercepts)
  p <- 0L; idx_a <- integer(n); idx_s <- integer(n); idx_d <- vector("list", n)
  for (i in seq_len(n)) {
    p <- p + 1L; idx_a[i] <- p
    if (bif) { p <- p + 1L; idx_s[i] <- p }
    idx_d[[i]] <- p + seq_len(K[i] - 1); p <- p + K[i] - 1L
  }
  if (bif && constrain_blocks)     # equal specific slopes within a block
    for (b in blocks) idx_s[b] <- idx_s[b[1]]
  unpack <- function(par) {
    a <- exp(par[idx_a])
    as <- if (bif) exp(par[idx_s]) else NULL
    d <- map(seq_len(n), \(i) { e <- par[idx_d[[i]]]; cumsum(c(e[1], -exp(e[-1]))) })
    list(a = a, as = as, d = d)
  }
  X <- as.matrix(dat) + 1L   # 1-based categories, NA allowed
  # per-person log-likelihood contributions
  ll_person <- function(par) {
    pr <- unpack(par)
    L <- matrix(1, N, gh_points)
    for (b in blocks) {
      if (!bif) {
        for (i in b) {
          P <- grm_probs(pr$a[i], pr$d[[i]], matrix(gh$nodes))
          ok <- !is.na(X[, i])
          L[ok, ] <- L[ok, ] * t(P)[X[ok, i], , drop = FALSE]
        }
      } else {
        Lb <- matrix(0, N, gh_points)
        for (q in seq_len(gh_points_u)) {
          Lq <- matrix(1, N, gh_points)
          Theta <- cbind(gh$nodes, ghu$nodes[q])
          for (i in b) {
            P <- grm_probs(c(pr$a[i], pr$as[i]), pr$d[[i]], Theta)
            ok <- !is.na(X[, i])
            Lq[ok, ] <- Lq[ok, ] * t(P)[X[ok, i], , drop = FALSE]
          }
          Lb <- Lb + ghu$weights[q] * Lq
        }
        L <- L * Lb
      }
    }
    log(as.vector(L %*% gh$weights))
  }
  penal <- function(par) {
    pr <- unpack(par)
    sum(par[unique(c(idx_a, idx_s[idx_s > 0]))]^2) / (2 * penalty["a"]^2) +
      sum(unlist(map(pr$d, \(dd) dd^2))) / (2 * penalty["d"]^2)
  }
  objective <- function(par) -(sum(ll_person(par)) - penal(par))
  # start values: log a = 0 (specific: log 0.5), thresholds from the marginal
  # category proportions
  start <- numeric(p)
  if (bif) start[idx_s] <- log(0.5)
  for (i in seq_len(n)) {
    pr <- prop.table(table(factor(X[, i], levels = seq_len(K[i]))))
    cum <- rev(cumsum(rev(pr)))[-1]
    d0 <- stats::qlogis(pmin(pmax(cum, 0.01), 0.99)) * 1.7
    start[idx_d[[i]]] <- c(d0[1], log(pmax(-diff(d0), 0.05)))
  }
  opt <- optim(start, objective, method = "BFGS",
               control = list(maxit = maxit, reltol = 1e-9))
  pr <- unpack(opt$par)
  a <- if (bif) {
    A <- matrix(0, n, 1 + length(blocks)); A[, 1] <- pr$a
    for (b in seq_along(blocks)) A[blocks[[b]], 1 + b] <- pr$as[blocks[[b]]]
    A
  } else matrix(pr$a, ncol = 1)
  d <- pr$d
  tab <- model_par_table(a, d)
  # free parameters on the (a, d) scale, in the order of tab (est == TRUE)
  first_in_block <- rep(TRUE, n)
  if (bif && constrain_blocks) for (b in blocks) first_in_block[b[-1]] <- FALSE
  to_ad <- function(par) {
    u <- unpack(par)
    unlist(map(seq_len(n), \(i) c(u$a[i], if (bif && first_in_block[i]) u$as[i], u$d[[i]])))
  }
  V <- NULL
  free_p <- length(unique(c(idx_a, idx_s[idx_s > 0], unlist(idx_d))))
  Vp <- tryCatch({
    if (free_p <= 40) {
      H <- numDeriv::hessian(objective, opt$par)
      solve((H + t(H)) / 2)
    } else {
      Sc <- numDeriv::jacobian(ll_person, opt$par)     # N x p score matrix
      solve(crossprod(Sc))
    }
  }, error = \(e) NULL)
  if (!is.null(Vp)) {
    J <- numDeriv::jacobian(to_ad, opt$par)
    V <- J %*% Vp %*% t(J)
    est_rows <- tab$est
    if (bif && constrain_blocks)
      est_rows <- est_rows & !(tab$item %in% which(!first_in_block) & grepl("^a[2-9]", tab$name))
    rownames(V) <- colnames(V) <- paste0(tab$name[est_rows], ".", tab$parnum[est_rows])
  }
  cons <- list()
  if (bif && constrain_blocks)
    cons <- map(seq_along(blocks), \(b) tab$parnum[tab$item %in% blocks[[b]] & tab$name == paste0("a", 1 + b)])
  new_scale_model(if (bif) "bifactor" else "grm", seq_len(n), K, values, a, d,
                  spec = spec, par_table = tab, vcov = V, constrain = cons,
                  engine = "internal",
                  info = list(converged = opt$convergence == 0, logLik = -opt$value,
                              parvec = to_ad(opt$par)))
}

# --- Parameter draws (Supplement S5.4) -----------------------------------------

#' Draw item-parameter sets from the asymptotic distribution of a model
#'
#' Multivariate normal draws on the free parameters with the stored covariance
#' matrix; draws whose intercepts are not strictly decreasing for every item
#' are rejected and redrawn. Equality constraints (mirt) are honoured by
#' drawing the first member of a group and copying it. Without a usable
#' covariance matrix the point estimates are returned for every draw.
draw_model_params <- function(model, n_draws, max_tries = 50) {
  tab <- model$par_table; V <- model$vcov
  if (is.null(tab) || is.null(V)) return(rep(list(model), n_draws))
  est_idx <- which(tab$est)
  cons <- model$constrain %||% list()
  drop_pn <- if (length(cons)) unlist(lapply(cons, \(g) g[-1])) else integer(0)
  est_idx <- est_idx[!(tab$parnum[est_idx] %in% drop_pn)]
  mu <- if (!is.null(model$info$parvec)) model$info$parvec else tab$value[est_idx]
  if (length(mu) != nrow(V)) {
    # mirt reports the vcov for the identifiable subset only; align by parnum
    vc_pn <- suppressWarnings(as.integer(sub("^.*\\.", "", colnames(V))))
    if (anyNA(vc_pn)) { message("draw_model_params: cannot align vcov; point estimates used."); return(rep(list(model), n_draws)) }
    sel <- match(vc_pn, tab$parnum[est_idx])
    if (anyNA(sel)) { message("draw_model_params: vcov names unknown; point estimates used."); return(rep(list(model), n_draws)) }
    draw_idx <- est_idx[sel]; mu <- mu[sel]
  } else draw_idx <- est_idx
  if (length(mu) / length(est_idx) < 0.5) {
    message("draw_model_params: information matrix covers only ", length(mu),
            " of ", length(est_idx), " parameters; point estimates used.")
    return(rep(list(model), n_draws))
  }
  V <- (V + t(V)) / 2
  ridge <- 0; L <- tryCatch(chol(V), error = \(e) NULL)
  while (is.null(L)) {
    ridge <- if (ridge == 0) 1e-10 else ridge * 10
    L <- tryCatch(chol(V + diag(ridge, nrow(V))), error = \(e) NULL)
  }
  d_rows <- split(which(grepl("^d[0-9]+$", tab$name)), tab$item[grepl("^d[0-9]+$", tab$name)])
  apply_cons <- function(tb) {
    for (g in cons) { r <- match(g, tb$parnum); tb$value[r] <- tb$value[r[1]] }
    tb
  }
  draw_valid <- function() {
    for (t in seq_len(max_tries)) {
      v <- as.numeric(mu + t(L) %*% rnorm(length(mu)))
      tb <- tab; tb$value[draw_idx] <- v; tb <- apply_cons(tb)
      ok <- all(map_lgl(d_rows, \(r) {
        dd <- tb$value[r][order(as.integer(sub("d", "", tb$name[r])))]
        all(diff(dd) < 0)
      }))
      if (ok) return(tb)
    }
    stop("draw_model_params: no valid draw within ", max_tries, " tries.")
  }
  map(seq_len(n_draws), \(k) apply_par_table(model, draw_valid()))
}

# --- Blocks: the common structure of all EAP computations ---------------------
#
# A block list describes how the items of a scale are grouped and how the
# latent variables of a block depend on the general trait g and a block
# deviation u ~ N(0, 1):
#   model      item model that holds the block's items
#   pos        item positions inside that model
#   theta_fun  function(g, u) -> Theta matrix for item_probs()
#   bw         weight of the block's summed value in the reported score
#              (times unit), so that the composite true score is
#              sum_b bw_b * E[block sum_b]
# grm: one block, theta = g, no block deviation.
# bifactor: one block per specific factor, Theta = (g, 0, .., u at the block's
#   column); composite: one block per component model, theta = lambda g +
#   sqrt(1 - lambda^2) u.

block_weights <- function(s, n_per_block) {
  N <- sum(n_per_block)
  switch(s$metric,
    sum = rep(1, length(n_per_block)),
    mean = rep(1 / N, length(n_per_block)),
    mean_of_means = 1 / (n_per_block * length(n_per_block)))
}

#' Block list of a scale given its model(s)
#' @param models for grm/bifactor one model; for composite a named list of
#'   component models (names = component ids) and `lambda` the loadings
#' @param subset for grm/bifactor: positions (within the model) of the items
#'   actually scored (short forms calibrated on the long form); NULL = all
make_blocks <- function(s, models, subset = NULL, lambda = NULL) {
  if (s$model == "composite") {
    stopifnot(!is.null(lambda), setequal(names(models), names(s$components)))
    ns <- lengths(s$components)
    bw <- block_weights(s, ns)
    map(seq_along(s$components), \(b) {
      m <- models[[names(s$components)[b]]]
      pos <- match(s$components[[b]], m$items)
      stopifnot(!anyNA(pos))
      lam <- lambda[[names(s$components)[b]]]
      list(model = m, pos = pos, bw = bw[b],
           theta_fun = \(g, u) matrix(lam * g + sqrt(1 - lam^2) * u),
           has_u = TRUE, id = names(s$components)[b])
    })
  } else if (s$model == "bifactor") {
    m <- models
    pos_all <- subset %||% seq_len(m$n_items)
    bl <- map(seq_along(m$blocks), \(b) intersect(m$blocks[[b]], pos_all))
    keep <- lengths(bl) > 0
    ns <- lengths(bl[keep]); bw <- block_weights(s, ns)
    nf <- ncol(m$a)
    map(which(keep), \(b) list(
      model = m, pos = bl[[b]], bw = bw[match(b, which(keep))],
      theta_fun = \(g, u) { Th <- matrix(0, length(g), nf); Th[, 1] <- g; Th[, 1 + b] <- u; Th },
      has_u = TRUE, id = paste0("block", b)))
  } else {
    m <- models
    pos <- subset %||% seq_len(m$n_items)
    list(list(model = m, pos = pos, bw = block_weights(s, length(pos)),
              theta_fun = \(g, u) matrix(g), has_u = FALSE, id = "all"))
  }
}

#' Lord-Wingersky recursion over integer item values
#' @param P list of G x K_i probability matrices, vals list of value vectors
#' @return G x (Smax + 1) likelihood matrix L(S | theta), S = 0..Smax
lw_recursion <- function(P, vals) {
  Smax <- sum(map_int(vals, max))
  G <- nrow(P[[1]])
  L <- matrix(0, G, Smax + 1); L[, 1] <- 1
  cur <- 0L
  for (i in seq_along(P)) {
    v <- vals[[i]]; Pi <- P[[i]]
    Ln <- matrix(0, G, Smax + 1)
    for (k in seq_along(v)) {
      cols <- seq_len(cur + 1) + v[k]
      Ln[, cols] <- Ln[, cols] + L[, seq_len(cur + 1)] * Pi[, k]
    }
    L <- Ln; cur <- cur + max(v)
  }
  L
}

#' Prepare one block: per node q the likelihood L_q(S_b | g) and the expected
#' block sum tau_b(g, q), on the g grid
block_prep <- function(bl, grid, gh) {
  G <- length(grid); Q <- length(gh$nodes)
  m <- bl$model
  vals <- m$values[bl$pos]
  Smax <- sum(map_int(vals, max))
  Lq <- array(0, c(G, Smax + 1, Q)); tau <- matrix(0, G, Q)
  for (q in seq_len(Q)) {
    Theta <- bl$theta_fun(grid, if (bl$has_u) gh$nodes[q] else 0)
    P <- map(bl$pos, \(i) item_probs(m, i, Theta))
    tau[, q] <- Reduce(`+`, map2(P, vals, \(Pi, v) as.vector(Pi %*% v)))
    Lq[, , q] <- lw_recursion(P, vals)
    if (!bl$has_u) { if (Q > 1) for (q2 in 2:Q) { Lq[, , q2] <- Lq[, , 1]; tau[, q2] <- tau[, 1] }; break }
  }
  Lmarg <- matrix(0, G, Smax + 1)
  for (q in seq_len(Q)) Lmarg <- Lmarg + gh$weights[q] * Lq[, , q]
  list(Lq = Lq, Lmarg = Lmarg, tau = tau, Smax = Smax, bw = bl$bw)
}

#' Expected reported score as a function of the general trait g
#'
#' @return list(theta = grid, mean_score = expected raw score), the format
#'   consumed by score_dist_to_theta_prior() and make_Tmap()
expected_score_curve <- function(blocks, s, grid = seq(-6, 6, length.out = 241),
                                 gh_points = 21, prep = NULL) {
  gh <- statmod::gauss.quad.prob(gh_points, "normal")
  B <- prep %||% map(blocks, block_prep, grid = grid, gh = gh)
  ms <- Reduce(`+`, map(B, \(b) b$bw * as.vector(b$tau %*% gh$weights))) * s$unit
  list(theta = grid, mean_score = ms)
}

#' Summed-score EAP table (Supplement S5.2)
#'
#' Exact for one block without block deviation (unidimensional GRM: posterior
#' of g per summed score); Monte Carlo over the block deviations otherwise
#' (bifactor and composite: the composite true score is sampled exactly given
#' the total summed score, as in the original eapsum_bifactor_table).
#'
#' @param blocks block list (make_blocks)
#' @param s scale
#' @param prior_theta,prior_w population prior on the g grid
#' @param sum_scores summed values (integer units) at which to tabulate;
#'   default all attainable sums of the blocks' items
#' @return tibble sum_score, raw, theta_eap, theta_lo, theta_hi, tau_eap,
#'   tau_sd, p_marg; with return_posterior = TRUE a list with the posterior
#'   over the g grid per summed score (for the integrated intervals)
eapsum_table <- function(blocks, s, prior_theta, prior_w, sum_scores = NULL,
                         gh_points = 21, probs = c(.025, .975), n_mc = 1500,
                         return_posterior = FALSE, prep = NULL) {
  gh <- statmod::gauss.quad.prob(gh_points, "normal")
  G <- length(prior_theta); Q <- gh_points
  B <- prep %||% map(blocks, block_prep, grid = prior_theta, gh = gh)
  nb <- length(B)
  unit <- s$unit
  mean_score <- Reduce(`+`, map(B, \(b) b$bw * as.vector(b$tau %*% gh$weights))) * unit
  convo <- function(A, Bm) {
    Ln <- matrix(0, G, ncol(A) + ncol(Bm) - 1)
    for (k in seq_len(ncol(Bm))) Ln[, k:(k + ncol(A) - 1)] <- Ln[, k:(k + ncol(A) - 1)] + A * Bm[, k]
    Ln
  }
  Csuf <- vector("list", nb); Csuf[[nb]] <- B[[nb]]$Lmarg
  if (nb > 1) for (b in (nb - 1):1) Csuf[[b]] <- convo(B[[b]]$Lmarg, Csuf[[b + 1]])
  Ltot <- Csuf[[1]]
  Smax_tot <- ncol(Ltot) - 1
  all_sums <- 0:Smax_tot
  joint <- sweep(Ltot, 1, prior_w, `*`)
  p_marg <- colSums(joint)
  attainable <- p_marg > 1e-300
  if (is.null(sum_scores)) sum_scores <- all_sums[attainable]
  stopifnot(all(sum_scores %in% all_sums))
  raw <- sum_to_raw(sum_scores, s)
  exact <- nb == 1 && !blocks[[1]]$has_u
  pseudo_theta <- function(tau) approx(mean_score, prior_theta, xout = tau, rule = 2)$y
  if (exact) {
    post <- sweep(joint, 2, pmax(p_marg, 1e-300), `/`)
    tau_theta <- B[[1]]$bw * B[[1]]$tau[, 1] * unit
    if (return_posterior)
      return(list(sum_scores = sum_scores, raw = raw, theta = prior_theta,
                  post = post[, sum_scores + 1, drop = FALSE]))
    return(map_dfr(sum_scores, \(sc) {
      p <- post[, sc + 1]; cdf <- cumsum(p)
      tau_eap <- sum(tau_theta * p)
      tibble(sum_score = sc, raw = sum_to_raw(sc, s),
             theta_eap = sum(prior_theta * p),
             theta_sd = sqrt(sum((prior_theta - sum(prior_theta * p))^2 * p)),
             theta_lo = prior_theta[which.max(cdf >= probs[1])],
             theta_hi = prior_theta[which.max(cdf >= probs[2])],
             tau_eap = tau_eap, tau_sd = sqrt(sum((tau_theta - tau_eap)^2 * p)),
             p_marg = p_marg[sc + 1])
    }))
  }
  post_g <- sweep(joint, 2, pmax(p_marg, 1e-300), `/`)
  cum_mat <- function(W) W %*% upper.tri(diag(ncol(W)), diag = TRUE)
  draw_rows <- function(W) {
    rs <- rowSums(W); W[rs <= 0, 1] <- 1; rs <- rowSums(W)
    rowSums(cum_mat(W) < runif(nrow(W)) * rs)
  }
  sample_tau <- function(sc) {
    gs <- sample.int(G, n_mc, replace = TRUE, prob = post_g[, sc + 1])
    rem <- rep(sc, n_mc); tau_r <- numeric(n_mc)
    for (b in seq_len(nb)) {
      Bb <- B[[b]]
      if (b < nb) {
        k <- 0:Bb$Smax; K <- length(k)
        P <- Bb$Lmarg[gs, , drop = FALSE]
        nxt <- outer(rem, k, `-`)
        Cn <- Csuf[[b + 1]]
        ok <- nxt >= 0 & nxt <= ncol(Cn) - 1
        Cv <- matrix(0, n_mc, K)
        Cv[ok] <- Cn[cbind(rep(gs, times = K)[ok], nxt[ok] + 1)]
        sb <- draw_rows(P * Cv)
      } else sb <- rem
      Wq <- matrix(Bb$Lq[cbind(rep(gs, Q), rep(sb + 1, Q), rep(seq_len(Q), each = n_mc))],
                   n_mc, Q) * rep(gh$weights, each = n_mc)
      qb <- draw_rows(Wq) + 1
      tau_r <- tau_r + Bb$bw * Bb$tau[cbind(gs, qb)] * unit
      rem <- rem - sb
    }
    tau_r
  }
  if (return_posterior) {
    post <- matrix(0, G, length(sum_scores))
    brk <- c(-Inf, (prior_theta[-1] + prior_theta[-G]) / 2, Inf)
    for (j in seq_along(sum_scores)) {
      pt <- pseudo_theta(sample_tau(sum_scores[j]))
      post[, j] <- tabulate(findInterval(pt, brk), nbins = G) / n_mc
    }
    return(list(sum_scores = sum_scores, raw = raw, theta = prior_theta, post = post))
  }
  map_dfr(sum_scores, \(sc) {
    tau_r <- sample_tau(sc)
    qs <- quantile(tau_r, probs)
    tibble(sum_score = sc, raw = sum_to_raw(sc, s),
           theta_eap = pseudo_theta(mean(tau_r)),
           theta_sd = NA_real_,
           theta_lo = pseudo_theta(qs[1]), theta_hi = pseudo_theta(qs[2]),
           tau_eap = mean(tau_r), tau_sd = sd(tau_r), p_marg = p_marg[sc + 1])
  })
}

#' Summed scores at which a long scale is tabulated: all attainable sums when
#' there are at most max_rows, otherwise about max_rows evenly spaced ones
eap_sum_rows <- function(s, max_rows = 61) {
  if (length(s$sums) <= max_rows) return(s$sums)
  s$sums[unique(round(seq(1, length(s$sums), length.out = max_rows)))]
}

# --- Composite route: second-order loadings (Supplement S5.3) -----------------

#' Model-implied reliability of a block's summed value in the calibration
#' population (theta ~ N(0, 1))
block_sum_reliability <- function(model, pos = NULL, gh_points = 41) {
  gh <- statmod::gauss.quad.prob(gh_points, "normal")
  pos <- pos %||% seq_len(model$n_items)
  Theta <- matrix(gh$nodes)
  m <- numeric(gh_points); v <- numeric(gh_points)
  for (i in pos) {
    P <- item_probs(model, i, Theta); k <- model$values[[i]]
    e1 <- as.vector(P %*% k); e2 <- as.vector(P %*% k^2)
    m <- m + e1; v <- v + (e2 - e1^2)
  }
  mu <- sum(gh$weights * m)
  var_true <- sum(gh$weights * (m - mu)^2)
  var_true / (var_true + sum(gh$weights * v))
}

#' Second-order loadings of the component traits of a composite scale,
#' estimated from the component summed values in the calibration data
#' (disattenuated correlations, principal-axis iteration)
total_loadings <- function(idx, s, inst, models, n_iter = 50) {
  sums <- map_dfc(set_names(names(s$components)), \(cp) {
    v <- item_values(idx, inst, s$components[[cp]])
    rowSums(as.matrix(v))
  })
  rel <- map_dbl(names(s$components), \(cp) {
    m <- models[[cp]]; block_sum_reliability(m, match(s$components[[cp]], m$items))
  })
  R <- cor(sums, use = "pairwise.complete.obs")
  R <- R / sqrt(outer(rel, rel)); diag(R) <- 1
  h2 <- tryCatch(1 - 1 / diag(solve(R)), error = \(e) rep(0.5, ncol(R)))
  h2 <- pmin(pmax(h2, 0.05), 0.98)
  for (it in seq_len(n_iter)) {
    Rr <- R; diag(Rr) <- h2
    e <- eigen(Rr, symmetric = TRUE)
    lam <- e$vectors[, 1] * sqrt(max(e$values[1], 1e-8))
    if (sum(lam) < 0) lam <- -lam
    h2_new <- pmin(lam^2, 0.98)
    if (max(abs(h2_new - h2)) < 1e-6) { h2 <- h2_new; break }
    h2 <- h2_new
  }
  lam <- pmin(pmax(lam, 0.10), 0.95)
  set_names(lam, names(s$components))
}
