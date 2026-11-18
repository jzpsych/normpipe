# ==============================================================================
# simulate.R -- synthetic data for tests and stress checks
# ==============================================================================

#' Simulate a reference sample for an instrument from graded response models
#'
#' One trait per version-"full" facet-level scale is drawn from a correlated
#' normal (shared general factor), items follow a GRM with random slopes and
#' thresholds on the configured categories. The sample has the variables
#' `age` and `gender`, which shift the traits (their range and codes follow
#' the norm groups when the instrument uses these variables), and uniformly
#' drawn values of every other variable of the norm groups. Intended for
#' tests, examples and simulation checks.
#'
#' @param inst Instrument ([read_instrument()]).
#' @param N Number of respondents.
#' @param seed Random seed.
#' @param general Loading of every trait on the shared general factor.
#' @param age_effect Linear age effect on the traits (per year, centred at 50).
#' @param gender_effect Shift of every second trait for `gender == "f"`.
#' @return List with `meta` (`id`, `age`, `gender`, `weight` and any other
#'   variable of the norm groups), `items` (data
#'   frame of response codes, columns named as in the item table), `theta`
#'   (true traits), `trait` (trait index per item) and `params` (true slopes
#'   and intercepts).
#' @examples
#' inst <- read_instrument(system.file("instruments/demo7/demo7.yml",
#'                                     package = "normpipe"))
#' sim <- simulate_reference(inst, N = 200)
#' str(sim$meta)
#' @export
simulate_reference <- function(inst, N = 2000, seed = 1, general = 0.6,
                               age_effect = 0.01, gender_effect = 0.3) {
  set.seed(seed)
  it <- inst$items
  n_items <- nrow(it)
  # assign every item to one trait: the first facet-level scale of the first
  # version that contains it (fallback: one trait for the rest)
  trait <- rep(NA_integer_, n_items); trait_names <- character(0)
  for (v in inst$versions) for (s in v$scales) if (is.null(s$components)) {
    pos <- match(s$items, it$item)
    if (any(is.na(trait[pos]))) {
      trait_names <- c(trait_names, paste0(v$id, "_", s$id))
      trait[pos][is.na(trait[pos])] <- length(trait_names)
    }
  }
  if (anyNA(trait)) { trait_names <- c(trait_names, "rest"); trait[is.na(trait)] <- length(trait_names) }
  Tn <- length(trait_names)
  age <- round(runif(N, sim_range(inst, "age")[1], sim_range(inst, "age")[2]))
  gender <- sample(inst$groups$vars$gender$levels %||% c("m", "f"), N, TRUE)
  g <- rnorm(N)
  theta <- sapply(seq_len(Tn), \(t) general * g + sqrt(1 - general^2) * rnorm(N) +
                    age_effect * (age - 50) + gender_effect * (gender == "f") * (t %% 2))
  theta <- matrix(theta, N, Tn)
  a <- runif(n_items, 0.8, 2.2)
  d <- map(it$K, \(K) sort(seq(2, -2, length.out = K - 1) + rnorm(K - 1, 0, 0.3),
                           decreasing = TRUE))
  items <- map(seq_len(n_items), \(i) {
    P <- grm_probs(a[i], d[[i]], matrix(theta[, trait[i]]))
    # inverse-CDF draw: category = number of cumulative probabilities below u
    u <- runif(N)
    k <- rowSums(t(apply(P, 1, cumsum)) < u)
    # reverse-keyed items: the data code falls with the trait
    if (isTRUE(it$reverse[i])) k <- (it$K[i] - 1L) - k
    it$categories[[i]][k + 1]
  })
  names(items) <- it$col
  items <- as_tibble(items)
  meta <- tibble(id = paste0("sim_", seq_len(N)), age = age, gender = gender,
                 weight = exp(rnorm(N, 0, 0.3)))
  for (v in setdiff(group_sources(inst), names(meta)))
    meta[[v]] <- if (v %in% categorical_sources(inst)) sample(inst$groups$vars[[v]]$levels, N, TRUE) else
      round(runif(N, sim_range(inst, v)[1], sim_range(inst, v)[2]))
  list(meta = meta[c("id", setdiff(names(meta), c("id", "weight")), "weight")],
       items = as.data.frame(items), theta = theta, trait = trait,
       params = list(a = a, d = d))
}

#' Range of a numeric variable of the norm groups for simulated data
#' (default 18 to 80)
#' @noRd
sim_range <- function(inst, v) {
  g <- inst$groups
  if (identical(g$continuous$variable, v)) return(range(g$continuous$values))
  br <- g$vars[[v]]$breaks
  br <- br[is.finite(br)]
  if (length(br) >= 2) return(range(br))
  c(18, 80)
}
