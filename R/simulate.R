# ==============================================================================
# simulate.R -- synthetic data for tests and stress checks
# ==============================================================================

#' Simulate a reference sample for an instrument from graded response models
#'
#' One trait per version-"full" facet-level scale is drawn from a correlated
#' normal (shared general factor), items follow a GRM with random slopes and
#' thresholds on the configured categories; age and gender shift the traits.
#' Returns list(meta, items, theta, params).
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
  age <- round(runif(N, 18, 80)); gender <- sample(c("m", "f"), N, TRUE)
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
  list(meta = tibble(id = paste0("sim_", seq_len(N)), age = age, gender = gender,
                     weight = exp(rnorm(N, 0, 0.3))),
       items = as.data.frame(items), theta = theta, trait = trait,
       params = list(a = a, d = d))
}
