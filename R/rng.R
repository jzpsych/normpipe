# ==============================================================================
# rng.R -- random numbers with a fixed seed per step, without side effects
#
# Every random step of the pipeline (bootstrap of a draw frame, imputation,
# projection of a transfer model, transfer fit, item-parameter draws,
# Monte Carlo of the summed-score posterior) runs with its own seed, derived
# from the name of the step and the scale, and restores the random-number
# state of the caller afterwards. A result therefore does not depend on what
# ran before it, on which results came from the cache, or on the state of
# the session.
#
# Steps that have had a seed since earlier versions keep its formula
# (sum of the character codes of a tag), so that their results stay the
# same. `uncertainty: seed` in the instrument replaces all seeds by seeds
# derived from that number.
# ==============================================================================

#' Evaluate `code` with a seed, then restore the caller's random-number state
#' @noRd
with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = genv) else
    if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv))
  set.seed(seed)
  code
}

#' A seed from strings: polynomial hash modulo 2^31 - 1, the same on every
#' platform and R version, reduced below 10^9 so that seed + chain number
#' stays an integer
#' @noRd
seed_of <- function(...) {
  x <- utf8ToInt(paste(c(...), collapse = "\r"))
  h <- 0
  for (ch in x) h <- (h * 131 + ch) %% 2147483647
  as.integer(h %% 999999937)
}

#' Seed of one random step of the pipeline
#'
#' `legacy = TRUE` for steps that have had a fixed seed since earlier
#' versions: without `uncertainty$seed` their formula is kept. Steps without
#' an earlier seed get one from seed_of().
#' @noRd
step_seed <- function(inst, tag, legacy = FALSE) {
  base <- inst$uncertainty$seed
  if (!is.null(base)) return(seed_of(as.character(base), tag))
  if (legacy) sum(utf8ToInt(tag)) else seed_of(tag)
}

#' Set a seed for the rest of the calling function; the caller's
#' random-number state is restored when that function exits
#' @noRd
local_seed <- function(seed, envir = parent.frame()) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  restore <- function() {
    if (had) assign(".Random.seed", old, envir = genv) else
      if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
  }
  do.call(on.exit, list(bquote(.(restore)()), add = TRUE), envir = envir)
  set.seed(seed)
  invisible(seed)
}
