# ==============================================================================
# cache.R -- stage caches that know what they were computed from
#
# Every cached result of run_norming() is stored as <name>.rds in out_dir
# together with a key in cache_manifest.csv. The key is a hash of everything
# the result depends on (configuration slice, input data, keys of the results
# it was built from). A result is reused only when its key is unchanged;
# otherwise it is recomputed, and the outdated file is moved to
# out_dir/superseded/ (nothing is deleted).
# ==============================================================================

#' Stage caches with dependency keys
#'
#' [run_norming()] caches every intermediate result as `<name>.rds` in
#' `out_dir` and records in `cache_manifest.csv` a key for it: a hash of the
#' part of the configuration, the data and the upstream results that the
#' result depends on. On the next run a result is reused when its key is
#' unchanged and recomputed otherwise, so that a change of the configuration,
#' the data or the package's computation only recomputes what depends on it.
#' An outdated file is not deleted but moved to `out_dir/superseded/`.
#'
#' `cache_key()` builds a key from any R objects. `cache_get()` returns the
#' cached result for a key or computes and stores it. `cache_status()` lists
#' the manifest of an output folder.
#'
#' @param ... Objects the result depends on (for `cache_key()`).
#' @param out_dir Output folder.
#' @param name Name of the result (file `<name>.rds`).
#' @param key Key of the result ([cache_key()]).
#' @param compute Function without arguments that computes the result.
#' @param say Function for progress messages (default: none).
#' @return `cache_key()`: a string. `cache_get()`: the result.
#'   `cache_status()`: a tibble with `name`, `key`, `version`, `time` and
#'   `exists`.
#' @examples
#' out <- file.path(tempdir(), "np-cache")
#' k <- cache_key(1:3, "a")
#' cache_get(out, "demo", k, function() sum(1:3))
#' cache_get(out, "demo", k, function() stop("not computed again"))
#' cache_status(out)
#' @export
cache_key <- function(...) rlang::hash(list(...))

#' @rdname cache_key
#' @export
cache_get <- function(out_dir, name, key, compute, say = NULL) {
  f <- file.path(out_dir, paste0(name, ".rds"))
  man <- read_manifest(out_dir)
  hit <- match(name, man$name)
  if (file.exists(f) && !is.na(hit) && identical(man$key[hit], key)) return(readRDS(f))
  if (file.exists(f)) {
    old <- file.path(out_dir, "superseded")
    dir.create(old, showWarnings = FALSE, recursive = TRUE)
    to <- file.path(old, paste0(name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
    file.rename(f, to)
    if (!is.null(say)) say(name, ": inputs changed, recomputing (previous result moved to superseded/)")
  }
  value <- compute()
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(f, ".tmp")
  saveRDS(value, tmp); file.rename(tmp, f)
  cache_record(out_dir, name, key)
  value
}

#' Record the key of a result that was written outside cache_get()
#' @noRd
cache_record <- function(out_dir, name, key) {
  man <- read_manifest(out_dir)
  man <- man[man$name != name, ]
  man <- rbind(man, data.frame(name = name, key = key,
                               version = as.character(utils::packageVersion("normpipe")),
                               time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                               stringsAsFactors = FALSE))
  f <- file.path(out_dir, "cache_manifest.csv")
  tmp <- paste0(f, ".tmp")
  utils::write.csv(man, tmp, row.names = FALSE); file.rename(tmp, f)
  invisible(key)
}

#' Is the cached result for this key present?
#' @noRd
cache_valid <- function(out_dir, name, key) {
  man <- read_manifest(out_dir)
  hit <- match(name, man$name)
  file.exists(file.path(out_dir, paste0(name, ".rds"))) && !is.na(hit) && identical(man$key[hit], key)
}

read_manifest <- function(out_dir) {
  f <- file.path(out_dir, "cache_manifest.csv")
  if (!file.exists(f))
    return(data.frame(name = character(0), key = character(0), version = character(0),
                      time = character(0), stringsAsFactors = FALSE))
  utils::read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
}

#' @rdname cache_key
#' @export
cache_status <- function(out_dir) {
  man <- read_manifest(out_dir)
  as_tibble(man) %>% mutate(exists = file.exists(file.path(out_dir, paste0(.data$name, ".rds"))))
}
