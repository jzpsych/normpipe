# ==============================================================================
# cache.R -- stage caches that know what they were computed from
#
# Every cached result of run_norming() is stored as <name>.rds in
# out_dir/cache together with a key in cache/cache_manifest.csv. The key is a
# hash of everything the result depends on (configuration slice, input data,
# keys of the results it was built from). A result is reused only when its
# key is unchanged; otherwise it is recomputed, and the outdated file is
# moved to cache/superseded/ (nothing is deleted). The tables for the user
# are written to out_dir itself (export.R).
# ==============================================================================

#' Stage caches with dependency keys
#'
#' [run_norming()] caches every intermediate result as `<name>.rds` in
#' `out_dir/cache` and records in `cache/cache_manifest.csv` a key for it: a
#' hash of the
#' part of the configuration, the data and the upstream results that the
#' result depends on. On the next run a result is reused when its key is
#' unchanged and recomputed otherwise, so that a change of the configuration,
#' the data or the package's computation only recomputes what depends on it.
#' An outdated file is not deleted but moved to `cache/superseded/`.
#'
#' `cache_key()` builds a key from any R objects. `cache_get()` returns the
#' cached result for a key or computes and stores it. `cache_status()` lists
#' the manifest of an output folder.
#'
#' @param ... Objects the result depends on (for `cache_key()`).
#' @param out_dir Output folder of a run (the cache is its subfolder
#'   `cache`).
#' @param name Name of the result (file `cache/<name>.rds`).
#' @param key Key of the result ([cache_key()]).
#' @param compute Function without arguments that computes the result.
#' @param say Function for progress messages (default: none).
#' @return `cache_key()`: a string. `cache_get()`: the result.
#'   `cache_status()`: a tibble with `name`, `key`, `version`, `time`,
#'   `adopted` (`"TRUE"` for a result taken over by [migrate_output()] that
#'   has not been used since) and `exists`.
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
  f <- cache_file(out_dir, name)
  man <- read_manifest(out_dir)
  hit <- match(name, man$name)
  if (file.exists(f) && !is.na(hit) && identical(man$key[hit], key)) return(readRDS(f))
  if (file.exists(f) && !is.na(hit) && identical(man$adopted[hit], "TRUE")) {
    # a result taken over by migrate_output(): accepted once, then keyed
    cache_record(out_dir, name, key)
    return(readRDS(f))
  }
  if (file.exists(f)) {
    old <- file.path(cache_dir(out_dir), "superseded")
    dir.create(old, showWarnings = FALSE, recursive = TRUE)
    to <- file.path(old, paste0(name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
    file.rename(f, to)
    if (!is.null(say)) say(name, ": inputs changed, recomputing (previous result moved to superseded/)")
  }
  value <- compute()
  dir.create(cache_dir(out_dir), showWarnings = FALSE, recursive = TRUE)
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
                               adopted = "", stringsAsFactors = FALSE))
  write_manifest(out_dir, man)
  invisible(key)
}

#' Is the cached result for this key present?
#' @noRd
cache_valid <- function(out_dir, name, key) {
  man <- read_manifest(out_dir)
  hit <- match(name, man$name)
  file.exists(cache_file(out_dir, name)) && !is.na(hit) &&
    (identical(man$key[hit], key) || identical(man$adopted[hit], "TRUE"))
}

cache_dir <- function(out_dir) file.path(out_dir, "cache")
cache_file <- function(out_dir, name) file.path(cache_dir(out_dir), paste0(name, ".rds"))

write_manifest <- function(out_dir, man) {
  dir.create(cache_dir(out_dir), showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache_dir(out_dir), "cache_manifest.csv")
  tmp <- paste0(f, ".tmp")
  utils::write.csv(man, tmp, row.names = FALSE); file.rename(tmp, f)
  invisible(f)
}

read_manifest <- function(out_dir) {
  f <- file.path(cache_dir(out_dir), "cache_manifest.csv")
  if (!file.exists(f))
    return(data.frame(name = character(0), key = character(0), version = character(0),
                      time = character(0), adopted = character(0), stringsAsFactors = FALSE))
  man <- utils::read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
  if (is.null(man$adopted)) man$adopted <- rep("", nrow(man))
  man$adopted[is.na(man$adopted)] <- ""
  man
}

#' @rdname cache_key
#' @export
cache_status <- function(out_dir) {
  man <- read_manifest(out_dir)
  as_tibble(man) %>% mutate(exists = file.exists(cache_file(out_dir, .data$name)))
}

# --- Output folders of normpipe < 0.4 ----------------------------------------------

#' Name of a cached result of normpipe < 0.4 in the current layout
#' @noRd
migrate_cache_name <- function(n) {
  n <- sub("^p0_reference$", "prepare_reference", n)
  n <- sub("^m1_model_", "calibrate_model_", n)
  n <- sub("^m1_loadings_", "calibrate_loadings_", n)
  n <- sub("^m3a_lodo_", "norms_crosscheck_", n)
  n <- sub("^m3a_norms_", "norms_", n)
  n <- sub("^m3b_continuous_", "continuous_norms_", n)
  n <- sub("^m3c_eap_", "intervals_conditional_", n)
  n <- sub("^m3d_integrated_", "intervals_integrated_", n)
  n <- sub("^m2_donor_frame$", "transfer_donor_frame", n)
  n <- sub("^(transfer_.*)_lodo((__[A-Za-z_]+)|(_selection))?$", "\\1_crosscheck\\2", n)
  n
}

#' Move an output folder of normpipe < 0.4 to the current layout
#'
#' Since normpipe 0.4 the stages have descriptive names, the cached
#' intermediate results live in `out_dir/cache`, and the result files are
#' named after their content. `migrate_output()` renames the files of an
#' output folder written by an earlier version, so that a rerun finds every
#' cached result and nothing is computed again:
#'
#' * cached results (`.rds`) move to `cache/` under their new names, e.g.
#'   `m1_model_<outcome>` to `calibrate_model_<outcome>`, `m3a_norms_` to
#'   `norms_`, `m3b_continuous_` to `continuous_norms_`, `m3c_eap_` to
#'   `intervals_conditional_`, `m3d_integrated_` to `intervals_integrated_`,
#'   `transfer_<outcome>_lodo` to `transfer_<outcome>_crosscheck`;
#' * `cache_manifest.csv` and `superseded/` move to `cache/`;
#' * result files are renamed: `norm_tables_long.csv` to `norm_tables.csv`,
#'   `norm_tables_age_continuous.csv` to `norm_tables_continuous.csv`,
#'   `p0_scale_registry.csv` to `scale_registry.csv`, `x1_reliability.csv`
#'   to `reliability.csv`, `m2_family_selection.csv` to
#'   `transfer_family_selection.csv`, `c1_crosscheck.csv` to
#'   `transfer_crosscheck.csv`, `p0_reference_scores.rds` to
#'   `reference_scores.rds`, `p0_rake_diagnostics.rds` to
#'   `rake_diagnostics.rds`, the folders `platform/` to `compact/` and
#'   `portal/` to `json/`.
#'
#' The moved results are marked as adopted in the manifest: the next run
#' accepts each of them once as current and records its key under the
#' current version. Migrate only folders whose results are up to date with
#' the configuration and data you will rerun them with; a change made at the
#' same time as the migration is not detected for adopted results. The
#' content of the exported tables is not changed; rerun the `export` stage
#' to write them in the current format.
#'
#' @param out_dir Output folder written by normpipe < 0.4.
#' @param dry_run If `TRUE`, nothing is moved; the planned moves are
#'   returned.
#' @return Tibble with `from` and `to` (paths relative to `out_dir`),
#'   invisibly unless `dry_run`.
#' @export
migrate_output <- function(out_dir, dry_run = FALSE) {
  stopifnot(dir.exists(out_dir))
  results <- c("norm_tables_long.csv" = "norm_tables.csv",
               "norm_tables_age_continuous.csv" = "norm_tables_continuous.csv",
               "p0_scale_registry.csv" = "scale_registry.csv",
               "x1_reliability.csv" = "reliability.csv",
               "m2_family_selection.csv" = "transfer_family_selection.csv",
               "c1_crosscheck.csv" = "transfer_crosscheck.csv",
               "p0_reference_scores.rds" = "reference_scores.rds",
               "p0_rake_diagnostics.rds" = "rake_diagnostics.rds",
               "platform" = "compact", "portal" = "json",
               "cache_manifest.csv" = file.path("cache", "cache_manifest.csv"),
               "superseded" = file.path("cache", "superseded"))
  top <- list.files(out_dir)
  plan <- tibble(from = intersect(names(results), top))
  plan$to <- unname(results[plan$from])
  rds_files <- setdiff(grep("\\.rds$", top, value = TRUE), names(results))
  rds_files <- rds_files[!file.info(file.path(out_dir, rds_files))$isdir]
  plan <- bind_rows(plan, tibble(from = rds_files,
                                 to = file.path("cache", paste0(migrate_cache_name(sub("\\.rds$", "", rds_files)), ".rds"))))
  if ("platform" %in% top) {
    age <- list.files(file.path(out_dir, "platform"), pattern = "_age\\.rds$")
    plan <- bind_rows(tibble(from = file.path("platform", age),
                             to = file.path("platform", sub("_age\\.rds$", "_continuous.rds", age))), plan)
  }
  clash <- plan$to[file.exists(file.path(out_dir, plan$to))]
  if (length(clash))
    stop("migrate_output(): target(s) already exist in '", out_dir, "': ",
         paste(utils::head(clash, 5), collapse = ", "), ". Has the folder been migrated already?")
  if (dry_run) return(plan)
  dir.create(cache_dir(out_dir), showWarnings = FALSE, recursive = TRUE)
  ok <- file.rename(file.path(out_dir, plan$from), file.path(out_dir, plan$to))
  if (!all(ok)) warning("migrate_output(): could not move ", paste(plan$from[!ok], collapse = ", "))
  man <- read_manifest(out_dir)
  if (nrow(man)) {
    man$name <- migrate_cache_name(man$name)
    man$adopted <- "TRUE"
    write_manifest(out_dir, man)
  }
  invisible(plan)
}
