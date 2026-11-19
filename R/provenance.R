# ==============================================================================
# provenance.R -- what a set of result tables was computed with
#
# The export stage writes run_info.json next to the tables: package and
# format version, the source of the package (commit when installed from
# GitHub), R and the versions of the packages that compute, a checksum of
# the instrument file and a hash of the configuration as read, the engine
# of the item models, the seed rule, the settings of `uncertainty`, the
# raking, the size of the reference sample and the routes of the scales.
# ==============================================================================

#' Versions of the packages a run may compute with (installed ones only)
#' @noRd
package_versions <- function() {
  pk <- c("normpipe", "mirt", "brms", "rstan", "cmdstanr", "loo", "survey", "statmod", "numDeriv",
          "dplyr", "tidyr", "purrr", "rlang", "jsonlite", "lavaan")
  out <- list()
  for (p in pk) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NULL)
    if (!is.null(v)) out[[p]] <- v
  }
  out
}

#' Commit of the installed package, when it was installed from GitHub
#' @noRd
package_commit <- function() {
  d <- tryCatch(utils::packageDescription("normpipe"), error = function(e) NULL)
  if (is.null(d) || !is.list(d)) return(NULL)
  d$RemoteSha %||% d$GithubSHA1
}

#' Write run_info.json into the output folder of a run
#' @noRd
write_run_info <- function(ctx) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    ctx$say("run_info.json not written: jsonlite is not installed")
    return(invisible(NULL))
  }
  inst <- ctx$inst
  src <- inst$source
  md5 <- if (is.character(src) && !is.na(src) && file.exists(src)) unname(tools::md5sum(src))
  dg <- ctx$rake_diagnostics
  reg <- ctx$reg
  info <- list(
    normpipe_version = as.character(utils::packageVersion("normpipe")),
    format_version = FORMAT_VERSION,
    commit = package_commit(),
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    instrument = list(id = inst$id, file = if (is.character(src) && !is.na(src)) basename(src),
                      file_md5 = md5,
                      configuration_hash = rlang::hash(unclass(inst)[setdiff(names(inst), c("source", "donors_data"))])),
    r_version = R.version.string,
    platform = R.version$platform,
    packages = package_versions(),
    measurement_engine = resolved_engine(inst),
    transfer_backend = if (isTRUE(inst$transfer$enabled)) transfer_backend(inst),
    seed = inst$uncertainty$seed,
    seed_rule = if (is.null(inst$uncertainty$seed)) "per step, derived from the names of the step and the scale" else
      "per step, derived from uncertainty$seed and the names of the step and the scale",
    uncertainty = inst$uncertainty,
    weighting = list(raked = !is.null(dg),
                     kish_neff_design = if (!is.null(dg)) dg$kish_neff_design,
                     kish_neff_raked = if (!is.null(dg)) dg$kish_neff_raked),
    reference_n = if (!is.null(ctx$ref_scores)) nrow(ctx$ref_scores),
    scales = if (!is.null(reg)) as.list(table(reg$route)),
    stages = as.list(union(ctx$done, "export")))
  info <- Filter(Negate(is.null), info)
  jsonlite::write_json(info, file.path(ctx$out_dir, "run_info.json"), auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null")
  invisible(file.path(ctx$out_dir, "run_info.json"))
}
