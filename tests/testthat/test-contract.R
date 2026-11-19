demo_inst <- function() read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))

# --- input -----------------------------------------------------------------------

test_that("item data that do not belong to the rows of meta are an error", {
  inst <- demo_inst(); inst$uncertainty$n_boot <- 2
  sim <- simulate_reference(inst, N = 120, seed = 3)
  out <- file.path(tempdir(), "np_align"); unlink(out, recursive = TRUE)
  expect_error(run_norming(inst, list(meta = sim$meta, items = sim$items[-1, ]), out, stages = "prepare",
                           verbose = FALSE), "120 in|rows, `meta` has 120")
  m2 <- sim$meta; m2$id[2] <- m2$id[1]
  expect_error(run_norming(inst, list(meta = m2, items = sim$items), out, stages = "prepare", verbose = FALSE),
               "not unique")
  # with an id column, the rows of items are matched to meta
  it <- cbind(id = sim$meta$id, sim$items)[sample(nrow(sim$items)), ]
  a <- align_sample(list(meta = sim$meta, items = sim$items), "ref")
  b <- align_sample(list(meta = sim$meta, items = it), "ref")
  expect_equal(as.data.frame(b$items), as.data.frame(a$items), ignore_attr = TRUE)
  expect_false("id" %in% names(b$items))
  it$id[1] <- "nobody"
  expect_error(align_sample(list(meta = sim$meta, items = it), "ref"), "not the ids")
  # an aligned sample is passed on unchanged (its cache keys stay the same)
  expect_identical(align_sample(sim, "ref"), sim)
})

# --- configuration ----------------------------------------------------------------

test_that("unknown configuration keys are an error, with a hint for typos", {
  cfg <- yaml::read_yaml(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
  bad <- cfg; bad$norm_scal <- "IQ"
  expect_error(read_instrument(bad), "`norm_scale` instead of `norm_scal`")
  bad <- cfg; bad$uncertainty$n_boots <- 1000
  expect_error(read_instrument(bad), "n_boots")
  old <- options(normpipe.unknown_keys = "warn")
  expect_warning(read_instrument(bad), "n_boots")
  options(old)
  bad <- cfg; bad$items$revers <- c(3, 7)
  expect_error(read_instrument(bad), "revers")
  bad <- cfg; bad$measurement$engin <- "mirt"
  expect_error(read_instrument(bad), "engin")
  # free keys of a scale are kept as metadata; a near miss of a known key warns
  ok <- cfg; ok$versions[[1]]$scales[[1]]$description <- "first scale"
  expect_silent(i <- read_instrument(ok))
  expect_equal(get_scale(i, "full_a")$extra$description, "first scale")
  bad <- cfg; bad$versions[[1]]$scales[[1]]$metrik <- "sum"
  expect_warning(read_instrument(bad), "`metric` instead of `metrik`")
  # notes and requires
  ok <- cfg; ok$notes <- "anything"; ok$requires <- "1.0"
  expect_equal(read_instrument(ok)$notes, "anything")
  bad <- cfg; bad$requires <- ">= 99.0"
  expect_error(read_instrument(bad), "needs normpipe >= 99.0")
  # the shipped configurations read without any message
  for (f in list.files(system.file("instruments", package = "normpipe"), "\\.yml$", recursive = TRUE,
                       full.names = TRUE))
    expect_silent(read_instrument(f))
})

test_that("the defaults are fixed", {
  inst <- read_instrument(list(id = "x", items = list(n_items = 2), versions = list(list(id = "v", scales = list(list(id = "s", items = 1:2))))))
  expect_identical(inst$uncertainty, list(prior_mass = 1, n_boot = 200, n_param_draws = 50, n_draws = 500,
                                          pctl_clamp = 1e-4, probs = c(.025, .975),
                                          pctl_convention = "midpoint", n_mc = 1500L))
  expect_identical(inst$weighting, list(trim_quantile = 0.99, rerake = TRUE, incomplete = "drop"))
  expect_identical(inst$measurement, list(source = "fit", engine = "auto", calibration = "reference"))
  expect_identical(inst$export, list(formats = c("long", "compact", "json"), language = "en"))
  expect_identical(inst$norm_scale[c("name", "mean", "sd")], list(name = "T", mean = 50, sd = 10))
  expect_identical(build_groups(list(continuous = list(variable = "age", range = c(18, 20))))$continuous$kernel,
                   list(half_width0 = 4, min_neff = 150, max_half_width = 12))
  s <- get_scale(inst, "v_s")
  expect_identical(s[c("metric", "model", "level", "min_answered")],
                   list(metric = "mean", model = "grm", level = "scale", min_answered = 1))
  expect_identical(missing_cfg(inst)$method, "prorate")
})

test_that("no exported function disappears", {
  # names exported by normpipe 1.0.0; later versions may add to them
  v1 <- c("add_age_group", "add_groups", "answered_share", "bridge_scores", "cache_get", "cache_key",
          "cache_status", "cell_vars", "conditional_tables", "crosscheck_table", "direct_frame",
          "draw_model_params", "export_compact", "export_json", "export_platform_tables",
          "export_portal_json", "fit_scale_model", "fit_transfer", "frequency_frame", "get_scale",
          "group_sources", "impute_items", "integrated_intervals", "invariance_decision",
          "invariance_flags", "invariance_screen", "invariance_sequence", "invariance_structure",
          "invariance_syntax", "item_indices", "item_values", "long_table_scale", "make_blocks",
          "make_margins", "marginal_reliability", "migrate_output", "norming_context",
          "norms_all_cells", "project_outcome", "rake_to_margins", "read_instrument", "run_invariance",
          "run_norming", "run_stage", "scale_registry", "scale_routes", "score_scales",
          "select_transfer_family", "simulate_reference", "snap_to_grid", "transfer_data",
          "transfer_donor_frame", "transfer_model", "transfer_model_multiple", "transfer_spec",
          "weighted_norms", "weighted_norms_continuous")
  expect_length(setdiff(v1, getNamespaceExports("normpipe")), 0)
  expect_equal(STAGES, c("prepare", "calibrate", "transfer", "norms", "continuous_norms",
                         "intervals_conditional", "intervals_integrated", "export"))
})

# --- random numbers ----------------------------------------------------------------

test_that("results do not depend on the random-number state or on the cache", {
  inst <- demo_inst()
  inst$uncertainty$n_boot <- 8; inst$uncertainty$n_param_draws <- 2; inst$uncertainty$n_mc <- 100L
  sim <- simulate_reference(inst, N = 400, seed = 21)
  st <- c("calibrate", "norms", "intervals", "export")
  o <- c("full_a", "full_b", "full_total")
  run <- function(out, rng) {
    unlink(out, recursive = TRUE)
    set.seed(rng)
    before <- .Random.seed
    long <- run_norming(inst, sim, out, stages = st, outcomes = o, verbose = FALSE)
    expect_identical(.Random.seed, before)                # the caller's state is untouched
    long
  }
  a <- run(file.path(tempdir(), "np_rng_a"), 1)
  b <- run(file.path(tempdir(), "np_rng_b"), 2)
  expect_identical(a, b)
  # frames from the cache, intervals recomputed: the same intervals
  out <- file.path(tempdir(), "np_rng_b")
  unlink(file.path(out, "cache", paste0("intervals_integrated_", o, ".rds")))
  set.seed(3); runif(10)
  c2 <- run_norming(inst, sim, out, stages = st, outcomes = o, verbose = FALSE)
  expect_identical(c2, a)
  # another seed gives other draws
  inst$uncertainty$seed <- 7
  d <- run_norming(inst, sim, file.path(tempdir(), "np_rng_c"), stages = st, outcomes = o, verbose = FALSE)
  expect_false(identical(d$T_lo, a$T_lo))
  expect_identical(with_seed(1, runif(1)), with_seed(1, runif(1)))
  expect_equal(seed_of("a", "b"), seed_of("a", "b")); expect_false(seed_of("ab") == seed_of("ba"))
})

# --- cache keys --------------------------------------------------------------------

test_that("labels, the kernel and new options at their default recompute nothing", {
  inst <- demo_inst(); inst$uncertainty$n_boot <- 4
  sim <- simulate_reference(inst, N = 200, seed = 6)
  out <- file.path(tempdir(), "np_keys"); unlink(out, recursive = TRUE)
  st <- c("norms", "continuous_norms")
  run_norming(inst, sim, out, stages = st, outcomes = c("full_a", "full_b"), verbose = FALSE)
  m1 <- cache_status(out)
  inst$versions$full$scales$b$label <- "another label"
  inst$versions$full$label <- "another version label"
  inst$uncertainty$n_mc <- 1500L
  run_norming(inst, sim, out, stages = st, outcomes = c("full_a", "full_b"), verbose = FALSE)
  expect_identical(cache_status(out)$key, m1$key)
  # the kernel changes the continuous norms only
  inst$groups$continuous$kernel$min_neff <- 40
  run_norming(inst, sim, out, stages = st, outcomes = c("full_a", "full_b"), verbose = FALSE)
  m2 <- cache_status(out)
  changed <- m2$name[m2$key != m1$key[match(m2$name, m1$name)]]
  expect_setequal(changed, c("continuous_norms_full_a", "continuous_norms_full_b"))
  # a new scale leaves the frames of the others alone
  inst2 <- inst
  inst2$versions$full$scales$b2 <- inst$versions$full$scales$b
  inst2$versions$full$scales$b2$id <- "b2"
  run_norming(inst2, sim, out, stages = st, outcomes = c("full_a", "full_b", "full_b2"), verbose = FALSE)
  m3 <- cache_status(out)
  keep <- c("frame_full_a", "frame_full_b", "norms_full_a", "norms_full_b")
  expect_identical(m3$key[match(keep, m3$name)], m2$key[match(keep, m2$name)])
})

test_that("results recorded under an earlier key format are taken over", {
  out <- file.path(tempdir(), "np_legacy"); unlink(out, recursive = TRUE)
  cache_get(out, "r", "old", function() 1)
  expect_equal(cache_get(out, "r", "new", function() stop("recomputed"), legacy_key = "old"), 1)
  expect_equal(cache_status(out)$key, "new")
  expect_equal(cache_get(out, "r", "newer", function() 2, legacy_key = "old"), 2)
  # accept(): a cached result that no longer fits is recomputed
  msg <- character(0)
  expect_equal(cache_get(out, "r", "newer", function() 3, accept = function(x) "too old",
                         say = \(...) msg <<- paste0(...)), 3)
  expect_match(msg, "too old")
})

test_that("an item model fitted with another engine is refitted", {
  inst <- demo_inst(); inst$measurement$engine <- "internal"
  sim <- simulate_reference(inst, N = 200, seed = 7)
  out <- file.path(tempdir(), "np_engine"); unlink(out, recursive = TRUE)
  run_norming(inst, sim, out, stages = "calibrate", outcomes = "full_b", verbose = FALSE)
  f <- file.path(out, "cache", "calibrate_model_full_b.rds")
  m <- readRDS(f); m$engine <- "mirt"; saveRDS(m, f)
  msg <- capture_messages(run_norming(inst, sim, out, stages = "calibrate", outcomes = "full_b"))
  expect_true(any(grepl("fitted with engine mirt", msg)))
  expect_equal(readRDS(f)$engine, "internal")
})

# --- outputs -------------------------------------------------------------------------

test_that("tables carry whole-number summed scores and the run is documented", {
  skip_if_not_installed("jsonlite")
  inst <- demo_inst(); inst$uncertainty$n_boot <- 3
  sim <- simulate_reference(inst, N = 150, seed = 8)
  out <- file.path(tempdir(), "np_info"); unlink(out, recursive = TRUE)
  long <- run_norming(inst, sim, out, stages = c("norms", "continuous_norms", "export"), verbose = FALSE)
  expect_equal(names(long), c("version", "scale", "level", "route", "cell", "gender", "age_group", "raw",
                              "pctl_med", "pctl_lo", "pctl_hi", "T", "T_lo", "T_hi", "T_cond_lo",
                              "T_cond_hi", "T_int_lo", "T_int_hi", "extrapolation", "n_eff", "sum_score"))
  a <- long[long$scale == "a" & long$version == "full", ]
  expect_equal(a$sum_score, as.integer(round(a$raw * 4)))
  b <- long[long$scale == "b", ]
  expect_equal(b$sum_score, as.integer(b$raw))
  ab <- long[long$scale == "ab", ]                    # mean of two equally long means
  expect_equal(ab$sum_score, as.integer(round(ab$raw * 8)))
  expect_true(all(is.na(sum_score_of(c(1, 1.01), list(dense = TRUE)))))   # reporting grid
  cmp <- read.csv(file.path(out, "compact", "norms_full.csv"))
  expect_equal(tail(names(cmp), 1), "sum_score")
  ct <- read.csv(file.path(out, "norm_tables_continuous.csv"))
  expect_true("sum_score" %in% names(ct))
  js <- jsonlite::read_json(file.path(out, "json", "demo7.json"))
  expect_equal(js$normpipe_version, as.character(utils::packageVersion("normpipe")))
  expect_equal(js$format_version, 1L)
  ri <- jsonlite::read_json(file.path(out, "run_info.json"))
  expect_equal(ri$normpipe_version, as.character(utils::packageVersion("normpipe")))
  expect_true(all(c("format_version", "instrument", "r_version", "packages", "measurement_engine",
                    "seed_rule", "uncertainty", "weighting", "reference_n", "scales", "stages") %in% names(ri)))
  expect_equal(ri$reference_n, 150L)
  expect_match(ri$instrument$file_md5, "^[0-9a-f]{32}$")
})

# --- reference values ------------------------------------------------------------------

test_that("the demo instrument gives the stored reference norms", {
  # A change of these numbers is a change of results: it needs a NEWS entry
  # under "Changes to results" and a bump of CACHE_SCHEMA for the stage.
  # Regenerate with tests/testthat/_ref/make_reference.R.
  ref_file <- test_path("_ref", "demo7_reference.csv")
  skip_if_not(file.exists(ref_file))
  long <- reference_run()
  ref <- utils::read.csv(ref_file, stringsAsFactors = FALSE)
  expect_equal(names(long), names(ref))
  num <- names(ref)[vapply(ref, is.numeric, logical(1))]
  for (v in setdiff(names(ref), num)) expect_equal(as.character(long[[v]]), as.character(ref[[v]]), info = v)
  for (v in num) expect_equal(as.numeric(long[[v]]), ref[[v]], tolerance = 1e-5, info = v)
})
