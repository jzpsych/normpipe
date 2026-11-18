demo_inst <- function() read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))

test_that("stage names of earlier versions and `intervals` resolve to the current stages", {
  expect_equal(resolve_stages(c("X1", "P0", "M3a")), c("prepare", "norms", "export"))
  expect_equal(resolve_stages(c("M1", "M2", "M3b", "M3c", "M3d")), STAGES[c(2, 3, 5, 6, 7)])
  expect_equal(resolve_stages("intervals"), c("intervals_conditional", "intervals_integrated"))
  expect_equal(resolve_stages(STAGES), STAGES)
  expect_error(resolve_stages("M9"), "Unknown stage")
})

test_that("the `cells` block of earlier versions gives the same groups", {
  g_old <- build_groups(cells = list(gender = c("m", "f"), age_breaks = list(18, 40, 65, "Inf"),
                                     ages = c(18, 80), kernel = list(min_neff = 60)))
  g_new <- build_groups(list(by = list(gender = c("m", "f"), age = list(breaks = list(18, 40, 65, "Inf"))),
                             continuous = list(variable = "age", range = c(18, 80),
                                               kernel = list(min_neff = 60))))
  expect_equal(g_old, g_new)
  expect_equal(demo_inst()$groups, g_new)
  expect_error(build_groups(list(), list()), "both")
  expect_error(build_groups(list(by = list(age = list(breaks = 5)))), "increasing")
  expect_error(build_groups(list(by = list(gender = c("m", "f")),
                                 continuous = list(variable = "age", range = c(18, 80), by = "region"))),
               "not a categorical variable")
})

test_that("norm groups can be built from any variables", {
  inst <- demo_inst()
  inst$groups <- build_groups(list(
    by = list(education = c("low", "mid", "high"), tenure = list(breaks = c(0, 5, 20))),
    continuous = list(variable = "tenure", range = c(0, 19), by = character(0),
                      kernel = list(min_neff = 40))))
  expect_equal(cell_vars(inst), c("education", "tenure_group"))
  expect_equal(group_sources(inst), c("tenure", "education"))
  sim <- simulate_reference(inst, N = 400, seed = 4)
  expect_true(all(c("education", "tenure") %in% names(sim$meta)))
  expect_true(all(sim$meta$tenure >= 0 & sim$meta$tenure <= 20))
  sim$meta$education[1:3] <- "other"                      # not a configured group
  s <- get_scale(inst, "full_b")
  ref <- dplyr::mutate(sim$meta, score = score_scales(sim$items, inst)$full_b)
  fr <- direct_frame(ref, s, inst, n_boot = 3)
  nm <- norms_all_cells(fr, s, inst)
  expect_equal(names(nm$cells)[1:2], c("education", "tenure_group"))
  expect_equal(nrow(dplyr::distinct(nm$cells, education, tenure_group)), 6)
  expect_false("other" %in% nm$cells$education)
  expect_equal(nrow(nm$overall), length(s$grid))
  ct <- weighted_norms_continuous(fr, s, inst)
  expect_equal(names(ct)[1], "tenure")
  expect_setequal(ct$tenure, 0:19)
  # frequency input with the same groups
  freq <- dplyr::count(add_groups(ref, inst), education, tenure_group, score, wt = weight, name = "count")
  fr0 <- direct_frame(ref, s, inst, n_boot = 0, prior_mass = 0)
  fq <- frequency_frame(freq, s, inst, n_boot = 0, prior_mass = 0)
  expect_equal(norms_all_cells(fr0, s, inst)$cells$T, norms_all_cells(fq, s, inst)$cells$T, tolerance = 1e-8)
  # the whole pipeline
  out <- file.path(tempdir(), "np_groups"); unlink(out, recursive = TRUE)
  inst$uncertainty$n_boot <- 5
  long <- run_norming(inst, sim, out, stages = c("norms", "continuous_norms", "export"),
                      outcomes = "full_b", verbose = FALSE)
  expect_true(all(c("education", "tenure_group") %in% names(long)))
  expect_true("low [0,5)" %in% long$cell)
  expect_true(all(file.exists(file.path(out, c("norm_tables.csv", "norm_tables_continuous.csv",
                                               "scale_registry.csv", "compact/norms_full.csv")))))
})

test_that("an instrument without groups gets the overall table only", {
  inst <- demo_inst()
  inst$groups <- build_groups(NULL)
  expect_length(cell_vars(inst), 0)
  sim <- simulate_reference(inst, N = 200, seed = 4)
  out <- file.path(tempdir(), "np_nogroups"); unlink(out, recursive = TRUE)
  inst$uncertainty$n_boot <- 5
  long <- run_norming(inst, list(meta = sim$meta["id"], items = sim$items), out,
                      stages = c("norms", "continuous_norms", "export"), outcomes = "full_b", verbose = FALSE)
  expect_equal(unique(long$cell), "overall")
  expect_false(file.exists(file.path(out, "norm_tables_continuous.csv")))
})

test_that("the norm scale sets the standard scores and their column names", {
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 800, seed = 5)
  s <- get_scale(inst, "full_b")
  ref <- dplyr::mutate(sim$meta, score = score_scales(sim$items, inst)$full_b)
  fr <- direct_frame(ref, s, inst, n_boot = 4)
  nt <- weighted_norms(fr, s, inst, cells = character(0))
  inst$norm_scale <- build_norm_scale("IQ")
  ni <- weighted_norms(fr, s, inst, cells = character(0))
  expect_true(all(c("IQ", "IQ_lo", "IQ_hi") %in% names(ni)))
  expect_false(any(c("T", "T_lo") %in% names(ni)))
  expect_equal(ni$IQ, 100 + 15 * (nt$T - 50) / 10)
  expect_equal(ni$pctl_med, nt$pctl_med)
  inst$norm_scale <- build_norm_scale("stanine")
  n9 <- weighted_norms(fr, s, inst, cells = character(0))
  expect_true(all(n9$stanine >= 1 & n9$stanine <= 9))
  sc <- build_norm_scale(list(name = "SW", mean = 100, sd = 10))
  expect_equal(pctl_to_norm(0.5, scale = sc), 100)
  expect_equal(norm_grid(NORM_SCALES$T), seq(10, 90, by = 0.1))
  expect_error(build_norm_scale("W"), "not a preset")
  expect_error(build_norm_scale(list(name = "raw", mean = 0, sd = 1)), "column name")
  # through the pipeline, with intervals
  inst$norm_scale <- build_norm_scale("z"); inst$uncertainty$n_boot <- 5; inst$uncertainty$n_param_draws <- 3
  out <- file.path(tempdir(), "np_scale"); unlink(out, recursive = TRUE)
  long <- run_norming(inst, sim, out, stages = c("calibrate", "norms", "intervals", "export"),
                      outcomes = "full_b", verbose = FALSE)
  expect_true(all(c("z", "z_lo", "z_hi", "z_cond_lo", "z_cond_hi", "z_int_lo", "z_int_hi") %in% names(long)))
  expect_false(any(grepl("^T($|_)", names(long))))
  expect_true(all(abs(long$z) < 4.5))
  ov <- long[long$cell == "overall", ]
  expect_true(all(ov$z_int_lo <= ov$z_int_hi))
})

test_that("the JSON export follows export$language", {
  skip_if_not_installed("jsonlite")
  inst <- demo_inst(); inst$uncertainty$n_boot <- 5
  sim <- simulate_reference(inst, N = 200, seed = 5)
  out <- file.path(tempdir(), "np_json"); unlink(out, recursive = TRUE)
  long <- run_norming(inst, sim, out, stages = c("norms", "export"), outcomes = "full_b", verbose = FALSE)
  en <- readLines(file.path(out, "json", "demo7.json"), warn = FALSE)
  f <- tempfile(fileext = ".json")
  export_json(long, inst, f, language = "de")
  de <- readLines(f, warn = FALSE)
  expect_false(identical(en, de))
  g <- tempfile(fileext = ".json")
  export_portal_json(long, inst, g)                        # name in normpipe < 0.4: German keys
  expect_identical(readLines(g, warn = FALSE), de)
  inst$export <- build_export(list(formats = "long"))
  out2 <- file.path(tempdir(), "np_json2"); unlink(out2, recursive = TRUE)
  run_norming(inst, sim, out2, stages = c("norms", "export"), outcomes = "full_b", verbose = FALSE)
  expect_false(dir.exists(file.path(out2, "json"))); expect_false(dir.exists(file.path(out2, "compact")))
  expect_true(file.exists(file.path(out2, "norm_tables.csv")))
})

test_that("migrate_output() moves an output folder of an earlier version and keeps its results", {
  out <- file.path(tempdir(), "np_migrate"); unlink(out, recursive = TRUE)
  dir.create(file.path(out, "platform"), recursive = TRUE)
  for (n in c("p0_reference", "m1_model_full_a", "m3a_norms_full_a", "m3b_continuous_full_a", "m3c_eap_full_a",
              "m3d_integrated_full_a", "frame_full_a", "transfer_full_b", "transfer_full_a_lodo",
              "transfer_full_a_lodo__gaussian", "m3a_lodo_full_a", "m2_donor_frame", "p0_reference_scores"))
    saveRDS(n, file.path(out, paste0(n, ".rds")))
  saveRDS(1, file.path(out, "platform", "norms_full_age.rds"))
  writeLines("x", file.path(out, "norm_tables_long.csv")); writeLines("x", file.path(out, "c1_crosscheck.csv"))
  utils::write.csv(data.frame(name = c("m3a_norms_full_a", "transfer_full_a_lodo"), key = c("k1", "k2"),
                              version = "0.3.1", time = "t"), file.path(out, "cache_manifest.csv"), row.names = FALSE)
  plan <- migrate_output(out, dry_run = TRUE)
  expect_true(file.exists(file.path(out, "m3a_norms_full_a.rds")))            # nothing moved yet
  migrate_output(out)
  expect_setequal(list.files(out), c("cache", "compact", "norm_tables.csv", "transfer_crosscheck.csv",
                                     "reference_scores.rds"))
  expect_setequal(sub("\\.rds$", "", list.files(file.path(out, "cache"), pattern = "rds$")),
                  c("prepare_reference", "calibrate_model_full_a", "norms_full_a", "continuous_norms_full_a",
                    "intervals_conditional_full_a", "intervals_integrated_full_a", "frame_full_a",
                    "transfer_full_b", "transfer_full_a_crosscheck", "transfer_full_a_crosscheck__gaussian",
                    "norms_crosscheck_full_a", "transfer_donor_frame"))
  expect_true(file.exists(file.path(out, "compact", "norms_full_continuous.rds")))
  st <- cache_status(out)
  expect_setequal(st$name, c("norms_full_a", "transfer_full_a_crosscheck"))
  expect_true(all(st$adopted == "TRUE")); expect_equal(st$key[st$name == "norms_full_a"], "k1")
  # an adopted result is accepted once and then keyed
  expect_equal(cache_get(out, "norms_full_a", "new", function() stop("recomputed")), "m3a_norms_full_a")
  expect_equal(cache_status(out)$key[cache_status(out)$name == "norms_full_a"], "new")
  expect_equal(cache_get(out, "norms_full_a", "newer", function() "fresh"), "fresh")
  expect_error(migrate_output(out), NA)                                          # nothing left to move
})
