demo_inst <- function() read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
pid5_inst <- function() read_instrument(system.file("instruments/pid5/pid5.yml", package = "normpipe"))

test_that("the factor structure follows the instrument", {
  st <- invariance_structure(pid5_inst())
  expect_length(st$first, 17)                       # the bridge blocks
  expect_true(all(lengths(st$first) == 2))
  expect_setequal(names(st$second), c("negative_affectivity", "detachment", "antagonism",
                                      "disinhibition", "psychoticism"))
  expect_setequal(setdiff(names(st$first), unlist(st$second)), c("perseveration", "rigid_perfect"))
  inst <- demo_inst()
  expect_error(invariance_structure(inst), "bridge")
  st <- invariance_structure(inst, c("full_a", "full_b"))
  expect_equal(st$first, list(a = 1:4, b = 5:8)); expect_length(st$second, 0)
  st <- invariance_structure(inst, "full_ab")
  expect_equal(st$second, list(ab = c("a", "b")))
  expect_error(invariance_structure(inst, c("full_a", "full_total")), "more than one")
})

test_that("the syntax fixes and frees what the levels say", {
  st <- list(first = list(a = c("x1", "x2"), b = c("x3", "x4"), c = c("x5", "x6"), z = c("x7", "x8")),
             second = list(d = c("a", "b", "c")))
  nthr <- setNames(c(3, 3, 3, 3, 3, 3, 3, 1), paste0("x", 1:8))
  syn <- function(...) strsplit(invariance_syntax(st, nthr, 2, ...), "\n")[[1]]
  cfg <- syn("configural")
  expect_true("x2 | c(NA, NA)*t1 + c(NA, NA)*t2 + c(NA, NA)*t3" %in% cfg)
  expect_true(all(c("x2 ~~ c(1, 1)*x2", "x2 ~ 0*1", "F_a ~ c(0, 0)*1", "D_d ~ c(0, 0)*1") %in% cfg))
  sca <- syn("scalar")
  expect_true(all(c("F_a =~ 1*x1 + c(l_x2, l_x2)*x2", "x2 ~~ c(1, NA)*x2", "F_a ~ c(0, NA)*1",
                    "x8 ~~ c(1, 1)*x8",             # two categories: residual variance stays 1
                    "D_d =~ 1*F_a + c(NA, NA)*F_b + c(NA, NA)*F_c", "D_d ~ c(0, 0)*1") %in% sca))
  expect_true("x2 ~~ c(1, 1)*x2" %in% syn("strict"))
  so <- syn("second_order")
  expect_true(all(c("F_a ~ c(0, 0)*1", "F_z ~ c(0, NA)*1", "D_d ~ c(0, NA)*1",
                    "D_d =~ 1*F_a + c(g_F_b, g_F_b)*F_b + c(g_F_c, g_F_c)*F_c") %in% so))
  dr <- syn("configural", "domain_residual")
  expect_true(all(c("x1 ~~ x2", "F_z =~ 1*x7 + c(NA, NA)*x8") %in% dr))
  expect_true(any(grepl("^D_d =~ 1\\*x1 \\+ .*x6$", dr)))
  expect_error(invariance_syntax(st, nthr, 2, "scalar", "domain_residual"), "configural")
})

test_that("the Chen rule needs the CFI and one of RMSEA or SRMR", {
  tab <- tibble::tibble(model = c("configural", "scalar", "strict"),
                        compared_with = c(NA, "configural", "scalar"),
                        cfi.scaled = c(.950, .940, .925), rmsea.scaled = c(.050, .052, .053),
                        srmr = c(.050, .060, .061))
  fl <- invariance_flags(tab)
  expect_equal(fl$d_cfi, c(NA, -.010, -.015))
  expect_equal(fl$flag, c(NA, TRUE, FALSE))       # scalar: CFI and SRMR; strict: CFI only
  expect_false(invariance_flags(tab, srmr = .030)$flag[2])
  d <- invariance_decision(dplyr::mutate(fl, grouping = "g"))
  expect_false(d$ok); expect_equal(nrow(d$flagged), 1)
  expect_true(any(grepl("FLAGGED", d$text)))
})

test_that("sparse categories are merged in all groups", {
  x <- c(rep(0:3, each = 20), 3, 3); g <- rep(c("a", "b"), length.out = 82)
  expect_identical(collapse_sparse(x, g, 5), x)
  x2 <- c(rep(0:2, each = 20), 3, 3, 3); g2 <- c(rep(c("a", "b"), 30), "a", "a", "b")
  expect_equal(max(collapse_sparse(x2, g2, 5)), 2)
  expect_identical(collapse_sparse(x2, g2, 0), x2)
  # merged once, even when the merged category is still sparse; empty ones always
  x3 <- c(rep(0:1, each = 20), 2, 2, 2, 2, 3, 3); g3 <- c(rep(c("a", "b"), 20), "a", "b", "a", "b", "a", "b")
  expect_equal(sort(unique(collapse_sparse(x3, g3, 5))), 0:2)
  x4 <- c(rep(0:2, each = 20), 3, 3); g4 <- c(rep(c("a", "b"), 30), "a", "a")
  expect_equal(max(collapse_sparse(x4, g4, 0)), 2)
})

test_that("the sequence runs and finds invariant items invariant", {
  skip_if_not_installed("lavaan")
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 1500, seed = 3)
  it <- sim$items; it[5, 2] <- NA
  st <- invariance_structure(inst, c("full_a", "full_b"))
  res <- invariance_sequence(it, sim$meta$gender, inst, st, label = "gender")
  expect_equal(res$model, c("configural", "scalar", "strict"))
  expect_true(all(res$converged)); expect_equal(res$n, rep(1500L, 3))
  # 2 groups, 8 items, 6 thresholds: configural 2 * 19, then + 6 loadings + 48
  # thresholds - 2 means - 8 residual variances, then + 8 residual variances
  expect_equal(res$df.scaled, c(38, 82, 90))
  expect_equal(res$flag, c(NA, FALSE, FALSE))
  expect_gt(min(res$d_cfi[-1]), -0.01)
  lw <- invariance_sequence(it, sim$meta$gender, inst, st, levels = "configural", missing = "listwise")
  expect_equal(lw$n, 1499L)
  # a shifted item in one group is noticed by the scalar model
  bad <- sim$items; f <- sim$meta$gender == "f"
  bad$it2[f] <- pmin(7, bad$it2[f] + 2L)
  rb <- invariance_sequence(bad, sim$meta$gender, inst, st, levels = c("configural", "scalar"))
  expect_lt(rb$d_cfi[2], res$d_cfi[2] - 0.01)
  sc <- invariance_screen(sim$items, sim$meta$gender, inst, scales = "full_a")
  expect_equal(unique(sc$outcome), "full_a"); expect_equal(nrow(sc), 3)
})

test_that("run_invariance writes aggregate files only", {
  skip_if_not_installed("lavaan")
  inst <- demo_inst()
  ref <- simulate_reference(inst, N = 900, seed = 1)
  don <- simulate_reference(inst, N = 900, seed = 2)
  don$meta$sample <- rep(c("d1", "d2"), length.out = 900)
  dir <- tempfile()
  out <- run_invariance(ref[c("meta", "items")], don[c("meta", "items")], inst, out_dir = dir,
                        structure = invariance_structure(inst, c("full_a", "full_b")))
  expect_setequal(unique(out$fit$grouping), c("samples", "samples_listwise", "gender", "age_group"))
  expect_setequal(list.files(dir), c("invariance_fit.csv", "invariance_screen.csv",
                                     "invariance_decision.txt"))
  expect_false(any(c("id", "age") %in% names(out$fit)))
  expect_true(is.logical(out$decision$ok))
})
