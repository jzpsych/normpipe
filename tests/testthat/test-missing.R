demo_inst <- function() read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))

test_that("skipped items are imputed within the response categories and reproducibly", {
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 300, seed = 4)
  expect_identical(impute_items(sim$items, inst, M = 2)[[1]], sim$items)   # complete data unchanged
  it <- sim$items
  set.seed(1); holes <- cbind(sample(300, 60), sample(8, 60, TRUE)); it[holes] <- NA
  it[7, ] <- NA                                                           # nothing answered
  imp <- impute_items(it, inst, M = 5, seed = 2)
  expect_length(imp, 5)
  expect_equal(sum(attr(imp, "n_imputed")), sum(is.na(it)))
  for (d in imp) {
    expect_false(anyNA(d))
    expect_true(all(as.matrix(d) %in% 1:7))
    expect_identical(d[!is.na(it)], it[!is.na(it)])                       # answered items stay
  }
  expect_identical(impute_items(it, inst, M = 5, seed = 2), imp)
  expect_false(identical(imp[[1]], imp[[2]]))
  expect_equal(answered_share(it, inst)[7], 0)
  expect_equal(answered_share(sim$items, inst), rep(1, 300))
})

test_that("the hot deck follows the answered items of the block and the covariates", {
  inst <- demo_inst()
  sim <- simulate_reference(inst, N = 600, seed = 4)
  # block of every item: the domain (items 1-8), the largest scale below the total
  expect_equal(imputation_blocks(inst, 1:8)[["3"]], 1:8)
  expect_equal(imputation_blocks(inst, c(1, 2, 5, 6))[["5"]], c(1L, 2L, 5L, 6L))  # no complete scale: all presented items
  expect_equal(imputation_blocks(inst, c(1, 2, 5, 6))[["1"]], c(1L, 2L))
  it <- sim$items; m <- rowMeans(it[, 2:4])
  hi <- which(m >= 5.5)[1:5]; lo <- which(m <= 2.5)[1:5]
  it$it1[c(hi, lo)] <- NA
  d <- do.call(cbind, lapply(impute_items(it, inst, M = 40), \(x) x$it1))
  expect_gt(mean(d[hi, ]), mean(d[lo, ]) + 1.5)
  # two skipped items of a block are copied from the same neighbour
  it2 <- sim$items; it2[1, c("it1", "it2")] <- NA
  imp <- impute_items(it2, inst, M = 30)
  pairs <- t(sapply(imp, \(x) unlist(x[1, c("it1", "it2")])))
  full <- paste(sim$items$it1, sim$items$it2)
  expect_true(all(paste(pairs[, 1], pairs[, 2]) %in% full))
  # a covariate separates otherwise identical candidates
  inst1 <- demo_inst()
  dat <- as.data.frame(matrix(4L, 40, 8, dimnames = list(NULL, paste0("it", 1:8))))
  g <- rep(c("m", "f"), each = 20); dat$it1 <- ifelse(g == "m", 2L, 6L); dat$it1[c(1, 21)] <- NA
  imp <- impute_items(dat, inst1, M = 20, covariates = data.frame(gender = g), k = 5)
  expect_true(all(sapply(imp, \(x) x$it1[1]) == 2) && all(sapply(imp, \(x) x$it1[21]) == 6))
})

test_that("direct frames take the scores of the completed data sets draw by draw", {
  inst <- demo_inst(); inst$versions$full$scales$b$min_answered <- 0.5
  s <- get_scale(inst, "full_b")
  sim <- simulate_reference(inst, N = 120, seed = 5)
  it <- sim$items; it$it5[1:10] <- NA
  sc <- score_scales(it, inst)
  imp <- impute_items(it, inst, M = 3)
  si <- sapply(imp, \(d) score_scales(d, inst)$full_b)
  ref <- dplyr::mutate(sim$meta, score = sc$full_b)
  ref$score[11] <- NA                                    # excluded respondent stays out
  f0 <- direct_frame(ref, s, inst, n_boot = 7)
  f1 <- direct_frame(ref, s, inst, n_boot = 7, scores_imp = si)
  expect_equal(f1$weight, f0$weight); expect_equal(f1$id, f0$id)
  id <- as.character(sim$meta$id)
  expect_false(id[11] %in% f1$id)
  for (r in 1:7) {
    x <- f1[f1$.draw == r, ]
    expect_equal(x$score[match(id[-11], x$id)], si[-11, (r - 1) %% 3 + 1])
  }
  complete <- f1$id %in% id[12:120]
  expect_equal(f1$score[complete], f0$score[complete])
  expect_equal(f1$score[!f1$id %in% id], f0$score[!f0$id %in% id])   # prior rows untouched
  expect_error(direct_frame(ref, s, inst, scores_imp = si[-1, ]))
})

test_that("run_norming() applies one inclusion rule and norms from completed data", {
  inst <- demo_inst(); inst$uncertainty$n_boot <- 10
  sim <- simulate_reference(inst, N = 300, seed = 6)
  set.seed(2); sim$items[cbind(sample(300, 40), sample(8, 40, TRUE))] <- NA
  sim$items[5, 1:5] <- NA; sim$items[6, ] <- NA             # below half of the items
  inst$missing <- list(method = "impute", n_imputations = 4)
  out <- file.path(tempdir(), "np_impute"); unlink(out, recursive = TRUE)
  expect_message(run_norming(inst, sim, out, stages = c("P0", "M3a"), outcomes = "full_b"),
                 "2 respondent\\(s\\) answered less than 50 %")
  p0 <- readRDS(file.path(out, "p0_reference.rds"))
  expect_length(p0$imputed, 4); expect_equal(nrow(p0$scores), 298)
  fr <- frame_long(readRDS(file.path(out, "frame_full_b.rds")))
  expect_false(anyNA(fr$score))
  expect_false(any(as.character(sim$meta$id[5:6]) %in% fr$id))
  expect_equal(sum(fr$.draw == 1 & fr$id %in% as.character(sim$meta$id)), 298)
  # the default keeps prorating with min_answered of each scale: here complete cases
  inst$missing <- list()
  out2 <- file.path(tempdir(), "np_prorate"); unlink(out2, recursive = TRUE)
  run_norming(inst, sim, out2, stages = c("P0", "M3a"), outcomes = "full_b", verbose = FALSE)
  fr2 <- frame_long(readRDS(file.path(out2, "frame_full_b.rds")))
  expect_lt(sum(fr2$.draw == 1 & fr2$id %in% as.character(sim$meta$id)), 298)
  inst$missing <- list(method = "listwise")
  expect_error(run_norming(inst, sim, file.path(tempdir(), "np_bad"), stages = "P0", verbose = FALSE), "prorate")
})

test_that("item models use the responses of respondents with skipped items", {
  inst <- demo_inst(); s <- get_scale(inst, "full_a")
  sim <- simulate_reference(inst, N = 400, seed = 8)
  it <- sim$items; it$it1[1:80] <- NA
  out <- file.path(tempdir(), "np_fiml"); unlink(out, recursive = TRUE)
  inst$measurement$engine <- "internal"
  expect_message(run_norming(inst, list(meta = sim$meta, items = it), out, stages = c("P0", "M1"), outcomes = "full_a"),
                 "400 respondents \\(80 with skipped items\\)")
  expect_s3_class(readRDS(file.path(out, "m1_model_full_a.rds")), "normpipe_model")
})
