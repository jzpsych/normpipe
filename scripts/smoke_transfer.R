# End-to-end test of the transfer module (stage M2, brms): the reference
# sample answers only the bridge items, the four full-form scales are
# projected from simulated donors, and the projected norms are compared with
# the norms of the (simulated, here known) full responses of the reference
# sample. The true percentile rank uses the same convention as the norms
# (uncertainty$pctl_convention, midpoint by default). About 25 minutes with
# rstan on four cores.
suppressPackageStartupMessages({library(normpipe); library(dplyr)})
inst <- read_instrument(system.file("instruments/demo7/demo7_transfer.yml", package = "normpipe"))
conv <- inst$uncertainty$pctl_convention
don <- simulate_reference(inst, N = 800, seed = 3)
don$meta$sample <- rep(c("quota", "clinical"), length.out = 800)
inst$donors_data <- list(meta = don$meta, items = don$items)
ref <- simulate_reference(inst, N = 1000, seed = 4)
out <- "demo-transfer-out"
res <- run_norming(inst, ref = list(meta = ref$meta, items = ref$items[, c("it1", "it2", "it5", "it6")]),
                   out_dir = out, stages = c("P0", "M2", "M3a", "X1"))
sc <- score_scales(ref$items, inst)
d <- add_age_group(mutate(ref$meta, a = sc$full_a, total = sc$full_total), inst)
true_pr <- function(x, at) {
  below <- mean(x < at - 1e-9); at_or_below <- mean(x <= at + 1e-9)
  if (conv == "cumulative") at_or_below else (below + at_or_below) / 2
}
nt <- readr::read_csv(file.path(out, "norm_tables_long.csv"), show_col_types = FALSE) %>%
  filter(version == "full", scale %in% c("a", "total"), !is.na(gender))
cmp <- nt %>% rowwise() %>%
  mutate(p_true = { sel <- d$gender == gender & as.character(d$age_group) == age_group
                    true_pr(d[[scale]][sel], raw) }) %>% ungroup() %>%
  filter(p_true > .05, p_true < .95) %>%
  mutate(T_true = 50 + 10 * qnorm(p_true), inside = T_true >= T_lo & T_true <= T_hi)
cat("Percentile convention:", conv, "\n")
print(cmp %>% group_by(scale) %>% summarise(n = n(), mean_abs_dT = mean(abs(T - T_true)), coverage = mean(inside)))
