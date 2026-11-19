# The reference run of the demo instrument (test-contract.R): small enough
# for the test suite, and with every random step and route of the direct
# norming. The internal estimator is used so that the numbers do not depend
# on whether mirt is installed.
reference_run <- function(out = file.path(tempdir(), "np_reference")) {
  inst <- read_instrument(system.file("instruments/demo7/demo7.yml", package = "normpipe"))
  inst$measurement$engine <- "internal"
  inst$uncertainty$n_boot <- 20; inst$uncertainty$n_param_draws <- 3; inst$uncertainty$n_mc <- 200L
  sim <- simulate_reference(inst, N = 400, seed = 21)
  pop <- expand.grid(age = 18:80, gender = c("m", "f")); pop$Freq <- 1000 + 5 * pop$age
  unlink(out, recursive = TRUE)
  long <- run_norming(inst, list(meta = sim$meta, items = sim$items), out,
                      margins = make_margins(pop, inst),
                      stages = c("calibrate", "norms", "intervals", "export"),
                      outcomes = c("full_a", "full_b", "full_total", "short_a"), verbose = FALSE)
  as.data.frame(long)
}
