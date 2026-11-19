# Writes the example output of the demo instrument to inst/examples/ (the
# tables the getting-started vignette reads). Run from the package folder:
#   Rscript scripts/make_examples.R
pkgload::load_all(".", quiet = TRUE)
inst <- read_instrument("inst/instruments/demo7/demo7.yml")
sim <- simulate_reference(inst, N = 1500, seed = 2)
out <- file.path(tempdir(), "normpipe-examples")
unlink(out, recursive = TRUE)
run_norming(inst, ref = list(meta = sim$meta, items = sim$items), out_dir = out,
            stages = c("calibrate", "norms", "continuous_norms", "intervals", "export"))
dir.create("inst/examples", showWarnings = FALSE)
file.copy(file.path(out, c("norm_tables.csv", "reliability.csv", "run_info.json")), "inst/examples",
          overwrite = TRUE)
file.copy(file.path(out, "json", "demo7.json"), "inst/examples", overwrite = TRUE)
