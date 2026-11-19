# Regenerates the reference norms of test-contract.R. Run from the package
# folder after a deliberate change of results (NEWS: "Changes to results"):
#   Rscript tests/testthat/_ref/make_reference.R
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-reference.R")
long <- reference_run()
utils::write.csv(long, "tests/testthat/_ref/demo7_reference.csv", row.names = FALSE)
