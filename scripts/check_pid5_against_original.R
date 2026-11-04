# Check the PID-5 configuration of normpipe against the scoring keys and the
# outcome registry of the original pipeline (pid5_norming_pipeline.R, Zenodo
# 10.5281/zenodo.22280403). Needs no data: item keys, item order, reverse
# keying, routes under the GESIS bridge and scoring of random responses.
#
# Usage: Rscript scripts/check_pid5_against_original.R path/to/pid5_norming_pipeline.R
suppressPackageStartupMessages({library(dplyr); library(purrr); library(tibble); library(normpipe)})
args <- commandArgs(TRUE)
orig_file <- if (length(args)) args[1] else "pid5_norming_pipeline.R"
yml <- if (length(args) > 1) args[2] else system.file("instruments/pid5/pid5.yml", package = "normpipe")

# evaluate only the top-level assignments needed from the original file
keep <- c("reverse_items", "facet_key", "domain_key", "sf_facet_key", "bf_domain_key",
          "bfplus_facet_key", "bfplus_domain_key", "bfplusm_observed", "bridge_items",
          "bfplusm_anankastia_doublets", "bfplusm_anankastia_items", "versions",
          "score_version", "total_items", "transfer_outcomes")
orig <- new.env()
for (e in parse(orig_file, keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, orig)
}
missing <- setdiff(keep, ls(orig))
if (length(missing)) stop("Not found in the original file: ", paste(missing, collapse = ", "))
attach(orig, name = "orig", warn.conflicts = FALSE)

inst <- read_instrument(yml)
reg <- scale_registry(inst)
cat("normpipe outcomes:", nrow(reg), "\n"); print(table(reg$version))
items_of <- function(o) as.numeric(get_scale(inst, o)$items)
res <- list()
chk <- function(name, ok) { res[[name]] <<- ok; if (!ok) cat("MISMATCH:", name, "\n") }
num <- function(x) as.numeric(unlist(x, use.names = FALSE))

rev_np <- inst$items$item[map_lgl(inst$items$values, \(v) v[1] > v[length(v)])]
chk("reverse items", setequal(rev_np, reverse_items))
for (f in names(facet_key))        chk(paste0("full_", f),   identical(items_of(paste0("full_", f)), num(facet_key[[f]])))
for (f in names(sf_facet_key))     chk(paste0("sf_", f),     identical(items_of(paste0("sf_", f)), num(sf_facet_key[[f]])))
for (f in names(bfplus_facet_key)) chk(paste0("bfplus_", f), identical(items_of(paste0("bfplus_", f)), num(bfplus_facet_key[[f]])))
for (d in names(bf_domain_key))    chk(paste0("bf_", d),     identical(items_of(paste0("bf_", d)), num(bf_domain_key[[d]])))
for (d in names(domain_key)) {     # order matters: it defines the bifactor blocks
  chk(paste0("full_", d),   identical(items_of(paste0("full_", d)), num(facet_key[domain_key[[d]]])))
  chk(paste0("sf_", d),     identical(items_of(paste0("sf_", d)), num(sf_facet_key[domain_key[[d]]])))
  chk(paste0("bfplus_", d), identical(items_of(paste0("bfplus_", d)), num(bfplus_facet_key[bfplus_domain_key[[d]]])))
}
chk("bfplus_m_anankastia", identical(items_of("bfplus_m_anankastia"), num(bfplusm_anankastia_items)))
for (v in c("full", "sf", "bf", "bfplus_m"))
  chk(paste0(v, "_total"), identical(items_of(paste0(v, "_total")), num(total_items(v))))

# routes when the reference sample answered only the 34 bridge items (GESIS)
rt <- scale_routes(inst, bridge_items)
to_np <- function(o) o |> sub("^fulldom_", "full_", x = _) |> sub("^sfdom_", "sf_", x = _) |>
  sub("^bfplusm_", "bfplus_m_", x = _)
orig_tr <- to_np(transfer_outcomes())
np_tr <- rt$outcome[rt$route == "transfer"]
chk("transfer outcomes", setequal(np_tr, orig_tr))
cat("\ntransfer: normpipe", length(np_tr), "| original", length(orig_tr), "\n")
cat("direct (normpipe):", paste(rt$outcome[rt$route == "direct"], collapse = ", "), "\n")
cat("original direct: PID5BF+M domains", paste(bfplusm_observed, collapse = ", "), "+ PID5BF+M total\n")

# scoring of random raw responses (reverse-keyed items unreversed)
set.seed(1); N <- 300
raw <- as_tibble(setNames(replicate(220, sample(0:3, N, TRUE), simplify = FALSE), paste0("PID_", 1:220)))
x <- raw; for (i in reverse_items) x[[paste0("PID_", i)]] <- 3 - x[[paste0("PID_", i)]]
np <- score_scales(raw, inst, prorate = FALSE)
diffs <- c()
for (v in c("full", "sf", "bf", "bfplus")) {
  o <- score_version(x, v)
  for (nm in names(o)) diffs[paste0(v, "_", nm)] <- max(abs(np[[paste0(v, "_", nm)]] - o[[nm]]))
}
for (v in c("full", "sf", "bf", "bfplus_m"))
  diffs[paste0(v, "_total")] <- max(abs(np[[paste0(v, "_total")]] - rowMeans(x[paste0("PID_", total_items(v))])))
diffs["bfplus_m_anankastia"] <- max(abs(np$bfplus_m_anankastia - rowMeans(x[paste0("PID_", bfplusm_anankastia_items)])))
chk("scoring", max(diffs) < 1e-12)
cat("\nscoring: ", length(diffs), " outcomes, max |difference| = ", signif(max(diffs), 3), "\n", sep = "")
cat("checks passed:", sum(unlist(res)), "of", length(res), "\n")
if (!all(unlist(res))) quit(status = 1)
