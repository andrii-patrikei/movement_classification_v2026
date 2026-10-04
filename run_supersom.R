#!/usr/bin/env Rscript
# Computes the results of supersom.Rmd outside the knit, so the knit takes a minute:
#   Rscript run_supersom.R [cores = all but two] [part = all] [passes = 1000]
# part: "main" (all data and held-out children, every method, with the missing-sensor test inside), "confirm" (the
# screening's winners against the best single-layer maps on fresh seeds; needs "main"), or "all" (both in order).
# Fewer than 1000 passes makes a quick run of 2 seeds into supersom/results_quick/ (for testing the pipeline).
# Run it from the project root. The results go to supersom/results/ under keys of the code, the jobs and the data.
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else max(1L, parallel::detectCores() - 2L)
part <- if (length(args) >= 2) args[2] else "all"
rlen <- if (length(args) >= 3) as.integer(args[3]) else 1000L
quick <- rlen < 1000L; dir <- if (quick) "supersom/results_quick" else "supersom/results"
if (!file.exists("supersom/ss_core.cpp")) stop("run this from the project root (the folder that holds supersom/)")
for (f in c("R/article_metrics.R", "R/protocols.R", "supersom/R/ss_maps.R", "supersom/R/ss_protocols.R", "supersom/R/ss_bench.R",
            "supersom/R/ss_raw.R")) source(f)
load_ss_core()
if (!file.exists("supersom/cache/article_27.rds") || !file.exists("supersom/cache/raw_layers.rds")) prepare_ss_cache()
S <- ss_data(rlen = rlen)
plan <- if (quick) ss_plan(1:2, 1:2) else ss_plan()
t0 <- Sys.time()
if (part %in% c("main", "all")) {
  cat(sprintf("%d cores, %d passes; main: %s jobs\n", cores, rlen, paste(sapply(plan, nrow), collapse = " + ")))
  res <- ss_compute(S, plan, cores, dir)
  for (p in names(res)) if (!is.null(attr(res[[p]], "failed"))) { cat(p, ": failed jobs\n"); print(attr(res[[p]], "failed")) }
}
if (part %in% c("confirm", "all")) {
  res <- ss_compute(S, plan, cores, dir)                   # from the cache
  print(unlist(pick_winners(res$heldout)))
  res_cf <- ss_compute(S, ss_plan_confirm(res, if (quick) 101:102 else 101:120), cores, dir)
}
cat("done in", format(round(Sys.time() - t0, 1)), "\n")
