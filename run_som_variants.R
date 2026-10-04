#!/usr/bin/env Rscript
# Computes the results of som_variants.Rmd outside the knit (about 1.5 hours on 18 cores), so the knit takes a minute:
#   Rscript run_som_variants.R [cores = all but two]
# Run it from the project root. The results go to results/ under keys made of the code, the jobs, the settings and
# the article's file; the knit finds them there. The article's saved models are downloaded to data/ if missing.
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else max(1L, parallel::detectCores() - 2L)
if (!file.exists("novel_som/ns_core.cpp")) stop("run this from the project root (the folder that holds novel_som/)")
for (f in c("R/article_metrics.R", "R/som_fits.R", "R/protocols.R", "R/bench.R")) source(f)
load_novel_som()
D <- article_data(get_article_rds())
plan <- bench_plan()
cat(sprintf("%d cores; jobs: %s; estimated %.1f h of CPU\n", cores,
            paste(sprintf("%s %d", names(plan), sapply(plan, nrow)), collapse = ", "), sum(sapply(plan, function(j) sum(j$cost))) / 3600))
t0 <- Sys.time()
res <- compute_all(D, plan, cores)
for (p in names(plan)) {
  failed <- attr(res[[p]], "failed")
  if (!is.null(failed)) { cat(p, ": failed jobs\n"); print(failed) }
}
cat(sprintf("%s rows in %s\n", paste(sprintf("%s %d", names(res), sapply(res, nrow)), collapse = ", "), format(round(Sys.time() - t0, 1))))
