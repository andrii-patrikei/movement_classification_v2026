# bench.R: the article's data, the jobs, the parallel runs and the cache of som_variants.Rmd

BENCH_CODE <- c("R/article_metrics.R", "R/som_fits.R", "R/protocols.R", "R/bench.R", "novel_som/ns_core.cpp", "novel_som/ns_maps.R")
# the article's five feature sets; Mix is RQA, autocorrelation and spectral features together, tsfresh all of tsfresh
FEATURE_SETS <- c(rqa = "RQA", autocorr = "Autocorrelation", spectr = "Spectral", mix = "Mix", tsfresh = "tsfresh (all)")

# the article's saved maps and feature matrices (trained_models/som_gain_2_35f_5x5.RDS of its repository,
# https://code.it4i.cz/ADAS/movement-classification), downloaded once when missing and checked by md5
ARTICLE_RDS_URL <- paste0("https://code.it4i.cz/api/v4/projects/ADAS%2Fmovement-classification/repository/files/",
                          "trained_models%2Fsom_gain_2_35f_5x5.RDS/raw?ref=main")
ARTICLE_RDS_MD5 <- "e7d4dcda41df3425163c1c8ae4bbe5e9"
get_article_rds <- function(path = "data/som_gain_2_35f_5x5.RDS") {
  if (!file.exists(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    message("Downloading the article's saved models (95 MB) from code.it4i.cz ...")
    old <- options(timeout = 3600); on.exit(options(old))
    utils::download.file(ARTICLE_RDS_URL, path, mode = "wb", quiet = TRUE)
  }
  if (!identical(unname(tools::md5sum(path)), ARTICLE_RDS_MD5)) warning(path, " is not the article's published file (md5 differs)")
  path
}

# the article's feature matrices, exactly as its SOMs saw them: for n features and every feature set, the n
# features with the highest XGBoost gain (selected by the article on all 648 segments), scaled, as
# train_som_cv() prepares them; plus the labels and the child of every segment
article_data <- function(rds, n_range = 2:35) {
  m <- readRDS(rds)
  first <- m[[max(n_range)]]$train_data$tsfresh$df_selected_features_xgboost
  X <- lapply(names(FEATURE_SETS), function(set) lapply(setNames(n_range, n_range), function(n) {
    d <- m[[n]]$train_data[[set]]$df_selected_features_xgboost
    Z <- as.matrix(d[, setdiff(names(d), c("student", "activity")), drop = FALSE])
    scale(Z[, colSums(is.na(Z)) == 0, drop = FALSE])
  }))
  names(X) <- names(FEATURE_SETS)
  # the article's own numbers for every saved map (seed 1245), for the comparison
  published <- do.call(rbind, lapply(n_range, function(n) {
    q <- as.data.frame(m[[n]]$quality[[1]])
    data.frame(set = names(FEATURE_SETS), n = n, QE = q[["Quantization error"]], EV = q[["(% explained variance)"]],
               TE = q[["Topographic error"]], KL = q[["Kaski-Lagus error"]], accuracy = q[["Accuracy rate all"]],
               mnLogLoss = unname(q[["mnLogLoss all"]]), kfold_train = q[["Accur mean k-fold tr"]],
               kfold_test = q[["Accur mean k-fold test"]], stringsAsFactors = FALSE)
  }))
  y8 <- as.integer(factor(first$activity, levels = unique(first$activity)))   # the order of the protocol
  stopifnot(identical(y8, rep(1:8, 81)))
  list(X = X, y8 = y8, y7 = task7_of(y8), child = as.character(first$student), published = published,
       importance = m[[27]]$train_data$tsfresh$df_important_feat_xgboost, models = m, md5 = unname(tools::md5sum(rds)))
}

# seconds per map on the 648 x 27 matrices (measured), to start the longest jobs first
fit_cost <- function(method) {
  v <- sub(":.*$", "", method); art <- grepl(":article$", method)
  ifelse(method == "kohonen", 1.1,
  ifelse(art, c(minkowski3 = 37, dtw = 20, dtw_anneal = 27, sbd = 7.4, anneal = 5.3)[v], c(softdtw = 45, dtw = 1.2, dtw_anneal = 1.2, minkowski3 = 1.1)[v]))
}

# one row per job: feature set, number of features, method, seed, and which protocols to run ("all", "kfold",
# "heldout", or "preds": the held-out children with the prediction of every segment kept)
make_jobs <- function(sets, n_range, methods, seeds, protocols = c("all", "kfold", "heldout")) {
  jobs <- expand.grid(set = sets, n = n_range, method = methods, seed = seeds, stringsAsFactors = FALSE)
  jobs$protocols <- paste(protocols, collapse = ",")
  cost <- fit_cost(jobs$method); cost[is.na(cost)] <- ifelse(grepl(":article$", jobs$method[is.na(cost)]), 2, 0.65)
  jobs$cost <- cost * (("all" %in% protocols) + 8 * ("kfold" %in% protocols) + 9 * any(c("heldout", "preds") %in% protocols))
  jobs[order(-jobs$cost), ]
}

# The plan: the methods are kohonen::som (the article) and every novel_SOM variant under both regimes, except
# soft-DTW under the article's 1000 passes (about an hour per map).
#   main     27 features, every feature set: all data with 20 seeds for every method; the article's k-fold and
#            the held-out children with the first 10 seeds, for the methods whose distance means something on a
#            feature vector (the elastic ones, and Minkowski under the article's schedule at 37 s a map, are
#            judged on all data only)
#   curve    2 to 35 features (the article's range; its Table 5 is 24 to 27): all data, 10 seeds, kohonen, the
#            variants meaningful on features under the article's schedule (without Minkowski), and the
#            Euclidean control under novel_SOM's own schedule
#   confirm  the variant the screening picks (Section 3 of the report: the best held-out accuracy over the five
#            feature sets, which is the Manhattan winner with the clipped pull) against kohonen, on 20 fresh seeds
#            (101 to 120, so new starts, new orders and new folds of children), with every prediction kept
bench_plan <- function(n_main = 27, n_curve = 2:35, seeds_main = 1:20, seeds_cv = 1:10, seeds_curve = 1:10,
                       seeds_confirm = 101:120, confirmed = "l1_huber:article") {
  sets <- names(FEATURE_SETS); v <- names(VARIANTS); flat <- setdiff(v, ELASTIC)
  m_all <- c("kohonen", paste0(setdiff(v, "softdtw"), ":article"), paste0(v, ":novel"))
  m_cv <- c("kohonen", paste0(setdiff(flat, "minkowski3"), ":article"), paste0(flat, ":novel"))
  m_curve <- c("kohonen", paste0(setdiff(flat, "minkowski3"), ":article"), "euclid:novel")
  list(main_all = make_jobs(sets, n_main, m_all, seeds_main, "all"),
       main_cv = make_jobs(sets, n_main, m_cv, seeds_cv, c("kfold", "heldout")),
       curve = make_jobs(sets, n_curve, m_curve, seeds_curve, "all"),
       confirm = make_jobs(sets, n_main, c("kohonen", confirmed), seeds_confirm, "preds"))
}

run_job <- function(job, D, rlen_novel) {
  X <- D$X[[job$set]][[as.character(job$n)]]
  pr <- strsplit(job$protocols, ",")[[1]]
  long <- job$method == "kohonen" || grepl(":article$", job$method)
  out <- cbind(job[, c("set", "n", "method", "seed")], rlen = if (long) 1000L else rlen_novel)
  if ("preds" %in% pr) return(cbind(job[, c("set", "n", "method", "seed")], heldout_predictions(X, D$y7, D$child, job$method, job$seed,
                                                                                            rlen_novel = rlen_novel), row.names = NULL))
  if ("all" %in% pr) out <- cbind(out, eval_all_data(X, D$y7, job$method, job$seed, rlen_novel))
  if ("kfold" %in% pr) out <- cbind(out, eval_article_kfold(X, D$y8, D$y7, job$method, job$seed, rlen_novel = rlen_novel))
  if ("heldout" %in% pr) out <- cbind(out, eval_heldout(X, D$y7, D$child, job$method, job$seed, rlen_novel = rlen_novel))
  out
}

# every part of the plan, each cached under a key of the code, the jobs, the settings and the article's file
compute_all <- function(D, plan, cores, rlen_novel = 10, rds_md5 = D$md5, dir = "results") {
  lapply(setNames(names(plan), names(plan)), function(part) {
    jobs <- plan[[part]]
    cached(paste0("somvar_", part), bench_key(jobs[, c("set", "n", "method", "seed", "protocols")], rlen_novel, rds_md5),
           function() run_jobs(jobs, D, rlen_novel, cores), dir)
  })
}

# a socket cluster (Windows has no fork); every worker compiles ns_core.cpp into its own cache, because
# sourceCpp() writes into its cache directory even when the build is there and the workers would collide
run_jobs <- function(jobs, D, rlen_novel, cores, root = getwd()) {
  if (cores <= 1) return(data.table::rbindlist(lapply(seq_len(nrow(jobs)), function(i) run_job(jobs[i, ], D, rlen_novel)), fill = TRUE))
  cl <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(root) {
    setwd(root)
    for (f in c("R/article_metrics.R", "R/som_fits.R", "R/protocols.R", "R/bench.R")) source(f)
    load_novel_som(cache_dir = file.path(tempdir(), "nscpp"))
    NULL
  }, root)
  D$models <- NULL                                                   # the workers need the matrices only
  parallel::clusterExport(cl, c("D", "jobs", "rlen_novel"), envir = environment())
  res <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), function(i)
    tryCatch(run_job(jobs[i, ], D, rlen_novel), error = function(e) conditionMessage(e)))
  ok <- !vapply(res, is.character, logical(1))                       # a failed job comes back as its error message
  out <- data.table::rbindlist(res[ok], fill = TRUE)
  if (any(!ok)) {
    warning(sum(!ok), " of ", nrow(jobs), " jobs failed, first: ", res[[which(!ok)[1]]])
    attr(out, "failed") <- data.frame(jobs[!ok, c("set", "n", "method", "seed")], error = unlist(res[!ok]))
  }
  out
}

bench_key <- function(...) {
  code <- lapply(BENCH_CODE, function(f) if (file.exists(f)) readLines(f, warn = FALSE) else NA)
  substr(digest::digest(list(code, list(...)), algo = "sha1"), 1, 12)
}
cached <- function(prefix, key, fn, dir = "results") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, sprintf("cache_%s_%s.rds", prefix, key))
  if (file.exists(path)) return(readRDS(path))
  t0 <- Sys.time(); value <- fn()
  attr(value, "seconds") <- as.numeric(Sys.time() - t0, units = "secs")
  saveRDS(value, path)
  value
}
