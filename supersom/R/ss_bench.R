# ss_bench.R: the data, the methods, the jobs, the parallel runs and the cache of supersom.Rmd
#
# The methods were fixed before the full runs (after a pilot of a few seeds on all data and one held-out seed,
# described in the report); the winner among the layered maps is chosen by a rule fixed in advance and then
# tested on fresh seeds against the best single-layer map.

SS_CODE <- c("supersom/ss_core.cpp", "supersom/R/ss_maps.R", "supersom/R/ss_protocols.R", "supersom/R/ss_bench.R",
             "R/article_metrics.R", "R/protocols.R")
F4 <- c("rqa", "autocorr", "spectr", "tsfresh")             # the article's feature sets (Mix is the first three together)
R3 <- c("posture", "accspec", "gyrospec")                  # the raw-signal curves (ss_raw.R)
POOL7 <- c(F4, R3)
LAYER_LABELS <- c(rqa = "RQA", autocorr = "Autocorrelation", spectr = "Spectral", mix = "Mix", tsfresh = "tsfresh",
                  posture = "Posture curve", accspec = "Acc. spectrum", gyrospec = "Gyro spectrum")

# every layer of the study, in the article's order of segments: the five feature matrices (27 features, the
# article's top XGBoost features, scaled) and the three raw curves; plus the 9 children without a gyroscope
ss_data <- function(article27 = "supersom/cache/article_27.rds", raw = "supersom/cache/raw_layers.rds", rlen = 1000L) {
  A <- readRDS(article27); RL <- readRDS(raw)
  stopifnot(identical(RL$segment, paste(A$child, rep(1:8, 81))))
  XL <- c(lapply(A$X, feature_layer), standardize_layers(RL[R3]))
  list(XL = XL, y7 = A$y7, child = A$child, published = A$published, seconds = RL$seconds, md5 = A$md5,
       no_gyro = RL$no_gyro, raw_scale = RL[R3], rlen = rlen)
}

M <- function(...) ss_method(...)
SS_METHODS <- list(
  # one layer: the controls
  M("kohonen_tsfresh", "kohonen::som on tsfresh (the article)", "tsfresh", "classic", group = "single"),
  M("rob_rqa", "Robust SOM, RQA", "rqa", "robust", group = "single"),
  M("rob_autocorr", "Robust SOM, autocorrelation", "autocorr", "robust", group = "single"),
  M("rob_spectr", "Robust SOM, spectral", "spectr", "robust", group = "single"),
  M("rob_mix", "Robust SOM, Mix", "mix", "robust", group = "single"),
  M("rob_tsfresh", "Robust SOM, tsfresh", "tsfresh", "robust", group = "single"),
  M("dtw_posture", "DTW-SOM, posture curve", "posture", "robust", group = "single"),
  M("euc_posture", "Euclidean SOM, posture curve", "posture", "classic", group = "single"),
  M("dtw_accspec", "DTW-SOM, acc. spectrum", "accspec", "robust", group = "single"),
  M("dtw_gyrospec", "DTW-SOM, gyro spectrum", "gyrospec", "robust", group = "single"),
  # layered: the article's four feature sets
  M("ss_classic_F4", "kohonen::supersom, 4 feature sets", F4, "classic", group = "features"),
  M("glue_classic_F4", "kohonen::som on the 4 sets glued (108 features)", F4, "classic", glue = TRUE, group = "features"),
  M("ss_rob_F4", "Robust SuperSOM, 4 feature sets", F4, "robust", group = "features"),
  M("glue_rob_F4", "Robust SOM on the 4 sets glued", F4, "robust", glue = TRUE, group = "features"),
  M("ss_rob_F4_agree", "Robust SuperSOM, 4 sets, agreement weights", F4, "robust", weights = "agree", group = "features"),
  M("ss_rob_F4_adapt", "Robust SuperSOM, 4 sets, inverse-error weights", F4, "robust", weights = "adapt", group = "features"),
  M("gss_rob_F4", "Growing robust SuperSOM, 4 sets", F4, "robust", grow = TRUE, float = 1, group = "features"),
  M("xyf_rob_F4", "Supervised robust SuperSOM (task layer), 4 sets", F4, "robust", label_weight = 0.5, group = "supervised"),
  # layered: the raw curves, with DTW
  M("ss_R3_dtw", "DTW SuperSOM, 3 raw curves", R3, "robust", group = "raw"),
  M("ss_R3_euc", "Euclidean SuperSOM, 3 raw curves", R3, "classic", group = "raw"),
  M("gss_R3_dtw", "Growing DTW SuperSOM, 3 raw curves", R3, "robust", grow = TRUE, float = 1, group = "raw"),
  M("ss_A2_dtw", "DTW SuperSOM, accelerometer curves only", c("posture", "accspec"), "robust", group = "raw"),
  # layered: all seven layers, and the layers each fold chooses for itself
  M("ss_all7", "Robust DTW SuperSOM, all 7 layers", POOL7, "robust", group = "all"),
  M("ss_all7_agree", "Robust DTW SuperSOM, all 7 layers, agreement weights", POOL7, "robust", weights = "agree", group = "all"),
  M("gss_all7", "Growing robust DTW SuperSOM, all 7 layers", POOL7, "robust", grow = TRUE, float = 1, group = "all"),
  M("ss_sel7", "Robust DTW SuperSOM, layers chosen on the training children", POOL7, "robust", select = TRUE, group = "selected"),
  M("ss_sel7_euc", "Robust SuperSOM, layers chosen, curves Euclidean", POOL7, "robust_euclid", select = TRUE, group = "selected"),
  M("gss_sel7", "Growing robust DTW SuperSOM, layers chosen", POOL7, "robust", select = TRUE, grow = TRUE, float = 1, group = "selected")
)
names(SS_METHODS) <- vapply(SS_METHODS, `[[`, "", "id")
UNSUPERVISED_LAYERED <- names(SS_METHODS)[vapply(SS_METHODS, function(m) m$group %in% c("features", "raw", "all", "selected"), TRUE)]

# seconds per map on all 648 segments (measured), to start the longest jobs first
ss_cost <- function(id) {
  m <- SS_METHODS[[id]]; layers <- m$layers
  dtw <- sum(layers %in% R3) * (m$modes %in% c("robust"))
  base <- 0.5 + 1.2 * sum(layers %in% F4 | layers == "mix") + 25 * dtw
  if (m$weights %in% c("agree", "adapt")) base <- base * 2
  if (m$grow) base <- base / 3
  base
}

make_ss_jobs <- function(methods, seeds, protocol) {
  jobs <- expand.grid(method = methods, seed = seeds, stringsAsFactors = FALSE)
  jobs$protocol <- protocol
  jobs$cost <- vapply(jobs$method, ss_cost, 0) * c(all = 1, heldout = 9, confirm = 9)[[protocol]]
  jobs[order(-jobs$cost), ]
}

# the plan: all data with 10 seeds, held-out children with 10 seeds (every prediction kept); the confirmation on
# fresh seeds and the missing-sensor test are separate parts (the confirmation needs the screening's winner)
ss_plan <- function(seeds_all = 1:10, seeds_cv = 1:10) {
  ids <- names(SS_METHODS)
  list(all = make_ss_jobs(ids, seeds_all, "all"), heldout = make_ss_jobs(ids, seeds_cv, "heldout"))
}

run_ss_job <- function(job, S, rlen = S$rlen) {
  meth <- SS_METHODS[[job$method]]
  if (job$protocol == "all") return(cbind(job[, c("method", "seed")], eval_all_data_ss(S$XL, S$y7, meth, job$seed, rlen = rlen, child = S$child,
                                                                                     no_gyro = S$no_gyro, raw_scale = S$raw_scale)))
  if (job$protocol %in% c("heldout", "confirm")) {
    p <- heldout_predictions_ss(S$XL, S$y7, S$child, meth, job$seed, rlen = rlen)
    return(cbind(method = job$method, seed = job$seed, p, layers = attr(p, "layers"), w = paste(sprintf("%.3f", attr(p, "w")), collapse = " ")))
  }
  stop("unknown protocol ", job$protocol)
}

# a socket cluster (Windows has no fork); every worker compiles ss_core.cpp into its own cache, one at a time
run_ss_jobs <- function(jobs, S, cores, root = getwd()) {
  if (cores <= 1) return(data.table::rbindlist(lapply(seq_len(nrow(jobs)), function(i) run_ss_job(jobs[i, ], S)), fill = TRUE))
  cl <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  for (i in seq_along(cl)) parallel::clusterCall(cl[i], function(root) {
    setwd(root)
    for (f in c("R/article_metrics.R", "R/protocols.R", "supersom/R/ss_maps.R", "supersom/R/ss_protocols.R", "supersom/R/ss_bench.R")) source(f)
    load_ss_core(cache_dir = file.path(tempdir(), "sscpp"))
    NULL
  }, root)
  parallel::clusterExport(cl, c("S", "jobs"), envir = environment())
  res <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), function(i) tryCatch(run_ss_job(jobs[i, ], S), error = function(e) conditionMessage(e)))
  ok <- !vapply(res, is.character, logical(1))
  out <- data.table::rbindlist(res[ok], fill = TRUE)
  if (any(!ok)) {
    warning(sum(!ok), " of ", nrow(jobs), " jobs failed, first: ", res[[which(!ok)[1]]])
    attr(out, "failed") <- data.frame(jobs[!ok, c("method", "seed", "protocol")], error = unlist(res[!ok]))
  }
  out
}

ss_key <- function(...) {
  code <- lapply(SS_CODE, function(f) if (file.exists(f)) readLines(f, warn = FALSE) else NA)
  substr(digest::digest(list(code, list(...)), algo = "sha1"), 1, 12)
}
ss_cached <- function(prefix, key, fn, dir = "supersom/results") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, sprintf("%s_%s.rds", prefix, key))
  if (file.exists(path)) return(readRDS(path))
  t0 <- Sys.time(); value <- fn()
  attr(value, "seconds") <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(attr(value, "failed"))) saveRDS(value, path)
  value
}
ss_compute <- function(S, plan, cores, dir = "supersom/results") {
  lapply(setNames(names(plan), names(plan)), function(part) {
    jobs <- plan[[part]]
    ss_cached(paste0("ss_", part), ss_key(jobs[, c("method", "seed", "protocol")], S$md5, S$rlen), function() run_ss_jobs(jobs, S, cores), dir)
  })
}

# the confirmation: the rule fixed in advance picks, from the screening (held-out seeds 1-10), the best unsupervised
# layered map and the best single-layer map, each by the mean held-out accuracy under the article's super-cluster
# rule and under the unit rule; these and the article's map run again on fresh seeds 101-120 (new starts, new
# orders, new folds of children), every prediction kept for the child-by-child test
pick_winners <- function(H) {
  acc <- H[, .(super = mean(pred_super == y7), unit = mean(pred_unit == y7)), by = .(method, seed)][, .(super = mean(super), unit = mean(unit)), by = method]
  lay <- acc[method %in% UNSUPERVISED_LAYERED]; one <- acc[method %in% names(SS_METHODS)[vapply(SS_METHODS, `[[`, "", "group") == "single"]]
  list(layered_super = lay[which.max(super)]$method, layered_unit = lay[which.max(unit)]$method,
       single_super = one[which.max(super)]$method, single_unit = one[which.max(unit)]$method)
}
ss_plan_confirm <- function(res, seeds = 101:120) {
  w <- pick_winners(res$heldout)
  jobs <- make_ss_jobs(unique(c(unlist(w), "kohonen_tsfresh")), seeds, "confirm")
  attr(jobs, "picked") <- w
  list(confirm = jobs)
}

# the two caches the study starts from: the article's 27-feature matrices (read from its saved models) and the raw
# curves (from the Zenodo file); both are rebuilt only when missing
prepare_ss_cache <- function(article_rds = "data/som_gain_2_35f_5x5.RDS", csv = "data/physical_exercise.csv") {
  dir.create("supersom/cache", showWarnings = FALSE, recursive = TRUE)
  if (!file.exists("supersom/cache/article_27.rds")) {
    source("R/bench.R", local = (e <- new.env())); D <- e$article_data(e$get_article_rds(article_rds))
    A <- list(X = lapply(D$X, function(s) s[["27"]]), y7 = D$y7, y8 = D$y8, child = D$child, md5 = D$md5, published = D$published,
              pub_codes = lapply(setNames(names(D$X), names(D$X)), function(s) D$models[[27]]$som[[s]]$SOM_results$som_model$codes[[1]]))
    saveRDS(A, "supersom/cache/article_27.rds")
  }
  if (!file.exists("supersom/cache/raw_layers.rds")) {
    A <- readRDS("supersom/cache/article_27.rds")
    saveRDS(build_raw_layers(csv, A$child), "supersom/cache/raw_layers.rds")
  }
  invisible(TRUE)
}
