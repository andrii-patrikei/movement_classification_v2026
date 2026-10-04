# checks of the SuperSOM core (run from the project root: Rscript supersom/tests/test_core.R)
#  1. one sum-of-squares layer, the article's schedule: kohonen::som, bit for bit (the article's call, aweSOM start)
#  2. four sum-of-squares layers: kohonen::supersom, bit for bit, including its layer weights
#  3. the DTW of a curve layer: the square of novel_SOM's DTW distance (ns_core.cpp, dependent multichannel DTW)
#  4. the weighted sum of the per-layer distances is the combined distance
#  5. the item warped onto a prototype by ss_warp_cpp is the warp of novel_SOM's item_seen_cpp
for (f in c("R/article_metrics.R", "R/som_fits.R", "supersom/R/ss_maps.R", "supersom/R/ss_protocols.R")) source(f)
load_ss_core()
A <- readRDS("supersom/cache/article_27.rds"); RL <- readRDS("supersom/cache/raw_layers.rds")
XL <- lapply(A$X, feature_layer)
ok <- function(what, cond) { cat(sprintf("%-78s %s\n", what, if (isTRUE(cond)) "ok" else "FAILED")); if (!isTRUE(cond)) quit(status = 1) }

d1 <- sapply(names(XL), function(s) sapply(1:3, function(seed) max(abs(fit_kohonen(A$X[[s]], seed)$M - fit_ss(XL, ss_method("k", "", s, "classic"), seed)$M))))
ok("1. one layer = kohonen::som (5 feature sets x 3 seeds, largest difference 0)", all(d1 == 0))
ok("   the article's published tsfresh map (seed 1245) comes back exactly",
   max(abs(fit_ss(XL, ss_method("k", "", "tsfresh", "classic"), 1245)$M - A$pub_codes$tsfresh)) == 0)

F4 <- c("rqa", "autocorr", "spectr", "tsfresh"); grid <- kohonen::somgrid(5, 5, "hexagonal")
d2 <- sapply(1:3, function(seed) {
  set.seed(seed); idx <- sample(648, 25, replace = TRUE)
  ks <- kohonen::supersom(A$X[F4], grid = grid, rlen = 1000, alpha = c(0.05, 0.01), init = lapply(A$X[F4], function(x) x[idx, ]),
                          dist.fcts = "sumofsquares", maxNA.fraction = 8)
  f <- fit_ss(XL, ss_method("s", "", F4, "classic"), seed)
  w <- ks$user.weights * ks$distance.weights
  c(codes = max(abs(do.call(cbind, ks$codes) - f$M)), weights = max(abs(w / sum(w) - f$w)))
})
ok("2. four layers = kohonen::supersom (3 seeds, codebooks and weights identical)", all(d2 == 0))

Sys.setenv(PKG_LIBS = "-lquadmath")
Rcpp::sourceCpp("novel_som/ns_core.cpp", cacheDir = file.path(tempdir(), "nscpp_test"), env = (ns <- new.env()))
P <- standardize_layers(RL["posture"])[[1]]
sp <- layer_spec(list(posture = P), "dtw")
mine <- ss_cross(P[1:40, ], P[41:60, ], sp, 1)
theirs <- ns$cross_dist_cpp(unclass(P)[1:40, ], unclass(P)[41:60, ], 6L, 64L, 3L, sp$band, 2, 0.1, sp$band)^2
ok("3. DTW cost of a 3-channel curve layer = novel_SOM's DTW squared (relative 1e-12)", max(abs(mine - theirs) / theirs) < 1e-12)

Xc <- cbind(unclass(P), A$X$spectr); sp2 <- layer_spec(list(posture = P, spectr = XL$spectr), c("dtw", "manhattan")); w <- c(0.3, 0.7)
Dl <- ss_cross(Xc[1:30, ], Xc[31:40, ], sp2, w, per_layer = TRUE)
ok("4. weighted sum of the per-layer distances = the combined distance", max(abs(w[1] * Dl[[1]] + w[2] * Dl[[2]] - ss_cross(Xc[1:30, ], Xc[31:40, ], sp2, w))) < 1e-9)

x <- as.numeric(P[3, ]); m <- as.numeric(P[200, ])
ok("5. the warped item = novel_SOM's item_seen_cpp", max(abs(ss_warp_cpp(x, m, 64L, 3L, sp$band) - ns$item_seen_cpp(x, m, 6L, 64L, 3L, sp$band, 2, 0.1, sp$band))) < 1e-12)
cat("all checks passed\n")
