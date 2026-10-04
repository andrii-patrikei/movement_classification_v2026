# som_fits.R: the article's kohonen SOM and the novel_SOM variants behind one interface
#
# A fit is a list: M (codebook, one unit per row), grid_dist and links (the lattice), dist (the map's own
# distance between the rows of two matrices), bmu (the unit of every training item, own distance), seconds.
# The variants are the seventeen of novel_SOM (github.com/andrii-patrikei/novel_SOM, next_steps/, commit 2ae89c5),
# vendored in novel_som/: ns_maps.R (lattices, schedule, variants, train_map) unchanged, and ns_core.cpp (the
# compiled online SOM) with one optional argument added, the bubble neighbourhood of kohonen (off by default).

load_novel_som <- function(dir = "novel_som", cache_dir = file.path(dir, ".cpp_cache")) {
  Sys.setenv(PKG_LIBS = "-lquadmath")
  Rcpp::sourceCpp(file.path(dir, "ns_core.cpp"), cacheDir = cache_dir, env = globalenv())
  sys.source(file.path(dir, "ns_maps.R"), envir = globalenv())
  invisible(TRUE)
}

# the variants whose distance needs a time axis: on a feature vector the "time" is the rank of the feature
# in the XGBoost gain list, so these align neighbouring ranks; run for completeness, not as candidates
ELASTIC <- c("dtw", "softdtw", "sbd", "dtw_anneal")

# a feature matrix as a novel_SOM data set: one item per row, a "series" of ncol points in one channel
feature_ds <- function(X) list(X = X, L = ncol(X), nch = 1L, band_share = 0.1)

# --------------------------------------------------------------------------------- the article's kohonen SOM
# train_som_cv() of the article: aweSOM::somInit(random) and kohonen::som, hexagonal, online, 1000 passes,
# alpha 0.05 to 0.01, sum of squares; set.seed() right before, as the article does
fit_kohonen <- function(X, seed, xdim = 5, ydim = 5, rlen = 1000) {
  set.seed(seed)
  init <- aweSOM::somInit(X, xdim, ydim, "random")
  t0 <- proc.time()[["elapsed"]]
  s <- kohonen::som(X, grid = kohonen::somgrid(xdim, ydim, "hexagonal"), rlen = rlen, alpha = c(0.05, 0.01),
                    dist.fcts = "sumofsquares", init = init, maxNA.fraction = 8)
  gd <- kohonen::unit.distances(s$grid)
  list(method = "kohonen", M = s$codes[[1]], grid_dist = gd, links = round(gd, 3) == 1,
       dist = function(A, B) sqrt(pmax(outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B), 0)),
       euclid_family = TRUE, bmu = s$unit.classif, seconds = proc.time()[["elapsed"]] - t0)
}

# ------------------------------------------------------------------------------------- a novel_SOM variant
# Two training regimes for every variant:
#   "novel"    train_map() of novel_SOM with its own conventions: Gaussian neighbourhood, alpha 0.5 to 0.01,
#              radius from the two-thirds quantile of the lattice distances to 0 (0.5 below 1), 10 passes,
#              the start a draw of K items without replacement
#   "article"  the schedule of the article's kohonen::som call: bubble neighbourhood, alpha 0.05 to 0.01,
#              radius from quantile(unit.distances(grid), 2/3) to 0 (0.5 below 1, so the last part updates the
#              winner only), 1000 passes, the start aweSOM::somInit(random) (K items with replacement)
# In both, the start and the order of items come from the seed, so every variant with one seed sees the same
# start and the same order and can be compared with its control seed by seed.
train_map_article <- function(ds, spec, size, seed, rlen = 1000, alpha = c(0.05, 0.01)) {
  grid <- make_grid(size, size, spec$topo, spec$torus)
  N <- nrow(ds$X); K <- grid$K; S <- rlen * N
  set.seed(seed)
  M0 <- ds$X[sample(N, K, replace = TRUE), , drop = FALSE]
  pick <- as.integer(floor(N * runif(S)) + 1)
  frac <- (seq_len(S) - 1) / S                                    # kohonen: curIter / totalIters
  radius <- unname(quantile(grid$dist, 2/3)) * (1 - frac); radius[radius < 1] <- 0.5
  alphas <- alpha[1] - (alpha[1] - alpha[2]) * frac
  sc <- scales_of(ds); gamma <- spec$gamma * sc$gamma_unit
  temp <- rep(0, S)
  if (spec$temp > 0) {                                            # train_map's cooling: geometric over the first 90%
    s2 <- scale2_under(ds, spec$mode, sc, spec$p, gamma); cool <- floor(0.9 * S)
    temp[1:cool] <- spec$temp * s2 * exp(log(1e-3) * (0:(cool - 1)) / cool)
  }
  t0 <- proc.time()[["elapsed"]]
  r <- som_train_cpp(ds$X, M0, grid$dist, pick, alphas, radius, temp, NS_MODES[[spec$mode]], ds$L, ds$nch,
                     sc$band, spec$p, gamma, sc$maxshift, spec$conscience, spec$huber * sc$huber_unit, spec$momentum,
                     NS_PREC[[spec$prec]], spec$kahan, integer(0), TRUE)
  r$seconds <- proc.time()[["elapsed"]] - t0
  r$grid <- grid; r$scales <- sc
  r
}

fit_novel <- function(X, variant, seed, size = 5, rlen = 10, regime = c("novel", "article")) {
  regime <- match.arg(regime)
  spec <- VARIANTS[[variant]]; ds <- feature_ds(X)
  r <- if (regime == "novel") train_map(ds, spec, size, seed, rlen = rlen) else train_map_article(ds, spec, size, seed)
  sc <- r$scales
  mode <- NS_MODES[[spec$mode]]; p <- spec$p; gamma <- spec$gamma * sc$gamma_unit
  own <- function(A, B) cross_dist_cpp(A, B, mode, ds$L, 1L, sc$band, p, gamma, sc$maxshift)
  fit <- list(method = variant, regime = regime, M = r$M, grid_dist = r$grid$dist, links = r$grid$links == 1, dist = own,
              euclid_family = spec$mode == "euclid", softdtw = spec$mode == "softdtw", seconds = r$seconds)
  fit$bmu <- map_bmu(fit, X)
  fit
}

# the unit of every row of Xnew, in the map's own distance (ties: the first unit)
map_bmu <- function(fit, Xnew) max.col(-fit$dist(Xnew, fit$M), ties.method = "first")

# the 7 super-clusters of the article: PAM on the codebook. On the Euclidean maps literally the article's
# call, cluster::pam(codes, 7); on the others PAM on the own distances between the prototypes (soft-DTW as
# its divergence, sdtw(a, b) - (sdtw(a, a) + sdtw(b, b)) / 2, which is 0 on the diagonal)
super_clusters <- function(fit, k = 7L) {
  if (isTRUE(fit$euclid_family)) return(cluster::pam(fit$M, k)$clustering)
  D <- fit$dist(fit$M, fit$M)
  if (isTRUE(fit$softdtw)) D <- D - outer(diag(D), diag(D), "+") / 2
  D <- pmax((D + t(D)) / 2, 0); diag(D) <- 0
  cluster::pam(stats::as.dist(D), k, diss = TRUE)$clustering
}

# one fitter for both families: method "kohonen", or "<variant>:<regime>" for a novel_SOM variant
# (regime "novel" or "article", see above)
fit_som <- function(X, method, seed, rlen_novel = 10) {
  if (method == "kohonen") return(fit_kohonen(X, seed))
  vr <- strsplit(method, ":", fixed = TRUE)[[1]]
  fit_novel(X, vr[1], seed, rlen = rlen_novel, regime = if (length(vr) > 1) vr[2] else "novel")
}
