# ss_maps.R: layers, layer weights, the fixed and the growing SuperSOM, and the maps of supersom.Rmd
#
# A data set is a list of layers, each a matrix with one row per item and the attributes
#   kind   "features" (a feature vector: standardised column by column) or "curves" (series of L points in nch
#          channels, stored channel by channel: standardised as a whole, so the shape of every curve survives)
#   L, nch the length and the number of channels of a curve layer (a feature layer is L = ncol, nch = 1)
# A map is fitted on the layers side by side (the columns of all layers, in order); every layer keeps its own
# distance, weight and pull (ss_core.cpp).

SS_MODES <- c(sumsq = 0L, manhattan = 1L, euclid = 2L, correlation = 3L, dtw = 4L)

load_ss_core <- function(dir = "supersom", cache_dir = file.path(dir, ".cpp_cache")) {
  Rcpp::sourceCpp(file.path(dir, "ss_core.cpp"), cacheDir = cache_dir, env = globalenv())
  invisible(TRUE)
}

feature_layer <- function(X) structure(as.matrix(X), kind = "features", L = ncol(X), nch = 1L)
curve_layer <- function(X, L, nch) { stopifnot(ncol(X) == L * nch); structure(as.matrix(X), kind = "curves", L = as.integer(L), nch = as.integer(nch)) }
keep_attr <- function(new, old) { for (a in c("kind", "L", "nch")) attr(new, a) <- attr(old, a); new }

# standardise every layer on the rows 'tr' and apply the same scaling to all rows (held-out children are
# scaled with the training children's means and sds); a feature layer column by column (the article's
# scale()), a curve layer channel by channel with one mean and sd over all its points
standardize_layers <- function(layers, tr = seq_len(nrow(layers[[1]]))) {
  lapply(layers, function(X) {
    if (attr(X, "kind") == "features") {
      mu <- colMeans(X[tr, , drop = FALSE]); sdv <- apply(X[tr, , drop = FALSE], 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
      keep_attr(scale(X, mu, sdv), X)
    } else {
      L <- attr(X, "L"); Z <- X
      for (ch in seq_len(attr(X, "nch"))) {
        cols <- (ch - 1) * L + seq_len(L); v <- X[tr, cols]
        Z[, cols] <- (X[, cols] - mean(v)) / sd(as.vector(v))
      }
      keep_attr(Z, X)
    }
  })
}

# the layout of the layers side by side, for ss_train_cpp() and ss_cross_cpp()
layer_spec <- function(layers, modes, clip = 0, band_share = 0.1) {
  len <- vapply(layers, ncol, 1L); nl <- length(layers)
  modes <- rep_len(modes, nl); clip <- rep_len(clip, nl)
  L <- vapply(layers, function(X) as.integer(attr(X, "L")), 1L); nch <- vapply(layers, function(X) as.integer(attr(X, "nch")), 1L)
  data.frame(layer = names(layers), off = as.integer(c(0, cumsum(len)[-nl])), len = as.integer(len), mode = unname(SS_MODES[modes]),
             mode_name = modes, L = L, nch = nch, band = as.integer(pmax(1, round(band_share * L))), clip = as.numeric(clip),
             stringsAsFactors = FALSE)
}

ss_cross <- function(X, M, spec, w, per_layer = FALSE)
  ss_cross_cpp(X, M, spec$off, spec$len, spec$mode, spec$L, spec$nch, spec$band, as.numeric(w), per_layer)

# kohonen's normalizeDataLayers: every layer weighted by 1 / (the mean distance between the start rows in that
# layer's own distance), times the user weights, normalised to sum 1 (supersom() in kohonen 3.0)
kohonen_weights <- function(X0, spec, user = rep(1, nrow(spec))) {
  user <- user / sum(user)
  scale <- vapply(seq_len(nrow(spec)), function(l) {
    s1 <- spec[l, ]; s1$off <- 0L
    D <- ss_cross(X0[, spec$off[l] + seq_len(spec$len[l]), drop = FALSE], X0[, spec$off[l] + seq_len(spec$len[l]), drop = FALSE], s1, 1)
    mean(D[lower.tri(D)])
  }, 0)
  w <- user / scale
  w / sum(w)
}

# --------------------------------------------------------------------------------------------- lattices
hex_grid <- function(xdim, ydim) {                      # kohonen's hexagonal somgrid and unit.distances
  g <- kohonen::somgrid(xdim, ydim, "hexagonal")
  list(pos = g$pts, dist = kohonen::unit.distances(g), K = nrow(g$pts))
}
# the growing lattice (from my earlier growing DTW-SOM): sites (column, row) of a hexagonal lattice, odd rows
# shifted half a unit; the map distance is the number of steps along the map
site_pos <- function(sites) cbind(x = sites[, 1] + 0.5 * (sites[, 2] %% 2), y = sites[, 2] * sqrt(3) / 2)
touching <- function(site, sites) abs(sqrt(colSums((t(site_pos(sites)) - c(site_pos(rbind(site))))^2)) - 1) < 1e-6
free_sites <- function(sites) {
  cand <- as.matrix(expand.grid(x = (min(sites[, 1]) - 1):(max(sites[, 1]) + 1), y = (min(sites[, 2]) - 1):(max(sites[, 2]) + 1)))
  cand <- cand[!(paste(cand[, 1], cand[, 2]) %in% paste(sites[, 1], sites[, 2])), , drop = FALSE]
  cand[apply(cand, 1, function(st) any(touching(st, sites))), , drop = FALSE]
}
hop_dist <- function(sites) {
  D <- as.matrix(dist(site_pos(sites))); D <- ifelse(abs(D - 1) < 1e-6, 1, Inf); diag(D) <- 0
  for (k in seq_len(nrow(D))) D <- pmin(D, outer(D[, k], D[k, ], "+"))
  unname(D)
}
start_sites <- function(n = 7) rbind(c(0, 0), free_sites(rbind(c(0, 0))))[seq_len(n), , drop = FALSE]

# --------------------------------------------------------------------------------------------- methods
# A method: which layers, the distance and the pull of every layer, how the layers are weighted, and whether
# the lattice grows. 'family' names the control a method must beat.
ss_method <- function(id, label, layers, modes = "sumsq", clip = 0, weights = "kohonen", grow = FALSE, float = 0,
                      n_start = 7, grow_share = 0.6, grow_radius = 1, birth = "medoid", label_weight = 0, glue = FALSE,
                      select = FALSE, alpha = c(0.05, 0.01), bubble = TRUE, group = "", what = "")
  list(id = id, label = label, layers = layers, modes = modes, clip = clip, weights = weights, grow = grow, float = float,
       n_start = n_start, grow_share = grow_share, grow_radius = grow_radius, birth = birth, label_weight = label_weight,
       glue = glue, select = select, alpha = alpha, bubble = bubble, group = group, what = what)

# the distance of every layer from its kind: "robust" = the Manhattan winner with the clipped pull on a feature layer
# (the screening winner of som_variants.Rmd) and DTW on a curve layer; "robust_euclid" = the same with the sum of
# squares on curves (the control for DTW); "classic" = kohonen's sum of squares on every layer
resolve_modes <- function(XL, modes, clip) {
  kinds <- vapply(XL, function(X) attr(X, "kind"), "")
  if (length(modes) == 1 && modes %in% c("robust", "robust_euclid", "classic")) {
    m <- switch(modes, robust = ifelse(kinds == "features", "manhattan", "dtw"),
                robust_euclid = ifelse(kinds == "features", "manhattan", "sumsq"), classic = rep("sumsq", length(kinds)))
    cl <- if (modes == "classic") rep(0, length(kinds)) else ifelse(kinds == "features", 1, 0)
    return(list(modes = unname(m), clip = unname(cl)))
  }
  list(modes = rep_len(modes, length(XL)), clip = rep_len(clip, length(XL)))
}

# ------------------------------------------------------------------------------ choosing the layers (nested)
# Greedy forward selection on the training children only: every layer alone, then the layer that most improves
# the leave-one-child-out 1-NN accuracy of the combined distance (every layer divided by its mean distance, as
# kohonen weighs layers), until no layer improves it. A cheap stand-in for training a map per subset; inside the
# held-out protocol it sees only the 72 training children of the fold.
select_layers <- function(Z, pool, modes, child, y) {
  rm <- resolve_modes(Z[pool], modes, 0)
  Dl <- lapply(seq_along(pool), function(i) {
    X <- Z[[pool[i]]]; sp <- layer_spec(setNames(list(X), pool[i]), rm$modes[i])
    D <- ss_cross(X, X, sp, 1); D / mean(D[upper.tri(D)])
  })
  same <- outer(child, child, "==")
  acc_of <- function(set) {
    D <- Reduce(`+`, Dl[set]); D[same] <- Inf
    mean(y[max.col(-D, ties.method = "first")] == y)
  }
  single <- vapply(seq_along(pool), acc_of, 0)
  chosen <- which.max(single); best <- single[chosen]; path <- data.frame(step = 1, add = pool[chosen], acc = best)
  repeat {
    rest <- setdiff(seq_along(pool), chosen)
    if (!length(rest)) break
    a <- vapply(rest, function(l) acc_of(c(chosen, l)), 0)
    if (max(a) <= best) break
    chosen <- c(chosen, rest[which.max(a)]); best <- max(a)
    path <- rbind(path, data.frame(step = nrow(path) + 1, add = pool[rest[which.max(a)]], acc = best))
  }
  list(layers = pool[sort(chosen)], single = setNames(single, pool), path = path)
}

# --------------------------------------------------------------------------------------------- one map
# The article's schedule (kohonen::som, the method's defaults): online, rlen passes, bubble neighbourhood, alpha 0.05 to 0.01, radius
# from the two-thirds quantile of the unit distances to 0 (0.5 below 1), the start K rows drawn with
# replacement (aweSOM::somInit "random"), the items drawn inside the loop as kohonen draws them.
# XL: the standardised layers; y: the task of every row (used only by a label layer); returns a fit with M (the
# codebook, X layers only), dist(A, B) (the map's own weighted distance), grid_dist, links, bmu, w, ...
fit_ss <- function(XL, meth, seed, xdim = 5, ydim = 5, rlen = 1000, y = NULL, trace = FALSE) {
  XL <- XL[meth$layers]
  if (meth$glue) {                                       # all layers as one feature vector (the control)
    XL <- list(glued = feature_layer(do.call(cbind, XL)))
  }
  rm <- resolve_modes(XL, meth$modes, meth$clip)
  spec <- layer_spec(XL, rm$modes, rm$clip)
  X <- do.call(cbind, lapply(XL, unclass)); N <- nrow(X); K <- xdim * ydim
  nl_x <- nrow(spec)
  if (meth$label_weight > 0) {                           # the label layer of a supervised map (xyf): one-hot tasks
    Y <- outer(y, 1:7, "==") * 1
    X <- cbind(X, Y)
    spec <- rbind(spec, data.frame(layer = "label", off = ncol(X) - 7L, len = 7L, mode = 0L, mode_name = "sumsq", L = 7L, nch = 1L,
                                   band = 1L, clip = 0, stringsAsFactors = FALSE))
  }
  S <- rlen * N; gd5 <- hex_grid(xdim, ydim)$dist; r0 <- unname(quantile(gd5, 2/3))
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  if (!meth$grow) {
    M0 <- X[sample(N, K, replace = TRUE), , drop = FALSE]
    w <- layer_weights(meth, M0, spec, nl_x)
    tr <- train_chunks(X, M0, gd5, spec, w, meth, S, r0, grow_steps = 0L, N = N, trace = trace, w0 = w)
    M <- tr$M; w <- tr$w; gdist <- gd5; pos <- kohonen::somgrid(xdim, ydim, "hexagonal")$pts; births <- NULL
  } else {
    g <- grow_ss(X, spec, meth, K, S, N, nl_x, trace = trace)
    M <- g$M; gdist <- g$map_dist; pos <- site_pos(g$sites); births <- g$births; tr <- g; w <- g$w
  }
  secs <- proc.time()[["elapsed"]] - t0
  keep <- seq_len(sum(spec$len[seq_len(nl_x)]))            # the X layers only
  sx <- spec[seq_len(nl_x), ]; wx <- w[seq_len(nl_x)] / sum(w[seq_len(nl_x)])
  fit <- list(method = meth$id, M = M[, keep, drop = FALSE], M_label = if (meth$label_weight > 0) M[, -keep, drop = FALSE] else NULL,
              spec = sx, w = wx, w_trace = tr$w_trace, grid_dist = gdist, links = abs(gdist - 1) < 1e-6, pos = pos, births = births,
              sumsq_only = all(sx$mode_name == "sumsq"), seconds = secs, layers = names(XL))
  fit$dist <- function(A, B) ss_cross(A, B, fit$spec, fit$w)
  fit$bmu <- max.col(-fit$dist(X[, keep, drop = FALSE], fit$M), ties.method = "first")
  fit
}

# the weights of the layers: kohonen's rule on the start rows (user weights over the distance scale of every
# layer); a label layer (xyf) gets the user weight label_weight and the X layers share the rest equally, as
# kohonen::xyf / supersom(user.weights = ...) would weigh them; "equal": no distance scaling at all
layer_weights <- function(meth, M0, spec, nl_x) {
  nl <- nrow(spec)
  user <- if (nl > nl_x) c(rep((1 - meth$label_weight) / nl_x, nl_x), meth$label_weight) else rep(1 / nl, nl)
  if (meth$weights == "equal") return(user)
  kohonen_weights(M0, spec, user)
}

# the loop of ss_train_cpp() over [s_from, s_to): one call per pass when the weights learn, else one call.
# Two rules that learn the layer weights, unsupervised, after every pass (both from a fresh mapping of all items):
#   "adapt"  w_l = len_l / Q_l, Q_l the mean distance of the items to their winner in layer l: the maximum-likelihood
#            scale of a layer whose residuals are Laplace (L1) or Gaussian (sum of squares) with a noise level of its
#            own; a layer the map fits well counts more (the inverse-variance weights of a weighted mean)
#   "agree"  w_l = w0_l / (1 + g_l), g_l the mean lattice distance between the layer's own nearest unit and the
#            winner, w0 kohonen's start weights: a layer that agrees with the consensus counts more, a layer of noise
#            disagrees with everything (an experimental weight rule of novel_SOM, not published)
# The label layer of a supervised map keeps its share.
LEARNED <- c("adapt", "agree")
train_chunks <- function(X, M, gdist, spec, w, meth, S, r0, grow_steps, N, s_from = 0L, s_to = S, owner = integer(nrow(X)),
                         err = numeric(nrow(X)), folds = numeric(nrow(M)), trace = FALSE, radius = c(r0, 0), w0 = w) {
  w_trace <- NULL
  cuts <- if (meth$weights %in% LEARNED) unique(c(seq(s_from, s_to, by = N), s_to)) else c(s_from, s_to)
  for (ci in seq_len(length(cuts) - 1)) {
    a <- cuts[ci]; b <- cuts[ci + 1]
    if (b <= a) next
    r <- ss_train_cpp(X, M, gdist, spec$off, spec$len, spec$mode, spec$L, spec$nch, spec$band, as.numeric(w), spec$clip,
                      as.integer(a), as.integer(b), as.integer(S), meth$alpha, radius, as.integer(grow_steps),
                      meth$bubble, TRUE, integer(0), as.integer(owner), as.numeric(err), as.numeric(folds), meth$float, 0.1)
    M <- r$M; owner <- r$owner; err <- r$err; folds <- r$folds
    if (meth$weights %in% LEARNED && b < s_to) {
      w <- learn_weights(X, M, spec, w, w0, gdist, meth$weights)
      if (trace) w_trace <- rbind(w_trace, c(step = b, w))
    }
  }
  list(M = M, owner = owner, err = err, folds = folds, w = w, w_trace = w_trace)
}

learn_weights <- function(X, M, spec, w, w0, gdist, rule) {
  Dl <- ss_cross(X, M, spec, w, per_layer = TRUE)
  Dc <- Reduce(`+`, Map(`*`, Dl, w)); win <- max.col(-Dc, ties.method = "first"); N <- nrow(X)
  lab <- spec$layer == "label"; ix <- which(!lab)
  wx <- if (rule == "adapt") {
    Q <- vapply(ix, function(l) mean(Dl[[l]][cbind(seq_len(N), win)]), 0)
    spec$len[ix] / pmax(Q, 1e-12)
  } else {
    g <- vapply(ix, function(l) mean(gdist[cbind(max.col(-Dl[[l]], ties.method = "first"), win)]), 0)
    w0[ix] / (1 + g)
  }
  out <- w; out[ix] <- wx / sum(wx) * (1 - sum(w[lab]))
  out
}

# --------------------------------------------------------------------------------------------- growing
# my earlier growing DTW-SOM with the SuperSOM's weighted distance: a unit and its six neighbours to start; births
# evenly spread over the first grow_share of the steps, each from the unit with the largest total distance to
# the items it last won; the child is the medoid (a real item, in every layer at once) of the half of the
# parent's items that the parent fits worst, placed on a free site next to the parent (the one with most
# occupied neighbours, then the closest prototypes); the parent's items that fit the child better move to it.
# The radius stays at grow_radius while the map grows, then falls to 0 (0.5 below 1); a unit whose items often
# have their runner-up off the lattice neighbours pulls with a wider radius (float).
grow_ss <- function(X, spec, meth, K, S, N, nl_x, trace = FALSE) {
  n0 <- min(meth$n_start, K); sites <- start_sites(n0)
  M <- X[sample(N, n0, replace = TRUE), , drop = FALSE]
  w <- w0 <- layer_weights(meth, M, spec, nl_x)
  map_dist <- hop_dist(sites)
  n_births <- K - n0; grow_steps <- if (n_births > 0) round(meth$grow_share * S) else 0L
  birth_at <- if (n_births > 0) pmax(1L, as.integer(round(grow_steps * seq_len(n_births) / n_births))) else integer(0)
  D0 <- ss_cross(X, M, spec, w); owner <- max.col(-D0, ties.method = "first"); err <- D0[cbind(seq_len(N), owner)]
  folds <- numeric(n0); births <- NULL; w_trace <- NULL; s <- 0L
  radius <- c(meth$grow_radius, 0)
  for (b in seq_along(birth_at)) {
    if (birth_at[b] > s) {
      tr <- train_chunks(X, M, map_dist, spec, w, meth, S, meth$grow_radius, grow_steps, N, s, birth_at[b], owner, err, folds, trace, radius, w0)
      M <- tr$M; owner <- tr$owner; err <- tr$err; folds <- tr$folds; w <- tr$w; w_trace <- rbind(w_trace, tr$w_trace); s <- birth_at[b]
    }
    q <- which.max(vapply(seq_len(nrow(M)), function(k) sum(err[owner == k]), 0))
    mine <- which(owner == q)
    child <- if (meth$birth == "copy" || length(mine) == 0) M[q, ] else {
      worse <- mine[err[mine] >= median(err[mine])]
      if (length(worse) == 1) X[worse, ] else { Dw <- ss_cross(X[worse, , drop = FALSE], X[worse, , drop = FALSE], spec, w); X[worse[which.min(rowSums(Dw))], ] }
    }
    free <- free_sites(sites)
    gap <- sqrt(colSums((t(site_pos(free)) - c(site_pos(sites[q, , drop = FALSE])))^2))
    free <- free[gap < min(gap) + 1e-6, , drop = FALSE]
    nb <- lapply(seq_len(nrow(free)), function(f) which(touching(free[f, ], sites)))
    fitn <- vapply(nb, function(k) mean(ss_cross(rbind(child), M[k, , drop = FALSE], spec, w)), 0)
    site <- free[order(-lengths(nb), fitn)[1], ]
    if (length(mine)) {
      d_new <- ss_cross(X[mine, , drop = FALSE], rbind(child), spec, w)[, 1]
      moved <- d_new < err[mine]; owner[mine[moved]] <- nrow(M) + 1L; err[mine[moved]] <- d_new[moved]
    } else moved <- logical(0)
    births <- rbind(births, data.frame(step = s, parent = q, moved = sum(moved)))
    M <- rbind(M, child, deparse.level = 0); sites <- rbind(sites, site, deparse.level = 0)
    map_dist <- hop_dist(sites); folds <- c(folds, 0)
  }
  tr <- train_chunks(X, M, map_dist, spec, w, meth, S, meth$grow_radius, grow_steps, N, s, S, owner, err, folds, trace, radius, w0)
  list(M = tr$M, sites = sites, map_dist = map_dist, births = births, w = tr$w, w_trace = rbind(w_trace, tr$w_trace))
}
