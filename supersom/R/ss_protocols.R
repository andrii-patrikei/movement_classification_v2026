# ss_protocols.R: how a SuperSOM is judged in supersom.Rmd; the same protocols as R/protocols.R of som_variants.Rmd
#
#   all data   one map on all 648 segments: the article's accuracy (7 PAM super-clusters, each named after its
#              majority task), ARI and NMI of the super-clusters, the unit purity, and the map's quality:
#              QE and TE in the map's own weighted distance, the topographic product, and aweSOM's Euclidean QE,
#              EV, TE and KL on the standardised features side by side (comparable between maps on the same layers)
#   held-out   9 folds of 9 children: map, super-clusters and their tasks from 72 children, every segment of the
#              other 9 predicted by the super-cluster of its winner (the article's rule) or by the winner itself
#              (unit majority); the layers are standardised on the 72 training children
# Uses article_accuracy(), som_quality(), ari(), balanced_accuracy(), child_folds() and unit_tasks() of the
# som_variants.Rmd code (R/article_metrics.R, R/protocols.R).

# the 7 super-clusters: PAM on the prototypes. A map of sum-of-squares layers is Euclidean once every layer is
# multiplied by the square root of its weight, so PAM runs on those scaled prototypes, which for one layer is
# literally the article's cluster::pam(codes, 7); otherwise PAM on the map's own distances between prototypes
super_clusters_ss <- function(fit, k = 7L) {
  if (fit$sumsq_only) {
    sc <- rep(sqrt(fit$w), fit$spec$len)
    Ms <- if (nrow(fit$spec) == 1) fit$M else sweep(fit$M, 2, sc, "*")
    return(cluster::pam(Ms, k)$clustering)
  }
  D <- fit$dist(fit$M, fit$M); D <- pmax((D + t(D)) / 2, 0); diag(D) <- 0
  cluster::pam(stats::as.dist(D), k, diss = TRUE)$clustering
}

nmi <- function(a, b) {
  tab <- table(a, b) / length(a); pa <- rowSums(tab); pb <- colSums(tab)
  h <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  mi <- sum(tab[tab > 0] * log(tab[tab > 0] / outer(pa, pb)[tab > 0]))
  den <- (h(pa) + h(pb)) / 2
  if (den == 0) 0 else mi / den
}
purity_of <- function(bmu, y) sum(tapply(y, bmu, function(v) max(table(v)))) / length(y)

# quality of a fit on the matrix it was trained on (Xc: the X layers side by side)
ss_quality <- function(fit, Xc, euclid_space = TRUE) {
  D <- fit$dist(Xc, fit$M); N <- nrow(Xc)
  bmu <- max.col(-D, ties.method = "first"); D2 <- D; D2[cbind(seq_len(N), bmu)] <- Inf
  bmu2 <- max.col(-D2, ties.method = "first")
  dV <- fit$dist(fit$M, fit$M); dV <- pmax((dV + t(dV)) / 2, 0); if (fit$sumsq_only) dV <- sqrt(dV)
  own <- list(QE_own = mean(D[cbind(seq_len(N), bmu)]), TE_own = mean(!fit$links[cbind(bmu, bmu2)]),
              TP = topographic_product_cpp(dV, fit$grid_dist), units_used = length(unique(bmu)) / nrow(fit$M))
  if (euclid_space) {                                    # aweSOM's measures in the Euclidean space of the features
    q <- som_quality(Xc, fit$M, fit$grid_dist, fit$links)
    own <- c(own, list(QE = q$QE, EV = q$EV, TE = q$TE, KL = q$KL))
  }
  own
}

# a method with select = TRUE picks its layers from meth$layers (the pool) first, on the rows it may see
with_selection <- function(meth, Z, child, y) {
  if (!isTRUE(meth$select)) return(meth)
  sel <- select_layers(Z, meth$layers, meth$modes, child, y)
  meth$layers <- sel$layers
  meth
}

eval_all_data_ss <- function(XL, y7, meth, seed, rlen = 1000, child = NULL, no_gyro = NULL, raw_scale = NULL) {
  meth <- with_selection(meth, XL, child, y7)
  f <- fit_ss(XL, meth, seed, rlen = rlen, y = y7)
  sup <- super_clusters_ss(f)
  a <- article_accuracy(f$bmu, sup, y7)
  ng <- data.frame(nogyro9_super = NA_real_, nogyro9_unit = NA_real_)
  use <- intersect(meth$layers, ACC_LAYERS)
  if (!is.null(no_gyro) && length(use)) {                 # the 9 children the article dropped: accelerometer layers only
    Zg <- lapply(setNames(use, use), function(l) scale_with(no_gyro[[l]], raw_scale[[l]]))
    b9 <- bmu_with(f, fill_layers(f, Zg), use)
    ng$nogyro9_super <- mean(a$task[sup[b9]] == no_gyro$y7); ng$nogyro9_unit <- mean(unit_tasks(f, y7)[b9] == no_gyro$y7)
  }
  Xc <- do.call(cbind, lapply(XL[meth$layers], unclass))
  q <- ss_quality(f, Xc, euclid_space = all(vapply(XL[meth$layers], function(X) attr(X, "kind") == "features", TRUE)))
  ut <- unit_tasks(f, y7)
  data.frame(accuracy = a$accuracy, mnLogLoss = a$mnLogLoss, ari_super = ari(sup[f$bmu], y7), nmi_super = nmi(sup[f$bmu], y7),
             purity = purity_of(f$bmu, y7), unit_acc_in = mean(ut[f$bmu] == y7), q, ng, K = nrow(f$M), seconds = f$seconds,
             w = paste(sprintf("%.3f", f$w), collapse = " "), layers = paste(meth$layers, collapse = "+"))
}

# the layer vote of a layered map: every layer finds its own nearest unit (that layer's distance alone, as
# kohonen's map(whatmap = l)), each such unit names its task (unit majority of the training items, mapped with
# all layers), and the layers vote; a tie goes to the task of the unit that all layers together choose
layer_vote <- function(f, Zte, unit_task, bmu_all) {
  nl <- nrow(f$spec)
  if (nl == 1) return(unit_task[bmu_all])
  Dl <- ss_cross(Zte, f$M, f$spec, f$w, per_layer = TRUE)
  V <- vapply(Dl, function(D) unit_task[max.col(-D, ties.method = "first")], integer(nrow(Zte)))
  V <- matrix(V, nrow = nrow(Zte))
  vapply(seq_len(nrow(V)), function(i) {
    t <- tabulate(V[i, ], 7L); best <- which(t == max(t))
    if (length(best) == 1) best else unit_task[bmu_all[i]]
  }, 1L)
}

heldout_predictions_ss <- function(XL, y7, child, meth, seed, k = 9, rlen = 1000) {
  folds <- child_folds(child, k, seed)
  N <- nrow(XL[[1]]); pred_super <- pred_unit <- pred_vote <- integer(N); acc_super <- acc_unit <- rep(NA_integer_, N)
  wsum <- list(); chosen <- character(0)
  for (te in folds) {
    tr <- setdiff(seq_len(N), te)
    Z <- standardize_layers(XL[meth$layers], tr)
    Ztr <- lapply(Z, function(X) keep_attr(X[tr, , drop = FALSE], X))
    mf <- with_selection(meth, Ztr, child[tr], y7[tr])     # the layers, chosen from the 72 training children only
    Z <- Z[mf$layers]; Ztr <- Ztr[mf$layers]; chosen <- c(chosen, paste(mf$layers, collapse = "+"))
    f <- fit_ss(Ztr, mf, seed, rlen = rlen, y = y7[tr])
    sup <- super_clusters_ss(f)
    task_super <- article_accuracy(f$bmu, sup, y7[tr])$task
    Zte <- do.call(cbind, lapply(Z, function(X) unclass(X)[te, , drop = FALSE]))
    bmu_te <- max.col(-f$dist(Zte, f$M), ties.method = "first")
    ut <- unit_tasks(f, y7[tr])
    pred_super[te] <- task_super[sup[bmu_te]]
    pred_unit[te] <- ut[bmu_te]
    pred_vote[te] <- layer_vote(f, Zte, ut, bmu_te)
    use <- intersect(mf$layers, ACC_LAYERS)
    if (length(use)) {                                    # the same children as if their gyroscope had failed
      b_acc <- bmu_with(f, Zte, use); acc_super[te] <- task_super[sup[b_acc]]; acc_unit[te] <- ut[b_acc]
    }
    wsum[[length(wsum) + 1]] <- f$w
  }
  out <- data.frame(segment = seq_len(N), child = child, y7 = y7, pred_super = pred_super, pred_unit = pred_unit, pred_vote = pred_vote,
                    acc_super = acc_super, acc_unit = acc_unit)
  attr(out, "w") <- if (length(unique(lengths(wsum))) == 1) colMeans(do.call(rbind, wsum)) else NA
  attr(out, "layers") <- paste(names(sort(table(chosen), decreasing = TRUE)), sort(table(chosen), decreasing = TRUE), sep = " x", collapse = "; ")
  out
}

eval_heldout_ss <- function(XL, y7, child, meth, seed, k = 9, rlen = 1000) {
  p <- heldout_predictions_ss(XL, y7, child, meth, seed, k, rlen)
  data.frame(heldout_super_acc = mean(p$pred_super == y7), heldout_super_bacc = balanced_accuracy(y7, p$pred_super, 1:7),
             heldout_unit_acc = mean(p$pred_unit == y7), heldout_unit_bacc = balanced_accuracy(y7, p$pred_unit, 1:7),
             heldout_vote_acc = mean(p$pred_vote == y7), w_heldout = paste(sprintf("%.3f", attr(p, "w")), collapse = " "),
             layers_heldout = attr(p, "layers"))
}

# ------------------------------------------------------------------------------------ a sensor goes missing
# kohonen's map(whatmap = ...): an item is mapped with some layers only, the other layers' weights set to zero
bmu_with <- function(f, Xfull, use) {
  w <- f$w * (f$spec$layer %in% use)
  max.col(-ss_cross(Xfull, f$M, f$spec, w), ties.method = "first")
}
# curves of new children scaled with the per-channel means and sds of the reference children's raw curves
scale_with <- function(Xn, Xref) {
  L <- attr(Xref, "L"); Z <- Xn
  for (ch in seq_len(attr(Xref, "nch"))) {
    cols <- (ch - 1) * L + seq_len(L); v <- as.vector(Xref[, cols])
    Z[, cols] <- (Xn[, cols] - mean(v)) / sd(v)
  }
  keep_attr(Z, Xref)
}
# the item matrix of a fit's layers, the layers in 'use' filled in and the others left at zero (never read)
fill_layers <- function(f, Zuse) {
  X <- matrix(0, nrow(Zuse[[1]]), sum(f$spec$len))
  for (l in names(Zuse)) { i <- match(l, f$spec$layer); X[, f$spec$off[i] + seq_len(f$spec$len[i])] <- unclass(Zuse[[l]]) }
  X
}

ACC_LAYERS <- c("posture", "accspec")                   # the accelerometer curves: what a child without a gyroscope still has
