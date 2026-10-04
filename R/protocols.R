# protocols.R: the three ways a map is judged in som_variants.Rmd
#
# (a) all data      the article's Table 4: one map on all 648 segments; QE, EV, TE, KL (Euclidean, as
#                   aweSOM::somQuality), the article's accuracy and mnLogLoss, and the ARI of the 7 super-clusters
# (b) article k-fold the article's own k-fold, as train_som_cv() does it: 4 folds stratified by the 8 segment
#                   labels; "train" = a map on 3/4 judged on those 3/4; "test" = a new map trained on the
#                   remaining 1/4 and judged on that 1/4 (so neither number is a prediction for unseen data)
# (c) held-out      what the article's test was meant to be: 9 folds of 9 children; the map, its super-clusters
#                   and their tasks come from the 72 training children; every segment of the 9 unseen children is
#                   sent to its winner and predicted as the task of the winner's super-cluster (the article's rule)
#                   or of the winner itself (unit majority; an empty unit takes the task of the nearest used unit)
# heldout_predictions() keeps the prediction of every segment, for the comparison child by child

standardize <- function(X) { Z <- scale(X); Z[!is.finite(Z)] <- 0; Z }

# accuracy uses the map's own distance for the winners and the super-clusters; accuracy_euclid reads the same
# codebook the article's way (Euclidean winners, cluster::pam(codes, 7)), so a gain cannot come from the
# clustering step alone (on the Euclidean maps the two are the same number)
eval_all_data <- function(X, y7, method, seed, rlen_novel = 10) {
  f <- fit_som(X, method, seed, rlen_novel)
  sup <- super_clusters(f)
  a <- article_accuracy(f$bmu, sup, y7)
  q <- som_quality(X, f$M, f$grid_dist, f$links)
  a_e <- article_accuracy(q$bmu, cluster::pam(f$M, 7)$clustering, y7)
  data.frame(QE = q$QE, EV = q$EV, TE = q$TE, KL = q$KL, accuracy = a$accuracy, mnLogLoss = a$mnLogLoss,
             accuracy_euclid = a_e$accuracy, ari_super = ari(sup[f$bmu], y7), units_used = length(unique(f$bmu)) / nrow(f$M),
             seconds = f$seconds)
}

eval_article_kfold <- function(X, y8, y7, method, seed, k = 4, rlen_novel = 10) {
  folds <- stratified_folds(y8, k, seed)
  acc <- t(sapply(folds, function(te) {
    tr <- setdiff(seq_len(nrow(X)), te)
    f_tr <- fit_som(standardize(X[tr, , drop = FALSE]), method, seed, rlen_novel)
    f_te <- fit_som(standardize(X[te, , drop = FALSE]), method, seed, rlen_novel)
    c(train = article_accuracy(f_tr$bmu, super_clusters(f_tr), y7[tr])$accuracy,
      test = article_accuracy(f_te$bmu, super_clusters(f_te), y7[te])$accuracy)
  }))
  data.frame(kfold_train = mean(acc[, "train"]), kfold_test = mean(acc[, "test"]))
}

# every unit named after the majority task of its training items; an empty unit takes the task of the
# nearest used unit (own distance between prototypes)
unit_tasks <- function(f, y7) {
  K <- nrow(f$M); task <- rep(NA_integer_, K)
  maj <- tapply(y7, f$bmu, function(v) as.integer(names(which.max(table(v)))))
  task[as.integer(names(maj))] <- maj
  if (anyNA(task)) {
    D <- f$dist(f$M, f$M); if (isTRUE(f$softdtw)) D <- D - outer(diag(D), diag(D), "+") / 2
    used <- which(!is.na(task))
    for (u in which(is.na(task))) task[u] <- task[used[which.min(D[u, used])]]
  }
  task
}

# the prediction of every segment when its child is held out (the training part is standardised on its own)
heldout_predictions <- function(X, y7, child, method, seed, k = 9, rlen_novel = 10) {
  folds <- child_folds(child, k, seed)
  pred_super <- pred_unit <- integer(nrow(X))
  for (te in folds) {
    tr <- setdiff(seq_len(nrow(X)), te)
    mu <- colMeans(X[tr, , drop = FALSE]); sdv <- apply(X[tr, , drop = FALSE], 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    Ztr <- scale(X[tr, , drop = FALSE], mu, sdv); Zte <- scale(X[te, , drop = FALSE], mu, sdv)
    f <- fit_som(Ztr, method, seed, rlen_novel)
    sup <- super_clusters(f)
    task_super <- article_accuracy(f$bmu, sup, y7[tr])$task
    bmu_te <- map_bmu(f, Zte)
    pred_super[te] <- task_super[sup[bmu_te]]
    pred_unit[te] <- unit_tasks(f, y7[tr])[bmu_te]
  }
  data.frame(segment = seq_len(nrow(X)), child = child, y7 = y7, pred_super = pred_super, pred_unit = pred_unit)
}

eval_heldout <- function(X, y7, child, method, seed, k = 9, rlen_novel = 10) {
  p <- heldout_predictions(X, y7, child, method, seed, k, rlen_novel)
  data.frame(heldout_super_acc = mean(p$pred_super == y7), heldout_super_bacc = balanced_accuracy(y7, p$pred_super, 1:7),
             heldout_unit_acc = mean(p$pred_unit == y7), heldout_unit_bacc = balanced_accuracy(y7, p$pred_unit, 1:7))
}
