# article_metrics.R: the measures of the Acta Gymnica article, ported without caret
#
# Patrikei et al. (2026), code at https://code.it4i.cz/ADAS/movement-classification, functions/functions_df_total_cv.R:
#   get_som_accuracy()    7 PAM super-clusters of the codebook, each named after the task with most members,
#                         every segment predicted as the task of its super-cluster (crawl 1 and 2 are one class)
#   aweSOM::somQuality()  quantization error, % explained variance, topographic error, Kaski-Lagus error
# The ports are checked against the numbers saved with the article's models (trained_models/*.RDS) in
# som_variants.Rmd, section 1.

# the 8 segments of a child, in the order of the data (the activity labels of the Zenodo data set): 4x10m Shuttle
# Run, Crawling 1, Progressive Cone Shuffle Run, Ring Collection and Placement Run, 10x Position Change Drill,
# 10x Gymnastic Hoops Passing, 10x Bench Jumping, Crawling 2; the target has 7 classes because the two crawls
# are one task
TASKS7 <- c("1" = "4x10 m shuttle run", "2" = "Crawling", "3" = "Progressive cone shuffle run",
            "4" = "Ring collection and placement run", "5" = "10x position change drill",
            "6" = "10x gymnastic hoops passing", "7" = "10x bench jumping")
task7_of <- function(seg8) { y <- as.integer(seg8); y[y == 8L] <- 2L; y }

# --------------------------------------------------------------------------------- the article's accuracy
# super: the super-cluster of every unit (1..n_super); bmu: the unit of every segment; y7: 1..7.
# Returns the super-cluster x class table, the task of every super-cluster, the accuracy and the article's
# mnLogLoss (see below).
article_accuracy <- function(bmu, super, y7, n_super = 7L) {
  tab <- table(factor(super[bmu], levels = seq_len(n_super)), factor(y7, levels = 1:7))
  task <- apply(tab, 1, which.max)                     # ties: the first class, as which.max in the article
  pred <- task[super[bmu]]
  list(table = tab, task = task, pred = unname(pred), accuracy = mean(pred == y7), mnLogLoss = article_mnlogloss(tab, task))
}

# The article's mnLogLoss: the super-clusters sorted by their task, every class column divided by the class
# total, and row i read as a prediction of class i (caret::mnLogLoss -> ModelMetrics::mlogLoss: clipped to
# [1e-15, 1 - 1e-15]). It is the mean of -log(share of class i that lands in the i-th sorted super-cluster),
# so a class that no super-cluster is named after costs log(1e15) / 7 = 4.9. Kept for the article's tables.
article_mnlogloss <- function(tab, task = apply(tab, 1, which.max)) {
  ord <- order(task)                                   # dplyr::arrange is stable, so is order()
  P <- sweep(unclass(tab)[ord, , drop = FALSE], 2, pmax(colSums(tab), 1), "/")
  k <- min(nrow(P), ncol(P))
  p <- pmin(pmax(P[cbind(seq_len(k), seq_len(k))], 1e-15), 1 - 1e-15)
  -mean(log(p))
}

# ------------------------------------------------------------------------------------- the map's quality
# aweSOM::somQuality on any map: X the (scaled) data, M the codebook, grid_dist the distances between the
# units on the lattice, links 1 where two units are direct neighbours. Euclidean throughout, as in the
# article: the winner and the runner-up of every item are the two nearest prototypes in Euclidean distance.
#   QE  mean squared distance to the winner (aweSOM's err.quant is the mean of squared distances)
#   EV  100 - 100 QE / total variance, rounded to 2 decimals as aweSOM does
#   TE  share of items whose winner and runner-up are not neighbours on the lattice
#   KL  Kaski-Lagus: distance to the winner + the shortest way from winner to runner-up along the lattice,
#       every step as long as the distance between the two prototypes
som_quality <- function(X, M, grid_dist, links = round(grid_dist, 3) == 1) {
  D2 <- pmax(outer(rowSums(X^2), rowSums(M^2), "+") - 2 * X %*% t(M), 0)
  bmu <- max.col(-D2, ties.method = "first")
  D2b <- D2; D2b[cbind(seq_len(nrow(X)), bmu)] <- Inf
  bmu2 <- max.col(-D2b, ties.method = "first")
  sq <- D2[cbind(seq_len(nrow(X)), bmu)]
  QE <- mean(sq)
  EV <- 100 - round(100 * QE / mean(rowSums(sweep(X, 2, colMeans(X))^2)), 2)
  TE <- mean(!links[cbind(bmu, bmu2)])
  dV <- as.matrix(dist(M))
  ways <- map_ways(dV, links)
  KL <- mean(ways[cbind(bmu, bmu2)] + sqrt(sq))
  list(QE = QE, EV = EV, TE = TE, KL = KL, bmu = bmu)
}

# shortest ways along the lattice (Floyd-Warshall), as e1071::allShortestPaths in aweSOM
map_ways <- function(dV, links) {
  K <- nrow(dV); P <- ifelse(links, dV, Inf); diag(P) <- 0
  for (k in seq_len(K)) P <- pmin(P, outer(P[, k], P[k, ], "+"))
  P
}

# --------------------------------------------------------------------------------------- class agreement
ari <- function(a, b) {
  tab <- table(a, b); pairs <- function(n) n * (n - 1) / 2
  both <- sum(pairs(tab)); in_a <- sum(pairs(rowSums(tab))); in_b <- sum(pairs(colSums(tab)))
  chance <- in_a * in_b / pairs(length(a)); den <- (in_a + in_b) / 2 - chance
  if (den == 0) return(0)
  (both - chance) / den
}
balanced_accuracy <- function(truth, pred, classes = sort(unique(truth))) {
  mean(sapply(classes, function(k) mean(pred[truth == k] == k)))
}

# --------------------------------------------------------------------------------------- the article's folds
# caret::createFolds(y, k) for a factor y: within every class the items are dealt to the k folds in a random
# order (balanced), so every fold holds about 1/k of every class. Same scheme, own random stream.
stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  for (cl in unique(y)) { idx <- which(y == cl); fold[idx] <- sample(rep_len(seq_len(k), length(idx))) }
  lapply(seq_len(k), function(f) which(fold == f))
}
# children dealt to k folds (the held-out protocol): no child is in the training and the test part at once
child_folds <- function(child, k, seed) {
  set.seed(seed)
  kids <- sample(unique(child)); f <- rep_len(seq_len(k), length(kids))
  lapply(seq_len(k), function(i) which(child %in% kids[f == i]))
}
