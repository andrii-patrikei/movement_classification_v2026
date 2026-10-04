# ss_report.R: the summaries, charts and tables of supersom.Rmd (plotly and DT, in the style of plots_interactive.R)
#
# Not part of the cache keys: editing this file never retrains a map. Colours: the reference palette's categorical
# slots in their fixed order, one per group of methods; the task colours of plots.R for the maps.

GROUPS <- c(single = "One layer", features = "Layered: the 4 feature sets", raw = "Layered: raw curves (DTW)",
            all = "Layered: all 7 layers", selected = "Layered: layers chosen per fold", supervised = "Layered + task layer (supervised)",
            separate = "Separate one-layer maps, voting")
GROUP_COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"), names(GROUPS))

method_table <- function() {
  m <- data.table::rbindlist(lapply(SS_METHODS, function(x) data.table::data.table(method = x$id, label = x$label, group = x$group,
                                                                               layers = paste(x$layers, collapse = "+"), grow = x$grow)))
  # the layered "Euclidean" controls put kohonen's sum of squares on the curves (one layer: the same winner as Euclidean)
  relabel <- c(ss_sel7_euc = "Robust SuperSOM, layers chosen, curves by sum of squares",
               ss_R3_euc = "Sum-of-squares SuperSOM, 3 raw curves")
  m[method %in% names(relabel), label := relabel[method]]
  rbind(m, data.table::data.table(method = c("vote_sep7", "vote_sep4"), label = c("Vote of 7 separate one-layer maps", "Vote of 4 separate feature-set maps"),
                                  group = "separate", layers = c(paste(POOL7, collapse = "+"), paste(F4, collapse = "+")), grow = FALSE))
}

# the vote of separate one-layer maps (late fusion), from the held-out predictions of the one-layer methods: the same
# seeds, so the same folds of children; a tie is broken at random (fixed seed)
separate_vote <- function(H, members, id) {
  set.seed(1)
  V <- data.table::dcast(H[method %in% members], seed + segment + child + y7 ~ method, value.var = c("pred_unit", "pred_super"))
  pick <- function(M) apply(M, 1, function(v) { t <- tabulate(v, 7L); b <- which(t == max(t)); if (length(b) == 1) b else sample(b, 1) })
  data.table::data.table(method = id, seed = V$seed, segment = V$segment, child = V$child, y7 = V$y7,
                         pred_unit = pick(as.matrix(V[, paste0("pred_unit_", members), with = FALSE])),
                         pred_super = pick(as.matrix(V[, paste0("pred_super_", members), with = FALSE])),
                         pred_vote = NA_integer_, acc_super = NA_integer_, acc_unit = NA_integer_)
}
with_votes <- function(H) {
  sep7 <- c("rob_rqa", "rob_autocorr", "rob_spectr", "rob_tsfresh", "dtw_posture", "dtw_accspec", "dtw_gyrospec")
  sep4 <- c("rob_rqa", "rob_autocorr", "rob_spectr", "rob_tsfresh")
  cols <- c("method", "seed", "segment", "child", "y7", "pred_unit", "pred_super", "pred_vote", "acc_super", "acc_unit")
  rbind(H[, cols, with = FALSE], separate_vote(H, sep7, "vote_sep7"), separate_vote(H, sep4, "vote_sep4"))
}

# per method and seed: the held-out accuracies (and, for maps with accelerometer layers, the same children mapped with
# those layers only); then per method: mean, sd, range and a 95% t interval over the seeds
seed_scores <- function(H) {
  H[, .(super = mean(pred_super == y7), unit = mean(pred_unit == y7), vote = mean(pred_vote == y7),
        super_bacc = balanced_accuracy(y7, pred_super, 1:7), unit_bacc = balanced_accuracy(y7, pred_unit, 1:7),
        acc_only_unit = if (all(is.na(acc_unit))) NA_real_ else mean(!is.na(acc_unit) & acc_unit == y7),
        acc_only_super = if (all(is.na(acc_super))) NA_real_ else mean(!is.na(acc_super) & acc_super == y7)), by = .(method, seed)]
}
summarise_seeds <- function(Sd, cols) {
  long <- data.table::melt(Sd, id.vars = c("method", "seed"), measure.vars = cols, variable.name = "measure")[!is.na(value)]
  long[, .(mean = mean(value), sd = sd(value), min = min(value), max = max(value), seeds = .N,
           lo = if (.N > 2 && sd(value) > 0) t.test(value)$conf.int[1] else mean(value),
           hi = if (.N > 2 && sd(value) > 0) t.test(value)$conf.int[2] else mean(value)), by = .(method, measure)]
}

pm <- function(x, d = 3) formatC(x, digits = d, format = "f")
leg_bottom <- function() list(orientation = "h", x = 0, xanchor = "left", y = 0, yref = "container", yanchor = "bottom",
                              font = list(color = INK2, size = 12))

# ------------------------------------------------------------------------------------- the verdict chart
# one row per map, the mean over the seeds with its range; buttons switch the rule (and re-sort the rows); the
# dashed line is the best one-layer map under that rule
p_verdict <- function(SM, mt, rules, rule_labels, button_labels = rule_labels, height = 860) {
  p <- plotly::plot_ly(height = height); vis <- list(); lay <- list(); n_per <- length(GROUPS); lo <- numeric(0)
  for (ri in seq_along(rules)) {
    d <- merge(SM[measure == rules[ri]], mt, by = "method")[order(mean)]
    best1 <- d[group == "single"][which.max(mean)]; lo[ri] <- floor(min(d$min) * 20) / 20
    for (g in names(GROUPS)) {
      x <- d[group == g]
      p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = x$mean, y = x$label, name = GROUPS[[g]], legendgroup = g,
                             showlegend = nrow(x) > 0,
                             error_x = list(type = "data", symmetric = FALSE, array = x$max - x$mean, arrayminus = x$mean - x$min,
                                            color = GROUP_COL[[g]], thickness = 1.4, width = 0),
                             marker = list(color = GROUP_COL[[g]], size = 11, line = list(color = "white", width = 2)),
                             text = sprintf("<b>%s</b><br>%s<br>mean %s over %d seeds (sd %s)<br>seeds from %s to %s<br>%+.1f points against %s",
                                            x$label, rule_labels[ri], pm(x$mean), x$seeds, pm(x$sd), pm(x$min), pm(x$max),
                                            100 * (x$mean - best1$mean), best1$label),
                             hovertemplate = "%{text}<extra></extra>", visible = ri == 1)
    }
    vis[[ri]] <- rep(seq_along(rules) == ri, each = n_per)
    lay[[ri]] <- list(yaxis = list(categoryarray = d$label, categoryorder = "array", type = "category", showgrid = FALSE,
                                   tickfont = list(color = INK2, size = 12), automargin = TRUE),
                      shapes = list(list(type = "line", x0 = best1$mean, x1 = best1$mean, y0 = 0, y1 = 1, yref = "paper",
                                         line = list(color = INK2, width = 1.5, dash = "dash"))),
                      annotations = list(list(x = best1$mean, y = 1, yref = "paper", yanchor = "bottom", xanchor = "right", showarrow = FALSE,
                                              text = sprintf("best one-layer map: %s", pm(best1$mean)), font = list(color = INK2, size = 11))))
  }
  xa <- function(i) p_axis(rule_labels[i], range = c(lo[i], 1.005), dtick = if (lo[i] >= 0.6) 0.05 else 0.1)
  p_finish(p, xaxis = xa(1), yaxis = lay[[1]]$yaxis, shapes = lay[[1]]$shapes, annotations = lay[[1]]$annotations,
           updatemenus = p_buttons(button_labels, vis, lapply(seq_along(rules), function(i) c(lay[[i]], list(xaxis = xa(i)))), y = 1.1),
           margin = list(t = 100, l = 10, r = 20, b = 150), legend = leg_bottom())
}

# ---------------------------------------------------------------------------- units pure, super-clusters not
# all data: in-sample unit accuracy (every unit names its task) against the article's super-cluster accuracy
p_purity_vs_super <- function(AS, mt) {
  d <- merge(AS, mt, by = "method")
  p <- plotly::plot_ly(height = 500)
  for (g in intersect(names(GROUPS), unique(d$group))) {
    x <- d[group == g]
    p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = x$unit_acc_in, y = x$accuracy, name = GROUPS[[g]],
                           marker = list(color = GROUP_COL[[g]], size = 11, line = list(color = "white", width = 2)),
                           text = sprintf("<b>%s</b><br>units name the task: %s<br>the article's 7 super-clusters: %s", x$label, pm(x$unit_acc_in), pm(x$accuracy)),
                           hovertemplate = "%{text}<extra></extra>")
  }
  p_finish(p, xaxis = p_axis("every unit names its task (all 648 segments)"), yaxis = p_axis("the article's rule: 7 super-clusters"),
           margin = list(t = 10, l = 60, r = 10, b = 140), legend = leg_bottom())
}

# --------------------------------------------------------------------------------------- child by child
# CC: child, a, b (per-child held-out accuracy of two maps, mean over the fresh seeds), one panel per rule (buttons)
p_child_bars <- function(CC, la, lb, rules, rule_labels, button_labels = rule_labels) {
  p <- plotly::plot_ly(height = 420); ann <- list(); vis <- list()
  for (ri in seq_along(rules)) {
    d <- CC[rule == rules[ri]][, diff := 100 * (b - a)][order(-diff, child)]
    col <- ifelse(d$diff > 0, GROUP_COL[["selected"]], ifelse(d$diff < 0, GROUP_COL[["single"]], GRID))
    p <- plotly::add_trace(p, type = "bar", x = seq_len(nrow(d)), y = d$diff, marker = list(color = col), showlegend = FALSE,
                           text = sprintf("<b>child %s</b><br>%s %s<br>%s %s<br>%+.1f points", d$child, la, pm(d$a), lb, pm(d$b), d$diff),
                           textposition = "none", hovertemplate = "%{text}<extra></extra>", visible = ri == 1)
    ann[[ri]] <- list(list(x = 0, y = 1.02, xref = "paper", yref = "paper", xanchor = "left", yanchor = "bottom", showarrow = FALSE,
                           text = sprintf("%s: layered better for <b>%d</b> children, worse for <b>%d</b>, equal for %d", rule_labels[ri],
                                          sum(d$diff > 0), sum(d$diff < 0), sum(d$diff == 0)), font = list(size = 13, color = INK)))
    vis[[ri]] <- seq_along(rules) == ri
  }
  p_finish(p, bargap = 0.15, xaxis = p_axis("children, sorted by the gain", showticklabels = FALSE, showgrid = FALSE),
           yaxis = p_axis(sprintf("%s minus %s (points)", lb, la), zeroline = TRUE, zerolinecolor = INK2),
           annotations = ann[[1]], margin = list(t = 90, l = 60, r = 10, b = 40),
           updatemenus = p_buttons(button_labels, vis, lapply(ann, function(a) list(annotations = a)), y = 1.25))
}

# ------------------------------------------------------------------------------------ a sensor goes missing
# one row per map that has accelerometer layers: held-out children with every layer, the same children with the
# accelerometer layers only, and the 9 children whose gyroscope failed (mapped with the accelerometer layers only)
p_missing <- function(MS) {
  p <- plotly::plot_ly(height = 140 + 34 * nrow(MS))
  parts <- list(c("all_layers", "held-out children, all layers", "circle", INK2),
                c("acc_only", "the same children, accelerometer layers only", "circle", GROUP_COL[["raw"]]),
                c("nogyro9", "the 9 children without a gyroscope (dropped by the article)", "diamond", GROUP_COL[["selected"]]))
  MS <- MS[order(nogyro9)]
  p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = as.vector(rbind(pmin(MS$acc_only, MS$all_layers, MS$nogyro9, na.rm = TRUE),
                                                                                 pmax(MS$acc_only, MS$all_layers, MS$nogyro9, na.rm = TRUE), NA)),
                         y = as.vector(rbind(MS$label, MS$label, NA)), line = list(color = GRID, width = 4), hoverinfo = "skip", showlegend = FALSE)
  for (pt in parts) {
    p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = MS[[pt[1]]], y = MS$label, name = pt[2],
                           marker = list(symbol = pt[3], color = pt[4], size = 12, line = list(color = "white", width = 2)),
                           text = sprintf("<b>%s</b><br>%s: %s", MS$label, pt[2], pm(MS[[pt[1]]])), hovertemplate = "%{text}<extra></extra>")
  }
  p_finish(p, xaxis = p_axis("accuracy, each unit names its task (mean over the seeds)"),
           yaxis = list(type = "category", categoryorder = "array", categoryarray = MS$label, showgrid = FALSE, automargin = TRUE,
                        tickfont = list(color = INK2, size = 12)),
           margin = list(t = 10, l = 10, r = 20, b = 120), legend = c(leg_bottom(), list(orientation = "v")))
}

# ------------------------------------------------------------------------------------ weights as they learn
p_weights <- function(WT, title_of) {
  p <- plotly::plot_ly(height = 460); ids <- unique(WT$method)
  dash <- setNames(c("solid", "dot")[seq_along(ids)], ids)
  layer_col <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7")[seq_along(POOL7)], POOL7)
  for (id in ids) for (l in unique(WT[method == id]$layer)) {
    d <- WT[method == id & layer == l]
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = d$pass, y = d$w, name = sprintf("%s (%s)", LAYER_LABELS[[l]], title_of[[id]]),
                           legendgroup = l, line = list(color = layer_col[[l]], width = 2, dash = dash[[id]]),
                           hovertemplate = paste0(LAYER_LABELS[[l]], ", ", title_of[[id]], "<br>pass %{x}: weight %{y:.3f}<extra></extra>"))
  }
  p_finish(p, xaxis = p_axis("pass (of 1000; 0 = kohonen's start weights)"), yaxis = p_axis("layer weight (sum 1)"),
           margin = list(t = 10, l = 60, r = 10, b = 150), legend = leg_bottom())
}

# ------------------------------------------------------------------------------------ the growing map, traced
# grow_ss() of ss_maps.R with a snapshot after every birth (the same code path, the same random numbers): the items
# each unit holds at that moment and where the units sit
grow_traced <- function(XL, meth, seed, rlen = 1000, K = 25) {
  XL <- XL[meth$layers]; rm <- resolve_modes(XL, meth$modes, meth$clip); spec <- layer_spec(XL, rm$modes, rm$clip)
  X <- do.call(cbind, lapply(XL, unclass)); N <- nrow(X); S <- rlen * N
  set.seed(seed)
  snaps <- list()
  env <- environment()
  hook <- function(M, sites, w, step, parent) {
    D <- ss_cross(X, M, spec, w); own <- max.col(-D, ties.method = "first")
    env$snaps[[length(env$snaps) + 1]] <- list(step = step, pos = site_pos(sites), owner = own, parent = parent, K = nrow(M))
  }
  g <- grow_ss_hooked(X, spec, meth, K, S, N, nrow(spec), hook)
  list(snaps = snaps, final = g)
}
# a copy of grow_ss() that calls hook() after every birth and at the end (kept apart so that ss_maps.R, which the
# cache keys hash, stays untouched); checked against grow_ss() in supersom.Rmd
grow_ss_hooked <- function(X, spec, meth, K, S, N, nl_x, hook) {
  n0 <- min(meth$n_start, K); sites <- start_sites(n0)
  M <- X[sample(N, n0, replace = TRUE), , drop = FALSE]
  w <- w0 <- layer_weights(meth, M, spec, nl_x)
  map_dist <- hop_dist(sites)
  n_births <- K - n0; grow_steps <- if (n_births > 0) round(meth$grow_share * S) else 0L
  birth_at <- if (n_births > 0) pmax(1L, as.integer(round(grow_steps * seq_len(n_births) / n_births))) else integer(0)
  D0 <- ss_cross(X, M, spec, w); owner <- max.col(-D0, ties.method = "first"); err <- D0[cbind(seq_len(N), owner)]
  folds <- numeric(n0); s <- 0L; radius <- c(meth$grow_radius, 0)
  hook(M, sites, w, 0L, NA)
  for (b in seq_along(birth_at)) {
    if (birth_at[b] > s) {
      tr <- train_chunks(X, M, map_dist, spec, w, meth, S, meth$grow_radius, grow_steps, N, s, birth_at[b], owner, err, folds, FALSE, radius, w0)
      M <- tr$M; owner <- tr$owner; err <- tr$err; folds <- tr$folds; w <- tr$w; s <- birth_at[b]
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
    }
    M <- rbind(M, child, deparse.level = 0); sites <- rbind(sites, site, deparse.level = 0)
    map_dist <- hop_dist(sites); folds <- c(folds, 0)
    hook(M, sites, w, s, q)
  }
  tr <- train_chunks(X, M, map_dist, spec, w, meth, S, meth$grow_radius, grow_steps, N, s, S, owner, err, folds, FALSE, radius, w0)
  hook(tr$M, sites, tr$w, S, NA)
  list(M = tr$M, sites = sites, map_dist = map_dist, w = tr$w)
}

# the growth as frames: every unit a hexagon coloured by the task most of its segments belong to (white: empty),
# its size by how many segments it holds; the slider steps through the births
p_growth <- function(G, y7, steps_total) {
  hexpts <- function(x, y, r = 0.5 / cos(pi / 6)) { a <- seq(30, 330, by = 60) * pi / 180; cbind(x + r * cos(a), y + r * sin(a)) }
  frames <- list(); xr <- range(unlist(lapply(G$snaps, function(s) s$pos[, 1]))) + c(-1, 1)
  yr <- range(unlist(lapply(G$snaps, function(s) s$pos[, 2]))) + c(-1, 1)
  rows <- list()
  for (fi in seq_along(G$snaps)) {
    s <- G$snaps[[fi]]; K <- s$K
    tab <- table(factor(s$owner, levels = seq_len(K)), factor(y7, levels = 1:7)); n <- rowSums(tab)
    maj <- ifelse(n > 0, apply(tab, 1, which.max), NA)
    pur <- ifelse(n > 0, apply(tab, 1, max) / pmax(n, 1), NA)
    for (u in seq_len(K)) {
      h <- hexpts(s$pos[u, 1], s$pos[u, 2])
      rows[[length(rows) + 1]] <- data.frame(frame = fi, unit = u, x = c(h[, 1], h[1, 1], NA), y = c(h[, 2], h[1, 2], NA),
                                             task = if (is.na(maj[u])) 0L else maj[u], n = n[u], pur = pur[u],
                                             newborn = !is.na(s$parent) && u == K)
    }
  }
  R <- do.call(rbind, rows)
  lab <- vapply(G$snaps, function(s) sprintf("%d units, step %s of %s", s$K, format(s$step, big.mark = ","), format(steps_total, big.mark = ",")), "")
  # one trace per task colour (fill), plus the centres with hover text; frames replace their data
  mk <- function(fi) {
    d <- R[R$frame == fi, ]
    tr <- lapply(0:7, function(k) {
      e <- d[d$task == k, ]
      list(type = "scatter", mode = "lines", x = if (nrow(e)) e$x else c(NA), y = if (nrow(e)) e$y else c(NA), fill = "toself",
           fillcolor = if (k == 0) "#ffffff" else TASK_COL[[k]], line = list(color = "white", width = 1.5), hoverinfo = "skip",
           name = if (k == 0) "empty unit" else TASKS7[[k]], showlegend = TRUE)
    })
    s <- G$snaps[[fi]]
    cen <- data.frame(x = s$pos[, 1], y = s$pos[, 2]); tab <- table(factor(s$owner, levels = seq_len(s$K)), factor(y7, levels = 1:7))
    txt <- vapply(seq_len(s$K), function(u) { n <- sum(tab[u, ]); if (!n) return(sprintf("unit %d: empty", u))
      sprintf("unit %d: %d segments<br>%s", u, n, paste(sprintf("%s %d", TASKS7[tab[u, ] > 0], tab[u, tab[u, ] > 0]), collapse = "<br>")) }, "")
    born <- if (!is.na(s$parent)) s$K else integer(0)
    c(tr, list(list(type = "scatter", mode = "markers", x = cen$x, y = cen$y, text = txt, hovertemplate = "%{text}<extra></extra>",
                    marker = list(size = 26, color = "rgba(0,0,0,0)"), showlegend = FALSE),
               list(type = "scatter", mode = "markers", x = I(if (length(born)) cen$x[born] else NA), y = I(if (length(born)) cen$y[born] else NA),
                    marker = list(symbol = "star", size = 14, color = INK, line = list(color = "white", width = 1)), name = "newborn unit",
                    hoverinfo = "skip", showlegend = TRUE)))
  }
  first <- mk(1)
  p <- plotly::plot_ly(height = 520)
  for (t in first) p <- do.call(plotly::add_trace, c(list(p), t))
  p$x$frames <- lapply(seq_along(G$snaps), function(fi) list(name = lab[fi], data = mk(fi)))
  steps <- lapply(seq_along(lab), function(fi) list(label = as.character(G$snaps[[fi]]$K), method = "animate",
                                                    args = list(list(lab[fi]), list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                                                                                    transition = list(duration = 0)))))
  p_finish(p, xaxis = list(visible = FALSE, range = xr, fixedrange = TRUE, constrain = "domain"),
           yaxis = list(visible = FALSE, range = yr, scaleanchor = "x", scaleratio = 1, fixedrange = TRUE, constrain = "domain"),
           sliders = list(list(active = 0, steps = steps, x = 0.05, len = 0.9, y = 0, currentvalue = list(prefix = "units on the map: ", font = list(color = INK2)),
                               pad = list(t = 30))),
           updatemenus = list(list(type = "buttons", showactive = FALSE, x = 0, y = 0, xanchor = "right", yanchor = "top", pad = list(t = 30, r = 10),
                                   buttons = list(list(label = "Play", method = "animate",
                                                       args = list(NULL, list(fromcurrent = TRUE, frame = list(duration = 450, redraw = TRUE), transition = list(duration = 0))))))),
           margin = list(t = 10, l = 10, r = 10, b = 90), legend = list(orientation = "h", x = 0, y = 1.08, font = list(color = INK2, size = 11)))
}

# ------------------------------------------------------------------------------------------ the layers alone
# how well every layer separates the tasks by itself, before any map: leave-one-child-out 1-NN accuracy (a segment
# takes the task of its nearest segment of another child), with the Euclidean distance and with the layer's own
# distance in the robust maps (Manhattan on features, DTW on curves)
layer_nn <- function(S) {
  pred_of <- function(X, mode) {
    sp <- layer_spec(list(a = X), mode); D <- ss_cross(X, X, sp, 1); D[outer(S$child, S$child, "==")] <- Inf
    S$y7[max.col(-D, ties.method = "first")]
  }
  ids <- c("rqa", "autocorr", "spectr", "mix", "tsfresh", R3)
  P <- lapply(setNames(ids, ids), function(l) { X <- S$XL[[l]]; own <- if (attr(X, "kind") == "features") "manhattan" else "dtw"
    list(euclid = pred_of(X, "sumsq"), own = pred_of(X, own)) })
  tab <- data.table::rbindlist(lapply(ids, function(l) data.table::data.table(layer = l, kind = attr(S$XL[[l]], "kind"), columns = ncol(S$XL[[l]]),
                                                                            euclid = mean(P[[l]]$euclid == S$y7), own = mean(P[[l]]$own == S$y7))))
  list(tab = tab, pred = lapply(P, `[[`, "own"))
}
