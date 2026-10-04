# plots.R: chart style and the figures of som_variants.Rmd
#
# One colour per training regime (the categorical slots 1-3 of the reference palette, validated as a set),
# an open marker for the elastic variants (no time axis in a feature vector), black for the article's own
# published map. Seven task colours for the maps (slots 1-7, in order).

COL <- c(kohonen = "#2a78d6", article = "#eb6834", novel = "#1baf7a", published = "#0b0b0b")
REGIME_LAB <- c(kohonen = "Article's SOM (kohonen)", article = "novel_SOM variant, article's schedule",
                novel = "novel_SOM variant, its own schedule")
TASK_COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"), TASKS7)
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e0"

theme_bench <- function(base = 11) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
                   axis.text = ggplot2::element_text(colour = INK2), axis.title = ggplot2::element_text(colour = INK2),
                   strip.text = ggplot2::element_text(colour = INK, face = "bold", hjust = 0),
                   plot.title = ggplot2::element_text(colour = INK, face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour = INK2), plot.caption = ggplot2::element_text(colour = INK2, hjust = 0),
                   legend.position = "top", legend.justification = "left", legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(colour = INK2))
}

# the method string of a job -> variant id, regime, and a readable label
method_parts <- function(method) {
  v <- sub(":.*$", "", method); regime <- ifelse(method == "kohonen", "kohonen", ifelse(grepl(":article$", method), "article", "novel"))
  label <- ifelse(v == "kohonen", "kohonen::som (the article)", vapply(v, function(x) if (x %in% names(VARIANTS)) VARIANTS[[x]]$label else x, ""))
  data.frame(variant = v, regime = regime, label = unname(label), elastic = v %in% ELASTIC, stringsAsFactors = FALSE)
}

# dot plot: one row per variant, one dot per regime (mean over the seeds) with the range of the seeds,
# one panel per feature set; the article's published value (seed 1245) as a black diamond.
# S: one measure, columns label, regime, elastic, set_label, mean, min, max
dot_plot <- function(S, published = NULL, xlab = NULL, title = NULL, subtitle = NULL, caption = NULL) {
  S <- as.data.frame(S); S <- S[!is.na(S$mean), ]
  ord <- aggregate(S$mean, list(label = S$label), mean)
  S$label <- factor(S$label, levels = ord$label[order(ord$x)])
  S$y <- as.integer(S$label) + dodge_y(S$label, S$regime)
  g <- ggplot2::ggplot(S, ggplot2::aes(y = y)) +
    ggplot2::geom_segment(ggplot2::aes(x = min, xend = max, yend = y, colour = regime), linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = mean, colour = regime, shape = elastic), size = 2.2, stroke = 0.9) +
    ggplot2::scale_colour_manual(values = COL[c("kohonen", "article", "novel")], labels = REGIME_LAB, breaks = names(REGIME_LAB)) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "distance meaningful on features", `TRUE` = "elastic distance (no time axis in features)")) +
    ggplot2::scale_y_continuous(breaks = seq_along(levels(S$label)), labels = levels(S$label), expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
    ggplot2::facet_wrap(~ set_label, nrow = 1, scales = "free_x") +
    ggplot2::labs(x = xlab, y = NULL, title = title, subtitle = subtitle, caption = caption) +
    theme_bench() + ggplot2::theme(legend.box = "vertical", panel.grid.major.y = ggplot2::element_blank(),
                                   panel.spacing.x = ggplot2::unit(1.4, "lines"))
  if (!is.null(published)) {
    published$y <- match("kohonen::som (the article)", levels(S$label))
    g <- g + ggplot2::geom_point(data = published, ggplot2::aes(x = value, y = y), shape = 23, size = 2.8, fill = COL[["published"]], colour = "white", stroke = 0.6)
  }
  g
}

# accuracy against the number of features: the article's published maps (one seed each, black), and for the
# methods in C the mean over the seeds (line) with the range of the seeds (band)
curve_plot <- function(C, published, cols, labels, title = NULL, subtitle = NULL, caption = NULL, mark_n = 27) {
  C <- as.data.frame(C)
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = mark_n, colour = GRID, linewidth = 0.8) +
    ggplot2::geom_ribbon(data = C, ggplot2::aes(x = n, ymin = min, ymax = max, fill = method), alpha = 0.14) +
    ggplot2::geom_line(data = C, ggplot2::aes(x = n, y = mean, colour = method), linewidth = 0.8) +
    ggplot2::geom_line(data = published, ggplot2::aes(x = n, y = accuracy), colour = COL[["published"]], linewidth = 0.4) +
    ggplot2::geom_point(data = published, ggplot2::aes(x = n, y = accuracy, shape = "The article's published maps (seed 1245)"),
                        colour = COL[["published"]], size = 1.3) +
    ggplot2::scale_colour_manual(values = cols, labels = labels, breaks = names(cols)) +
    ggplot2::scale_fill_manual(values = cols, labels = labels, breaks = names(cols)) +
    ggplot2::scale_shape_manual(values = 16) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, ncol = 1), fill = ggplot2::guide_legend(order = 1, ncol = 1),
                    shape = ggplot2::guide_legend(order = 2)) +
    ggplot2::facet_wrap(~ set_label, nrow = 2) +
    ggplot2::labs(x = "number of features (top XGBoost gain)", y = "accuracy, all data (article's measure)", title = title, subtitle = subtitle, caption = caption) +
    theme_bench() + ggplot2::theme(legend.box = "vertical")
}

# rows with both regimes: article's schedule above, own schedule below; rows with one: on the label
dodge_y <- function(label, regime) {
  both <- tapply(regime, label, function(r) length(unique(r)) > 1)[as.character(label)]
  ifelse(both, ifelse(regime == "article", 0.18, -0.18), 0)
}

# difference from the article's SOM: mean over seeds minus kohonen's mean, with a 95% Welch interval
diff_plot <- function(Dd, title = NULL, subtitle = NULL, caption = NULL, xlab = "difference in accuracy from kohonen::som (percentage points)") {
  ord <- aggregate(Dd$diff, list(label = Dd$label), max)
  Dd$label <- factor(Dd$label, levels = ord$label[order(ord$x)])
  Dd$y <- as.integer(Dd$label) + dodge_y(Dd$label, Dd$regime)
  ggplot2::ggplot(Dd, ggplot2::aes(y = y)) +
    ggplot2::geom_vline(xintercept = 0, colour = INK2, linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = 100 * lo, xend = 100 * hi, yend = y, colour = regime), linewidth = 0.7, alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = 100 * diff, colour = regime, shape = elastic), size = 2.1, stroke = 0.9) +
    ggplot2::scale_colour_manual(values = COL[c("article", "novel")], labels = REGIME_LAB[c("article", "novel")]) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    ggplot2::scale_y_continuous(breaks = seq_along(levels(Dd$label)), labels = levels(Dd$label), expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::facet_grid(protocol ~ set_label) +
    ggplot2::labs(x = xlab, y = NULL, title = title, subtitle = subtitle, caption = caption) +
    theme_bench() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(angle = 0),
                                   panel.spacing.x = ggplot2::unit(1.2, "lines"))
}

# ------------------------------------------------------------------------------------------- the maps
# the article's Fig. 3 (aweSOM "Pie"): every unit a hexagon with a pie of the tasks it wins, the 7
# super-clusters as thick boundaries between units
hex_pie_map <- function(grid_pos, bmu, y7, super, title = NULL, subtitle = NULL, square = FALSE) {
  K <- nrow(grid_pos)
  if (square) { ang <- seq(45, 315, by = 90) * pi / 180; R <- sqrt(0.5) }       # the cells of a rectangular lattice
  else { ang <- seq(30, 330, by = 60) * pi / 180; R <- 0.5 / cos(pi / 6) }      # pointy-top hexagons, as kohonen
  nv <- length(ang)
  hex <- data.frame(unit = rep(seq_len(K), each = nv), x = rep(grid_pos[, 1], each = nv) + R * cos(ang), y = rep(grid_pos[, 2], each = nv) + R * sin(ang))
  tab <- table(factor(bmu, levels = seq_len(K)), factor(y7, levels = 1:7))
  n <- rowSums(tab)
  pies <- do.call(rbind, lapply(seq_len(K), function(u) {
    if (n[u] == 0) return(NULL)
    r <- 0.42 * sqrt(n[u] / max(n)) + 0.06
    cum <- c(0, cumsum(tab[u, ])) / n[u]
    do.call(rbind, lapply(which(tab[u, ] > 0), function(k) {
      th <- 2 * pi * seq(cum[k], cum[k + 1], length.out = max(3, ceiling(60 * (cum[k + 1] - cum[k]))))
      whole <- sum(tab[u, ] > 0) == 1                          # a whole pie is a circle: no seam to the centre
      data.frame(id = paste(u, k), task = TASKS7[k], whole = whole,
                 x = c(if (!whole) grid_pos[u, 1], grid_pos[u, 1] + r * sin(th)),
                 y = c(if (!whole) grid_pos[u, 2], grid_pos[u, 2] + r * cos(th)))
    }))
  }))
  pies$task <- factor(pies$task, levels = TASKS7)
  # boundaries: the edge two neighbouring hexagons share, where their super-clusters differ
  V <- data.frame(unit = hex$unit, x = round(hex$x, 6), y = round(hex$y, 6))
  edges <- list()
  for (a in seq_len(K - 1)) for (b in (a + 1):K) {
    if (super[a] == super[b]) next
    va <- V[V$unit == a, ]; vb <- V[V$unit == b, ]
    common <- merge(va[, c("x", "y")], vb[, c("x", "y")])
    if (nrow(common) == 2) edges[[length(edges) + 1]] <- data.frame(x = common$x[1], y = common$y[1], xend = common$x[2], yend = common$y[2])
  }
  edges <- do.call(rbind, edges)
  lab <- data.frame(x = grid_pos[, 1], y = grid_pos[, 2] - 0.47, label = ifelse(n > 0, n, ""))
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = hex, ggplot2::aes(x, y, group = unit), fill = "#f5f4f1", colour = "white", linewidth = 0.8) +
    ggplot2::geom_polygon(data = pies[pies$whole, ], ggplot2::aes(x, y, group = id, fill = task), colour = NA, show.legend = FALSE) +
    ggplot2::geom_polygon(data = pies[!pies$whole, ], ggplot2::aes(x, y, group = id, fill = task), colour = "white", linewidth = 0.25, show.legend = FALSE) +
    # the legend comes from this invisible layer alone, so it is the same on every map and patchwork can merge it
    ggplot2::geom_point(data = data.frame(task = factor(TASKS7, levels = TASKS7), x = grid_pos[1, 1], y = grid_pos[1, 2]),
                        ggplot2::aes(x, y, fill = task), shape = 22, size = 4, colour = NA, alpha = 0) +
    { if (!is.null(edges)) ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.3, colour = INK, lineend = "round") } +
    ggplot2::geom_text(data = lab, ggplot2::aes(x, y, label = label), size = 2.5, colour = INK2) +
    ggplot2::scale_fill_manual(values = TASK_COL, drop = FALSE) +
    ggplot2::coord_equal() + ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(colour = INK2),
                   plot.title = ggplot2::element_text(face = "bold", colour = INK), plot.subtitle = ggplot2::element_text(colour = INK2)) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 5))) +
    ggplot2::labs(title = title, subtitle = subtitle)
}

# ------------------------------------------------------------------------------- static figures for the README
# the headline as a static dumbbell: one row per feature set, the classical and the robust SOM, the gain beside
dumbbell_plot <- function(S2, title = NULL, subtitle = NULL, caption = NULL) {
  w <- data.table::dcast(data.table::as.data.table(S2), set_label ~ who, value.var = "mean")
  lab_c <- "Classical SOM (kohonen::som, the article)"; lab_r <- "Robust SOM (L1 winner + clipped pull)"
  ggplot2::ggplot(w, ggplot2::aes(y = set_label)) +
    ggplot2::geom_segment(ggplot2::aes(x = classic, xend = robust, yend = set_label), colour = GRID, linewidth = 3.2, lineend = "round") +
    ggplot2::geom_point(ggplot2::aes(x = classic, colour = lab_c), size = 4.4) +
    ggplot2::geom_point(ggplot2::aes(x = robust, colour = lab_r), size = 4.4) +
    ggplot2::geom_text(ggplot2::aes(x = pmax(classic, robust), label = sprintf("%+.1f", 100 * (robust - classic))),
                       hjust = -0.5, size = 3.9, fontface = "bold", colour = INK) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(COL[["kohonen"]], COL[["article"]]), c(lab_c, lab_r)), breaks = c(lab_c, lab_r)) +
    ggplot2::scale_y_discrete(limits = rev(levels(w$set_label))) +
    ggplot2::scale_x_continuous(limits = c(0.3, 1.03), breaks = seq(0.3, 1, 0.1)) +
    ggplot2::labs(x = "accuracy (share of segments named with the right task)", y = NULL, title = title, subtitle = subtitle, caption = caption) +
    theme_bench() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# child by child: the gain of every child (robust minus classical), sorted, one panel per feature set
children_plot <- function(PC, title = NULL, subtitle = NULL) {
  PC <- data.table::as.data.table(PC)[, diff := 100 * (robust - classic)]
  PC <- PC[order(set_label, -diff)][, rank := seq_len(.N), by = set_label]
  PC[, side := factor(ifelse(diff > 0, "robust SOM better", ifelse(diff < 0, "classical SOM better", "equal")),
                      levels = c("robust SOM better", "classical SOM better", "equal"))]
  ggplot2::ggplot(PC, ggplot2::aes(x = rank, y = diff, fill = side)) +
    ggplot2::geom_col(width = 0.85) + ggplot2::geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(`robust SOM better` = COL[["article"]], `classical SOM better` = COL[["kohonen"]], equal = GRID), drop = FALSE) +
    ggplot2::facet_wrap(~ set_label, nrow = 1) +
    ggplot2::labs(x = "children, sorted by the gain", y = "robust minus classical (percentage points)", title = title, subtitle = subtitle) +
    theme_bench() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank())
}
