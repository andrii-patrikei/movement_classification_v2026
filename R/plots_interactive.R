# plots_interactive.R: the interactive figures and tables of som_variants.Rmd (plotly and DT)
#
# The same roles as the static figures of plots.R: blue the article's SOM (kohonen), orange the robust SOM,
# aqua novel_SOM's own schedule, black the article's published maps. Every chart has a hover layer; the
# buttons and menus switch what is shown without changing the colours.

P_FONT <- list(family = "Lato, 'Helvetica Neue', Helvetica, Arial, sans-serif", size = 13, color = INK)
COL_ROBUST <- COL[["article"]]; COL_CLASSIC <- COL[["kohonen"]]; COL_OWN <- COL[["novel"]]; COL_PUB <- COL[["published"]]
rgba <- function(hex, a) { v <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", v[1], v[2], v[3], a) }

p_axis <- function(title = NULL, ...) utils::modifyList(list(title = list(text = title, font = list(size = 13, color = INK2)), gridcolor = GRID,
                                                              zeroline = FALSE, linecolor = GRID, tickfont = list(color = INK2, size = 12)), list(...))
p_finish <- function(p, ..., legend_y = -0.18, legend = NULL) {
  if (is.null(legend)) legend <- list(orientation = "h", x = 0, y = legend_y, xanchor = "left", yanchor = "top", font = list(color = INK2))
  p <- plotly::layout(p, font = P_FONT, paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
                      hoverlabel = list(bgcolor = "white", bordercolor = GRID, font = list(color = INK, size = 13)),
                      legend = legend, ...)
  p <- plotly::config(p, displaylogo = FALSE, toImageButtonOptions = list(format = "png", scale = 2),
                      modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "toggleSpikelines", "hoverCompareCartesian"))
  # as wide as the text column on any screen (knitr would otherwise give the widget the figure width in pixels)
  p$sizingPolicy$knitr$figure <- FALSE
  p$sizingPolicy$defaultWidth <- "100%"
  p
}
p_buttons <- function(labels, visible, extra = NULL, x = 0, y = 1.22, type = "buttons", active = 0) {
  list(list(type = type, direction = "right", x = x, y = y, xanchor = "left", yanchor = "top", showactive = TRUE, active = active,
            bgcolor = "white", bordercolor = GRID, font = list(color = INK, size = 12), pad = list(r = 4, t = 2, b = 2),
            buttons = lapply(seq_along(labels), function(i)
              list(label = labels[i], method = "update",
                   args = c(list(list(visible = visible[[i]])), if (!is.null(extra)) list(extra[[i]]))))))
}
pm <- function(x, d = 3) formatC(x, digits = d, format = "f")

# ------------------------------------------------------------------------------------------ the headline
# S2: set_label, protocol, who ("classic"/"robust"), mean, sd, min, max, seeds; one row per feature set,
# protocol and map. A dumbbell per feature set: the article's map, the robust SOM, the gain written beside it.
p_dumbbell <- function(S2, protocols, labels) {
  sets <- levels(S2$set_label); yrev <- rev(sets)
  p <- plotly::plot_ly(height = 470)
  vis <- list(); ann <- list(); k <- 0
  for (pi in seq_along(protocols)) {
    d <- S2[S2$protocol == protocols[pi], ]
    cl <- d[d$who == "classic", ][match(sets, d$set_label[d$who == "classic"]), ]
    rb <- d[d$who == "robust", ][match(sets, d$set_label[d$who == "robust"]), ]
    tip <- function(x, who) sprintf("<b>%s</b><br>%s<br>accuracy %s &plusmn; %s<br>seeds %s to %s (%d seeds)", x$set_label, who,
                                    pm(x$mean), pm(x$sd), pm(x$min), pm(x$max), x$seeds)
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = as.vector(rbind(cl$mean, rb$mean, NA)),
                           y = as.vector(rbind(sets, sets, NA)), line = list(color = GRID, width = 5), hoverinfo = "skip",
                           showlegend = FALSE, visible = pi == 1)
    p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = cl$mean, y = sets, name = "Classical SOM (the article's kohonen::som)",
                           marker = list(color = COL_CLASSIC, size = 14, line = list(color = "white", width = 2)),
                           text = tip(cl, "Classical SOM (kohonen::som)"), hovertemplate = "%{text}<extra></extra>",
                           legendgroup = "classic", visible = pi == 1)
    p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = rb$mean, y = sets, name = "Robust SOM (L1 winner + clipped pull)",
                           marker = list(color = COL_ROBUST, size = 14, line = list(color = "white", width = 2)),
                           text = tip(rb, "Robust SOM"), hovertemplate = "%{text}<extra></extra>",
                           legendgroup = "robust", visible = pi == 1)
    vis[[pi]] <- rep(seq_along(protocols) == pi, each = 3)
    ann[[pi]] <- lapply(seq_along(sets), function(i) list(x = max(cl$mean[i], rb$mean[i]), y = sets[i], xref = "x", yref = "y",
                                                         text = sprintf("<b>%+.1f</b>", 100 * (rb$mean[i] - cl$mean[i])),
                                                         xanchor = "left", xshift = 14, showarrow = FALSE, font = list(color = INK, size = 13)))
  }
  p_finish(p, xaxis = p_axis("accuracy (share of segments named with the right task)", range = c(0.3, 1.04), dtick = 0.1),
           yaxis = p_axis(NULL, type = "category", categoryorder = "array", categoryarray = yrev, showgrid = FALSE),
           annotations = ann[[1]], updatemenus = p_buttons(labels, vis, lapply(ann, function(a) list(annotations = a)), y = 1.16),
           margin = list(t = 60, l = 10, r = 40, b = 120), legend_y = -0.22)
}

# ------------------------------------------------------------------------------------- seed by seed
# Rs: set_label, protocol, who, seed, value; pub: set_label, value (the article's published map, all data only)
p_seed_boxes <- function(Rs, pub, protocols, labels) {
  p <- plotly::plot_ly(height = 500); vis <- list(); n_tr <- 0
  for (pi in seq_along(protocols)) {
    for (w in c("classic", "robust")) {
      d <- Rs[Rs$protocol == protocols[pi] & Rs$who == w, ]
      col <- if (w == "classic") COL_CLASSIC else COL_ROBUST
      p <- plotly::add_trace(p, type = "box", x = d$set_label, y = d$value, name = if (w == "classic") "Classical SOM" else "Robust SOM",
                             boxpoints = "all", jitter = 0.45, pointpos = 0, fillcolor = rgba(col, 0.12), line = list(color = col, width = 1.4),
                             marker = list(color = col, size = 5, opacity = 0.75), legendgroup = w,
                             text = sprintf("seed %d", d$seed), hovertemplate = "%{x}<br>%{text}: %{y:.3f}<extra></extra>",
                             offsetgroup = w, visible = pi == 1)
    }
    n_tr <- n_tr + 2
  }
  p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = pub$set_label, y = pub$value, name = "The article's published map (seed 1245)",
                         marker = list(symbol = "diamond", color = COL_PUB, size = 11, line = list(color = "white", width = 1.5)),
                         hovertemplate = "%{x}<br>the article's published map: %{y:.3f}<extra></extra>", visible = protocols[1] == "accuracy")
  vis <- lapply(seq_along(protocols), function(pi) c(rep(seq_along(protocols) == pi, each = 2), protocols[pi] == "accuracy"))
  p_finish(p, boxmode = "group",
           xaxis = p_axis(NULL, type = "category", categoryorder = "array", categoryarray = levels(Rs$set_label)),
           yaxis = p_axis("accuracy"), updatemenus = p_buttons(labels, vis, y = 1.14), margin = list(t = 60, l = 60, r = 10, b = 100),
           legend_y = -0.12)
}

# --------------------------------------------------------------------------------- beyond the seed noise
# dd2: set_label, protocol_label, diff, lo, hi (robust minus classical, accuracy units)
p_forest <- function(dd2) {
  sets <- levels(dd2$set_label); protos <- levels(dd2$protocol_label)
  off <- setNames(c(0.15, -0.15)[seq_along(protos)], protos)
  sym <- setNames(c("circle", "circle-open")[seq_along(protos)], protos)
  p <- plotly::plot_ly(height = 430)
  for (pr in protos) {
    d <- dd2[dd2$protocol_label == pr, ]; y <- match(as.character(d$set_label), rev(sets)) + off[[pr]]
    p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = 100 * d$diff, y = y, name = pr,
                           error_x = list(type = "data", symmetric = FALSE, array = 100 * (d$hi - d$diff), arrayminus = 100 * (d$diff - d$lo),
                                          color = COL_ROBUST, thickness = 1.6, width = 0),
                           marker = list(color = COL_ROBUST, symbol = sym[[pr]], size = 11, line = list(color = COL_ROBUST, width = 2)),
                           text = sprintf("<b>%s</b>, %s<br>%+.1f points, 95%% interval %+.1f to %+.1f", d$set_label, pr, 100 * d$diff, 100 * d$lo, 100 * d$hi),
                           hovertemplate = "%{text}<extra></extra>")
  }
  p_finish(p, xaxis = p_axis("robust SOM minus classical SOM (percentage points of accuracy)", zeroline = TRUE, zerolinecolor = INK2, zerolinewidth = 1.2),
           yaxis = p_axis(NULL, tickvals = seq_along(sets), ticktext = rev(sets), range = c(0.4, length(sets) + 0.6), showgrid = FALSE),
           margin = list(t = 20, l = 10, r = 10, b = 120), legend_y = -0.22)
}

# ----------------------------------------------------------------------------------- child by child
# PC: set_label, child, classic, robust (per-child accuracy on held-out children, mean over the seeds);
# one bar per child, sorted by the gain; a menu switches the feature set
p_children <- function(PC, first = levels(PC$set_label)[1]) {
  sets <- levels(PC$set_label); p <- plotly::plot_ly(height = 430); ann <- list()
  for (s in sets) {
    d <- PC[PC$set_label == s, ]; d$diff <- 100 * (d$robust - d$classic); d <- d[order(-d$diff, d$child), ]
    col <- ifelse(d$diff > 0, COL_ROBUST, ifelse(d$diff < 0, COL_CLASSIC, GRID))
    p <- plotly::add_trace(p, type = "bar", x = seq_len(nrow(d)), y = d$diff, marker = list(color = col), showlegend = FALSE,
                           text = sprintf("<b>child %s</b><br>classical SOM %s<br>robust SOM %s<br>%+.1f points", d$child, pm(d$classic),
                                          pm(d$robust), d$diff), textposition = "none", hovertemplate = "%{text}<extra></extra>", visible = s == first)
    nb <- sum(d$diff > 0); nw <- sum(d$diff < 0); ne <- sum(d$diff == 0)
    ann[[s]] <- list(list(x = 0, y = 1.02, xref = "paper", yref = "paper", xanchor = "left", yanchor = "bottom", showarrow = FALSE,
                          text = sprintf("<b>%s</b>: robust SOM better for <b><span style='color:%s'>%d</span></b> children, classical SOM better for <b><span style='color:%s'>%d</span></b>, equal for %d",
                                         s, COL_ROBUST, nb, COL_CLASSIC, nw, ne), font = list(size = 13, color = INK)))
  }
  vis <- lapply(sets, function(s) sets == s)
  p_finish(p, bargap = 0.15, xaxis = p_axis("children, sorted by the gain", showticklabels = FALSE, showgrid = FALSE),
           yaxis = p_axis("robust minus classical (percentage points)", zeroline = TRUE, zerolinecolor = INK2),
           annotations = ann[[first]], margin = list(t = 95, l = 60, r = 10, b = 40),
           updatemenus = p_buttons(sets, vis, lapply(ann, function(a) list(annotations = a)), y = 1.27, active = match(first, sets) - 1))
}

# ------------------------------------------------------------------------------------ the leaderboard
# L: label, regime_label, and for every measure a matrix of set means; one heatmap, a menu switches the measure
p_leaderboard <- function(mats, row_labels, col_labels, measure_labels, lower_better, bold = integer(0)) {
  row_labels[bold] <- sprintf("<b>%s</b>", row_labels[bold])
  blues <- list(c(0, "#f0f6fe"), c(0.5, "#6da7ec"), c(1, "#0d366b"))
  first <- mats[[1]]
  txt <- function(m) matrix(ifelse(is.na(m), "", formatC(m, digits = 2, format = "f")), nrow(m))
  p <- plotly::plot_ly(height = 26 * length(row_labels) + 140) |>
    plotly::add_trace(type = "heatmap", z = first, x = col_labels, y = row_labels, text = txt(first), texttemplate = "%{text}",
                      textfont = list(size = 10), colorscale = blues, reversescale = lower_better[1], xgap = 2, ygap = 2,
                      hovertemplate = "%{y}<br>%{x}: <b>%{z:.3f}</b><extra></extra>", hoverongaps = FALSE,
                      colorbar = list(thickness = 10, outlinewidth = 0, tickfont = list(color = INK2)))
  buttons <- lapply(seq_along(mats), function(i)
    list(label = measure_labels[i], method = "restyle",
         args = list(list(z = list(mats[[i]]), text = list(txt(mats[[i]])), reversescale = lower_better[i]))))
  p_finish(p, xaxis = p_axis(NULL, side = "top", showgrid = FALSE, tickangle = 0, tickfont = list(color = INK2, size = 11)),
           yaxis = p_axis(NULL, autorange = "reversed", showgrid = FALSE),
           margin = list(t = 120, l = 10, r = 10, b = 10),
           updatemenus = list(list(type = "dropdown", x = 1, y = 1, yref = "paper", xanchor = "right", yanchor = "bottom", pad = list(b = 40),
                                   bgcolor = "white", bordercolor = GRID, font = list(color = INK, size = 12), buttons = buttons)))
}

# -------------------------------------------------------------------------------------- the curve
# C: set_label, n, who, mean, min, max; pub: set_label, n, accuracy; one feature set at a time
p_curve <- function(C, pub, who_cols, who_labels, first = levels(C$set_label)[1]) {
  sets <- levels(C$set_label); p <- plotly::plot_ly(height = 480); per <- 0
  for (s in sets) {
    for (w in names(who_cols)) {
      d <- C[C$set_label == s & C$who == w, ]; d <- d[order(d$n), ]
      p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = c(d$n, rev(d$n)), y = c(d$max, rev(d$min)), fill = "toself",
                             fillcolor = rgba(who_cols[[w]], 0.13), line = list(width = 0), hoverinfo = "skip", showlegend = FALSE,
                             legendgroup = w, visible = s == first)
      p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = d$n, y = d$mean, name = who_labels[[w]], line = list(color = who_cols[[w]], width = 2.5),
                             text = sprintf("range of the seeds %s to %s", pm(d$min), pm(d$max)),
                             hovertemplate = paste0(who_labels[[w]], ": <b>%{y:.3f}</b> (%{text})<extra></extra>"),
                             legendgroup = w, visible = s == first)
    }
    d <- pub[pub$set_label == s, ]; d <- d[order(d$n), ]
    p <- plotly::add_trace(p, type = "scatter", mode = "lines+markers", x = d$n, y = d$accuracy, name = "The article's published maps (one seed each)",
                           line = list(color = COL_PUB, width = 1), marker = list(color = COL_PUB, size = 6),
                           hovertemplate = "the article's published map: <b>%{y:.3f}</b><extra></extra>",
                           legendgroup = "pub", visible = s == first)
    per <- 2 * length(who_cols) + 1
  }
  vis <- lapply(seq_along(sets), function(i) rep(seq_along(sets) == i, each = per))
  p_finish(p, hovermode = "x unified",
           xaxis = p_axis("number of features (the top features by XGBoost gain)", dtick = 5, range = c(1.5, 35.5)),
           yaxis = p_axis("accuracy, all 648 segments"),
           shapes = list(list(type = "line", x0 = 27, x1 = 27, y0 = 0, y1 = 1, yref = "paper", line = list(color = GRID, width = 3))),
           annotations = list(list(x = 27, y = 0, yref = "paper", text = "27 features (the article)", showarrow = FALSE, xanchor = "left",
                                   yanchor = "bottom", xshift = 6, font = list(color = INK2, size = 11))),
           updatemenus = p_buttons(sets, vis, y = 1.12, active = match(first, sets) - 1), margin = list(t = 60, l = 60, r = 10, b = 130),
           legend = list(orientation = "h", x = 0, y = -0.22, xanchor = "left", yanchor = "top", font = list(color = INK2)))
}

# ----------------------------------------------------------------------------------------------- the maps
# The article's Fig. 2 (aweSOM "Pie"): every unit a cell with a pie of the tasks it wins (area by the number of
# segments), the 7 super-clusters as thick boundaries. All maps share one coordinate system, side by side, so
# the cells keep their shape; hovering a cell lists its segments. maps: a list of list(pos, bmu, super, acc,
# title, subtitle, square), acc being the result of article_accuracy(); y7: the task of every segment.
p_maps <- function(maps, y7, gap = 1.6) {
  p <- plotly::plot_ly(height = 430)
  ann <- list(); cells <- list(); wedges <- vector("list", 7); edges <- list(); hov <- list(); counts <- list()
  x_shift <- 0; nmax <- max(sapply(maps, function(m) max(table(m$bmu))))
  for (mi in seq_along(maps)) {
    m <- maps[[mi]]; pos <- m$pos; pos[, 1] <- pos[, 1] - min(pos[, 1]) + x_shift; K <- nrow(pos)
    if (isTRUE(m$square)) { ang <- seq(45, 315, by = 90) * pi / 180; R <- sqrt(0.5) }
    else { ang <- seq(30, 330, by = 60) * pi / 180; R <- 0.5 / cos(pi / 6) }
    hx <- lapply(seq_len(K), function(u) cbind(pos[u, 1] + R * cos(ang), pos[u, 2] + R * sin(ang)))
    cells[[mi]] <- do.call(rbind, lapply(hx, function(h) rbind(h, h[1, ], c(NA, NA))))
    tab <- table(factor(m$bmu, levels = seq_len(K)), factor(y7, levels = 1:7)); n <- rowSums(tab)
    for (u in which(n > 0)) {
      r <- 0.42 * sqrt(n[u] / nmax) + 0.06; cum <- c(0, cumsum(tab[u, ])) / n[u]
      whole <- sum(tab[u, ] > 0) == 1                              # a whole pie is a circle: no seam to the centre
      for (k in which(tab[u, ] > 0)) {
        th <- 2 * pi * seq(cum[k], cum[k + 1], length.out = max(3, ceiling(60 * (cum[k + 1] - cum[k]))))
        xy <- cbind(c(if (!whole) pos[u, 1], pos[u, 1] + r * sin(th)), c(if (!whole) pos[u, 2], pos[u, 2] + r * cos(th)))
        wedges[[k]] <- rbind(wedges[[k]], xy, xy[1, ], c(NA, NA))
      }
    }
    for (a in seq_len(K - 1)) for (b in (a + 1):K) if (m$super[a] != m$super[b]) {
      common <- merge(data.frame(x = round(hx[[a]][, 1], 6), y = round(hx[[a]][, 2], 6)),
                      data.frame(x = round(hx[[b]][, 1], 6), y = round(hx[[b]][, 2], 6)))
      if (nrow(common) == 2) edges[[length(edges) + 1]] <- rbind(as.matrix(common), c(NA, NA))
    }
    hov[[mi]] <- data.frame(x = pos[, 1], y = pos[, 2], text = sapply(seq_len(K), function(u) {
      parts <- if (n[u] > 0) paste(sprintf("%s: %d", TASKS7[tab[u, ] > 0], tab[u, tab[u, ] > 0]), collapse = "<br>") else "no segments"
      sprintf("<b>%s</b><br>unit %d, super-cluster %d (named %s)<br>%d segments<br>%s", m$title, u, m$super[u],
              TASKS7[m$acc$task[m$super[u]]], n[u], parts) }))
    counts[[mi]] <- data.frame(x = pos[, 1], y = pos[, 2] - 0.5, label = ifelse(n > 0, n, ""))
    ann[[mi]] <- list(x = mean(range(pos[, 1])), y = max(pos[, 2]) + 0.95, xref = "x", yref = "y", showarrow = FALSE, xanchor = "center",
                      text = sprintf("<b>%s</b><br><span style='color:%s'>%s</span>", m$title, INK2, m$subtitle), font = list(size = 13, color = INK))
    x_shift <- max(pos[, 1]) + gap + 0.5
  }
  cl <- do.call(rbind, cells)
  p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = cl[, 1], y = cl[, 2], fill = "toself", fillcolor = "#f5f4f1",
                         line = list(color = "white", width = 1.5), hoverinfo = "skip", showlegend = FALSE)
  for (k in 1:7) if (!is.null(wedges[[k]]))
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = wedges[[k]][, 1], y = wedges[[k]][, 2], fill = "toself",
                           fillcolor = TASK_COL[[k]], line = list(color = "white", width = 0.6), name = TASKS7[[k]], hoverinfo = "skip")
  if (length(edges)) {
    e <- do.call(rbind, edges)
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = e[, 1], y = e[, 2], line = list(color = INK, width = 3),
                           hoverinfo = "skip", showlegend = FALSE)
  }
  h <- do.call(rbind, hov); cn <- do.call(rbind, counts)
  p <- plotly::add_trace(p, type = "scatter", mode = "text", x = cn$x, y = cn$y, text = cn$label, textfont = list(size = 9, color = INK2),
                         hoverinfo = "skip", showlegend = FALSE)
  p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = h$x, y = h$y, text = h$text, marker = list(size = 34, color = "rgba(0,0,0,0)"),
                         hovertemplate = "%{text}<extra></extra>", showlegend = FALSE)
  # fixed ranges with constrain = "domain": plotly shrinks the drawing area to keep the cells' shape instead of
  # cutting the range, so every map stays in view at any width
  p_finish(p, xaxis = list(visible = FALSE, fixedrange = TRUE, range = range(cl[, 1], na.rm = TRUE) + c(-0.3, 0.3), constrain = "domain"),
           yaxis = list(visible = FALSE, fixedrange = TRUE, range = range(cl[, 2], na.rm = TRUE) + c(-0.75, 1.45), scaleanchor = "x",
                        scaleratio = 1, constrain = "domain"),
           annotations = ann, margin = list(t = 10, l = 0, r = 0, b = 10),
           legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.02, yanchor = "top", font = list(color = INK2, size = 12)))
}

# -------------------------------------------------------------------------------------------- tables
dt_table <- function(df, digits = 3, round_cols = NULL, page = 10, filter = "top", order = NULL, ...) {
  opts <- list(pageLength = page, dom = "Bfrtip", buttons = c("copy", "csv"), autoWidth = FALSE, scrollX = TRUE,
               lengthMenu = c(10, 25, 50, 100))
  if (!is.null(order)) opts$order <- order
  w <- DT::datatable(df, rownames = FALSE, filter = filter, extensions = "Buttons", options = opts, class = "compact stripe hover", ...)
  if (length(round_cols)) w <- DT::formatRound(w, round_cols, digits)
  w
}
