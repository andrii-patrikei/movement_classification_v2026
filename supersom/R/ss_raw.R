# ss_raw.R: layers made from the raw IMU signals of the 648 segments, for the DTW layers of supersom.Rmd
#
# The Zenodo file (Cuberek et al. 2024, data/physical_exercise.csv, 100 Hz, back-mounted IMU) holds 90 children;
# the 9 whose gyroscope is missing are dropped (as the article's na.omit drops them), which leaves the article's
# 81 children. Every (child, activity) is one contiguous block; in the article's order (children as they first
# appear in the file, activities 1 to 8) the segments line up with the rows of the article's feature matrices.
#
# Three curve layers per segment, each a curve along an axis where an elastic (DTW) alignment means something:
#   posture   the three accelerometer axes, low-passed (three passes of a 1 s moving average, about 0.3 Hz) and
#             put on a normalised time axis of 64 points: how the trunk's orientation to gravity changes during
#             the task (lying, crawling, standing, jumping); DTW lets the phases of a task come earlier or later
#   accspec   the log power spectrum (Welch: 4 s Hann windows, half overlap) of the three accelerometer axes,
#             0.25 to 10 Hz in 40 bins; DTW along frequency lets a child step faster or slower
#   gyrospec  the same for the three gyroscope axes
# Resampling to a fixed length throws the duration away; the low-pass first keeps the resampling from aliasing
# the 2 to 3.5 Hz step rhythm into the posture curve (the spectra keep that rhythm instead).

RAW_CSV_MD5 <- "5ea4d80cbcfdaa1e54491821db57ab74"

smooth3 <- function(x, w) {                              # three passes of a centred moving average, edges reflected
  n <- length(x); h <- w %/% 2
  for (k in 1:3) {
    xp <- c(rev(x[seq_len(min(h, n))]), x, rev(x[seq.int(max(1, n - h + 1), n)]))
    f <- stats::filter(xp, rep(1 / w, w), sides = 2)
    x <- as.numeric(f[h + seq_len(n)])
  }
  x
}
on_grid <- function(x, L) stats::approx(seq(0, 1, length.out = length(x)), x, seq(0, 1, length.out = L))$y

welch_logspec <- function(x, fs = 100, win = 400, step = 200, bins = 1:40) {
  n <- length(x)
  if (n < win) x <- c(x, rep(mean(x), win - n))
  starts <- seq(1, length(x) - win + 1, by = step)
  hann <- 0.5 - 0.5 * cos(2 * pi * (0:(win - 1)) / (win - 1))
  P <- rowMeans(vapply(starts, function(s) {
    v <- x[s:(s + win - 1)]; v <- (v - mean(v)) * hann
    Mod(stats::fft(v))[bins + 1]^2 / (fs * sum(hann^2))
  }, numeric(length(bins))))
  log10(P + 1e-6)                                        # bins k = 1..40 are k * fs / win = 0.25 .. 10 Hz
}

build_raw_layers <- function(csv = "data/physical_exercise.csv", child_order, L_post = 64L) {
  if (!identical(unname(tools::md5sum(csv)), RAW_CSV_MD5)) warning(csv, " is not the Zenodo file the report was built on (md5 differs)")
  d_all <- data.table::fread(csv, select = c("accelerometer_x", "accelerometer_y", "accelerometer_z", "gyroscope_x", "gyroscope_y",
                                             "gyroscope_z", "student_id", "activity_id"), showProgress = FALSE)
  d <- d_all[!is.na(gyroscope_x)]
  # the children without a gyroscope (the article drops them): their accelerometer layers, for the missing-sensor test
  e <- d_all[is.na(gyroscope_x) & !is.na(accelerometer_x)]
  e_kids <- unique(e$student_id); e_ord <- paste(rep(e_kids, each = 8), rep(1:8, length(e_kids)))
  e_key <- paste(e$student_id, e$activity_id); e_segs <- split(seq_len(nrow(e)), factor(e_key, levels = unique(e_key)))
  stopifnot(setequal(names(e_segs), e_ord))
  EA <- as.matrix(e[, c("accelerometer_x", "accelerometer_y", "accelerometer_z")])
  ER <- t(vapply(e_segs[e_ord], function(ix) c(unlist(lapply(1:3, function(k) on_grid(smooth3(EA[ix, k], 100), L_post))),
                                               unlist(lapply(1:3, function(k) welch_logspec(EA[ix, k]))), length(ix) / 100),
                 numeric(3 * L_post + 120 + 1)))
  rm(d_all, e, EA)
  art <- unique(d$student_id)                            # the article's order of children: first appearance in the file
  stopifnot(identical(art, unique(child_order)))
  key <- paste(d$student_id, d$activity_id)
  segs <- split(seq_len(nrow(d)), factor(key, levels = unique(key)))
  ord <- paste(rep(art, each = 8), rep(1:8, length(art)))
  stopifnot(setequal(names(segs), ord))
  segs <- segs[ord]
  A <- as.matrix(d[, c("accelerometer_x", "accelerometer_y", "accelerometer_z")])
  G <- as.matrix(d[, c("gyroscope_x", "gyroscope_y", "gyroscope_z")])
  one <- function(ix) {
    a <- A[ix, , drop = FALSE]; g <- G[ix, , drop = FALSE]
    c(posture = unlist(lapply(1:3, function(k) on_grid(smooth3(a[, k], 100), L_post))),
      accspec = unlist(lapply(1:3, function(k) welch_logspec(a[, k]))),
      gyrospec = unlist(lapply(1:3, function(k) welch_logspec(g[, k]))),
      seconds = length(ix) / 100)
  }
  R <- t(vapply(segs, one, numeric(3 * L_post + 3 * 40 + 3 * 40 + 1)))
  p <- 3 * L_post
  list(posture = curve_layer(R[, seq_len(p)], L_post, 3L),
       accspec = curve_layer(R[, p + 1:120], 40L, 3L),
       gyrospec = curve_layer(R[, p + 120 + 1:120], 40L, 3L),
       seconds = R[, ncol(R)], segment = ord,
       no_gyro = list(posture = curve_layer(ER[, seq_len(p)], L_post, 3L), accspec = curve_layer(ER[, p + 1:120], 40L, 3L),
                      seconds = ER[, ncol(ER)], segment = e_ord, child = rep(e_kids, each = 8), y7 = rep(c(1:7, 2L), length(e_kids))))
}
