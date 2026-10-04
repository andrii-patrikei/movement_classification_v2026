// ss_core.cpp: the compiled core of supersom.Rmd, a multi-layer ("super") SOM on a fixed or a growing lattice
//
// Every row of X is one item with its layers side by side: layer l takes the columns off[l] .. off[l] + len[l] - 1
// and has its own distance, its own weight w[l] and its own pull. The winner minimises sum_l w[l] d_l(x, m).
//
//   mode 0 SUMSQ      sum (x - m)^2              kohonen's "sumofsquares"     m += a h (x - m)
//   mode 1 MANHATTAN  sum |x - m|                kohonen's "manhattan"        same step
//   mode 2 EUCLID     sqrt(sum (x - m)^2)        kohonen's "euclidean"        same step
//   mode 3 CORREL     1 - Pearson r                                           same step
//   mode 4 DTW        the DTW cost: the sum of squared local differences along the best path in a Sakoe-Chiba
//                     band of 'band' points; a layer of nch channels of L points (stored channel by channel)
//                     shares one alignment (dependent DTW); before the pull the item is warped onto the
//                     prototype's time axis (every prototype point takes the mean of the item points matched to it)
//
// clip[l] > 0 clips every component of layer l's pull to [-clip, clip] (the robust step of novel_SOM's l1_huber).
//
// With one SUMSQ layer, the bubble neighbourhood, the random numbers drawn inside and no growth, the loop is
// kohonen 3.0's RcppSupersom step for step: the same draw of the item, the same rule for near ties (distances
// within a relative 1e-8: an extra uniform draw), the same linear decay of alpha and radius with kohonen's 0.5
// below 1, the same update. With several SUMSQ layers it is kohonen::supersom.
//
// The growing map (supersom.Rmd, after Step 12 of som_vs_dtw_som.Rmd): the loop runs in chunks between births;
// it keeps, for every item, the unit that won it last and that distance (owner, err), and for every unit a
// running share of folds (the runner-up of its items is not a lattice neighbour); a unit with many folds pulls
// with a wider radius, r (1 + float * folds). The births themselves happen in R (ss_grow.R).
//
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cfloat>
using namespace Rcpp;

enum { SUMSQ = 0, MANHATTAN = 1, EUCLID = 2, CORREL = 3, DTW = 4 };
static const double KEPS = 1e-8;                       // kohonen's EPS (distance-functions.h)
static const double BIG = 1e300;

struct Layer { int off, len, mode, L, nch, band; double w, clip; };

static std::vector<Layer> make_layers(const IntegerVector& off, const IntegerVector& len, const IntegerVector& mode,
                                      const IntegerVector& L, const IntegerVector& nch, const IntegerVector& band,
                                      const NumericVector& w, const NumericVector& clip, int D) {
  const int nl = off.size();
  if (len.size() != nl || mode.size() != nl || L.size() != nl || nch.size() != nl || band.size() != nl || w.size() != nl || clip.size() != nl)
    stop("every layer argument needs one value per layer");
  std::vector<Layer> ly(nl);
  for (int l = 0; l < nl; l++) {
    Layer& a = ly[l];
    a.off = off[l]; a.len = len[l]; a.mode = mode[l]; a.L = L[l]; a.nch = nch[l]; a.band = band[l] < 0 ? L[l] : band[l];
    a.w = w[l]; a.clip = clip[l];
    if (a.off < 0 || a.off + a.len > D) stop("layer %d lies outside the columns of X", l + 1);
    if (a.mode < 0 || a.mode > 4) stop("unknown distance mode %d", a.mode);
    if (a.mode == DTW && a.L * a.nch != a.len) stop("DTW layer %d: L * nch must equal its number of columns", l + 1);
  }
  return ly;
}

// ------------------------------------------------------------------------------------------- distances
// R is a (L + 1) x (L + 1) table that was filled with BIG once (new_table) and is only ever written inside the
// band of this layer, so the cells outside the band stay BIG and the table can be reused without a reset
static std::vector<double> new_table(const Layer& a) {
  std::vector<double> R((size_t) (a.L + 1) * (a.L + 1), BIG);
  R[0] = 0;
  return R;
}
static double dtw_cost(const double* x, const double* m, const Layer& a, std::vector<double>& R) {
  const int L = a.L, W = L + 1, b = a.band;
  for (int i = 1; i <= L; i++) {
    int lo = std::max(1, i - b), hi = std::min(L, i + b);
    for (int j = lo; j <= hi; j++) {
      double best = R[(size_t) (i - 1) * W + (j - 1)];
      double up = R[(size_t) (i - 1) * W + j], left = R[(size_t) i * W + (j - 1)];
      if (up < best) best = up;
      if (left < best) best = left;
      double c = 0;
      for (int ch = 0; ch < a.nch; ch++) { double d = x[ch * L + i - 1] - m[ch * L + j - 1]; c += d * d; }
      R[(size_t) i * W + j] = c + best;
    }
  }
  return R[(size_t) L * W + L];
}

// walk the path of the last dtw_cost() back; every prototype point takes the mean of the item points matched to it
static void dtw_warp(const double* x, const Layer& a, const std::vector<double>& R, double* xw) {
  const int L = a.L, W = L + 1;
  std::vector<int> cnt(L, 0);
  for (int k = 0; k < L * a.nch; k++) xw[k] = 0;
  int i = L, j = L;
  while (i > 0 && j > 0) {
    for (int ch = 0; ch < a.nch; ch++) xw[ch * L + (j - 1)] += x[ch * L + (i - 1)];
    cnt[j - 1]++;
    if (i == 1 && j == 1) break;
    double d = R[(size_t) (i - 1) * W + (j - 1)], u = R[(size_t) (i - 1) * W + j], l = R[(size_t) i * W + (j - 1)];
    if (d <= u && d <= l) { i--; j--; }               // ties: the diagonal first (as ns_core.cpp)
    else if (u <= l) i--;
    else j--;
  }
  for (int j2 = 0; j2 < L; j2++) if (cnt[j2] > 0) for (int ch = 0; ch < a.nch; ch++) xw[ch * L + j2] /= (double) cnt[j2];
}

// the distance of one layer, in the arithmetic order of kohonen's distance functions
static double layer_dist(const double* x, const double* m, const Layer& a, std::vector<double>& R) {
  const int n = a.len;
  switch (a.mode) {
    case SUMSQ:     { double d = 0.0, t; for (int i = 0; i < n; ++i) { t = x[i] - m[i]; d += t * t; } return d; }
    case MANHATTAN: { double d = 0.0; for (int i = 0; i < n; ++i) d += std::fabs(x[i] - m[i]); return d; }
    case EUCLID:    { double d = 0.0, t; for (int i = 0; i < n; ++i) { t = x[i] - m[i]; d += t * t; } return std::sqrt(d); }
    case CORREL:    { double mx = 0, mm = 0; for (int i = 0; i < n; i++) { mx += x[i]; mm += m[i]; } mx /= n; mm /= n;
                      double xy = 0, xx = 0, yy = 0;
                      for (int i = 0; i < n; i++) { double p = x[i] - mx, q = m[i] - mm; xy += p * q; xx += p * p; yy += q * q; }
                      double den = std::sqrt(xx) * std::sqrt(yy);
                      return den > 0 ? 1.0 - xy / den : 1.0; }
    case DTW:       return dtw_cost(x, m, a, R);
  }
  return 0;
}

// --------------------------------------------------------------------------------------------- training
// Steps s_from .. s_to - 1 (0-based) of a run of s_total steps. Schedule as kohonen: tmp = s / s_total,
// alpha = alpha0 - (alpha0 - alpha1) tmp; the radius stays at radius0 while s < grow_steps and then falls
// linearly to radius1 over the remaining steps (grow_steps = 0: kohonen's decay); below 1 it is 0.5 (half_rule).
// pick: the items of these steps (1-based), or empty to draw them inside with R's unif_rand() as kohonen does
// (then near ties are settled by an extra draw, as kohonen; with pick given, a tie goes to the first unit).
// [[Rcpp::export]]
List ss_train_cpp(NumericMatrix X, NumericMatrix M0, NumericMatrix grid_dist,
                  IntegerVector off, IntegerVector len, IntegerVector mode, IntegerVector L, IntegerVector nch,
                  IntegerVector band, NumericVector w, NumericVector clip,
                  int s_from, int s_to, int s_total, NumericVector alpha, NumericVector radius, int grow_steps,
                  bool bubble, bool half_rule, IntegerVector pick,
                  IntegerVector owner, NumericVector err, NumericVector folds, double floatf, double memory) {
  const int N = X.nrow(), D = X.ncol(), K = M0.nrow();
  if (M0.ncol() != D) stop("M0 must have the columns of X");
  if (grid_dist.nrow() != K || grid_dist.ncol() != K) stop("grid_dist must be K x K");
  if (owner.size() != N || err.size() != N || folds.size() != K) stop("owner and err need one value per item, folds one per unit");
  const bool draw = pick.size() == 0;
  if (!draw && pick.size() != s_to - s_from) stop("pick needs one item per step");
  std::vector<Layer> ly = make_layers(off, len, mode, L, nch, band, w, clip, D);
  const int nl = ly.size();
  // row-major copies: item n at Xt[n * D], unit k at M[k * D] (kohonen's layout)
  std::vector<double> Xt((size_t) N * D), M((size_t) K * D), G((size_t) K * K);
  for (int n = 0; n < N; n++) for (int d = 0; d < D; d++) Xt[(size_t) n * D + d] = X(n, d);
  for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) M[(size_t) k * D + d] = M0(k, d);
  for (int k = 0; k < K; k++) for (int j = 0; j < K; j++) G[(size_t) k * K + j] = grid_dist(k, j);
  IntegerVector own = clone(owner); NumericVector er = clone(err); NumericVector fo = clone(folds);
  IntegerVector counts(K);
  NumericMatrix lerr(K, nl);                                    // per unit and layer: sum of w_l d_l over its wins
  std::vector<double> dist(K), dl((size_t) K * nl), xw(D);
  // one DTW table per layer and unit: filled in the winner search, read again by the update of that unit
  std::vector<std::vector<double>> tab((size_t) nl * K);
  for (int l = 0; l < nl; l++) if (ly[l].mode == DTW) for (int k = 0; k < K; k++) tab[(size_t) l * K + k] = new_table(ly[l]);
  if (draw) GetRNGstate();
  for (int s = s_from; s < s_to; s++) {
    int i = draw ? (int) (N * unif_rand()) : pick[s - s_from] - 1;
    if (i < 0 || i >= N) { if (draw) PutRNGstate(); stop("item outside X"); }
    const double* x = &Xt[(size_t) i * D];
    // 1. the winner: kohonen's FindBestMatchingUnit, over the weighted sum of the layer distances
    int nearest = -1, nind = 1; double best = DBL_MAX;
    for (int k = 0; k < K; k++) {
      double dd = 0.0;
      for (int l = 0; l < nl; l++) {
        double v = layer_dist(&x[ly[l].off], &M[(size_t) k * D + ly[l].off], ly[l], tab[(size_t) l * K + k]);
        dl[(size_t) k * nl + l] = v;
        dd += ly[l].w * v;
      }
      dist[k] = dd;
      if (dd <= best * (1 + KEPS)) {
        if (dd < best * (1 - KEPS)) { nind = 1; nearest = k; }
        else if (draw) { if (++nind * unif_rand() < 1.0) nearest = k; }
        best = dd;
      }
    }
    // the reported distance is the winner's own (kohonen keeps the last tied value; the same up to 1e-8)
    const double dwin = dist[nearest];
    // 2. schedule
    double tmp = (double) s / (double) s_total;
    double a = alpha[0] - (alpha[0] - alpha[1]) * tmp;
    double r;
    if (s < grow_steps) r = radius[0];
    else { double t2 = (double) (s - grow_steps) / (double) (s_total - grow_steps); r = radius[0] - (radius[0] - radius[1]) * t2; }
    if (half_rule && r < 1.0) r = 0.5;
    // 3. the floating radius of the winner (Step 12): a running share of folds at this unit
    if (floatf > 0 && K > 1) {
      int second = -1; double sd = DBL_MAX;
      for (int k = 0; k < K; k++) if (k != nearest && dist[k] < sd) { sd = dist[k]; second = k; }
      double fold = G[(size_t) nearest * K + second] > 1.0 + 1e-6 ? 1.0 : 0.0;
      fo[nearest] = (1 - memory) * fo[nearest] + memory * fold;
      r = r * (1 + floatf * fo[nearest]);
    }
    own[i] = nearest + 1; er[i] = dwin; counts[nearest]++;
    for (int l = 0; l < nl; l++) lerr(nearest, l) += ly[l].w * dl[(size_t) nearest * nl + l];
    // 4. the update, layer by layer
    for (int k = 0; k < K; k++) {
      double g = G[(size_t) nearest * K + k];
      double h = bubble ? (g <= r ? 1.0 : 0.0) : std::exp(-(g * g) / (2 * r * r));
      if (!(h > 0)) continue;
      double* m = &M[(size_t) k * D];
      for (int l = 0; l < nl; l++) {
        const Layer& ay = ly[l];
        const double* xl = &x[ay.off];
        double* ml = &m[ay.off];
        if (ay.mode == DTW) { dtw_warp(xl, ay, tab[(size_t) l * K + k], &xw[ay.off]); xl = &xw[ay.off]; }   // the table of the search
        if (ay.clip > 0) {
          for (int j = 0; j < ay.len; j++) {
            double dlt = xl[j] - ml[j];
            if (dlt > ay.clip) dlt = ay.clip; else if (dlt < -ay.clip) dlt = -ay.clip;
            ml[j] += h * a * dlt;
          }
        } else {
          for (int j = 0; j < ay.len; j++) ml[j] += h * a * (xl[j] - ml[j]);
        }
      }
    }
  }
  if (draw) PutRNGstate();
  NumericMatrix Mout(K, D);
  for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) Mout(k, d) = M[(size_t) k * D + d];
  return List::create(_["M"] = Mout, _["owner"] = own, _["err"] = er, _["folds"] = fo, _["counts"] = counts, _["layer_err"] = lerr);
}

// -------------------------------------------------------------------------------- distances for judging
// every row of X against every row of M: the weighted sum over the layers (per_layer = FALSE, an N x K matrix),
// or a list of the unweighted N x K matrices of every layer (per_layer = TRUE)
// [[Rcpp::export]]
SEXP ss_cross_cpp(NumericMatrix X, NumericMatrix M, IntegerVector off, IntegerVector len, IntegerVector mode,
                  IntegerVector L, IntegerVector nch, IntegerVector band, NumericVector w, bool per_layer = false) {
  const int N = X.nrow(), K = M.nrow(), D = X.ncol();
  if (M.ncol() != D) stop("X and M need the same columns");
  NumericVector clip(off.size());
  std::vector<Layer> ly = make_layers(off, len, mode, L, nch, band, w, clip, D);
  const int nl = ly.size();
  std::vector<double> x(D), m((size_t) K * D);
  std::vector<std::vector<double>> tab(nl);
  for (int l = 0; l < nl; l++) if (ly[l].mode == DTW) tab[l] = new_table(ly[l]);
  for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) m[(size_t) k * D + d] = M(k, d);
  std::vector<NumericMatrix> per; for (int l = 0; l < nl; l++) per.push_back(NumericMatrix(per_layer ? N : 0, per_layer ? K : 0));
  NumericMatrix out(per_layer ? 0 : N, per_layer ? 0 : K);
  for (int n = 0; n < N; n++) {
    for (int d = 0; d < D; d++) x[d] = X(n, d);
    for (int k = 0; k < K; k++) {
      double dd = 0.0;
      for (int l = 0; l < nl; l++) {
        double v = layer_dist(&x[ly[l].off], &m[(size_t) k * D + ly[l].off], ly[l], tab[l]);
        if (per_layer) per[l](n, k) = v; else dd += ly[l].w * v;
      }
      if (!per_layer) out(n, k) = dd;
    }
  }
  if (!per_layer) return out;
  List res(nl); for (int l = 0; l < nl; l++) res[l] = per[l];
  return res;
}

// the item warped onto a prototype's time axis in one DTW layer (for pictures and checks)
// [[Rcpp::export]]
NumericVector ss_warp_cpp(NumericVector x, NumericVector m, int L, int nch, int band) {
  Layer a; a.off = 0; a.len = L * nch; a.mode = DTW; a.L = L; a.nch = nch; a.band = band < 0 ? L : band; a.w = 1; a.clip = 0;
  std::vector<double> xx(x.begin(), x.end()), mm(m.begin(), m.end()), R = new_table(a), xw(a.len);
  dtw_cost(xx.data(), mm.data(), a, R);
  dtw_warp(xx.data(), a, R, xw.data());
  return NumericVector(xw.begin(), xw.end());
}

// the distance of row i of A to row i of B, every layer apart (unweighted): an N x nlayers matrix
// [[Rcpp::export]]
NumericMatrix ss_rowdist_cpp(NumericMatrix A, NumericMatrix B, IntegerVector off, IntegerVector len, IntegerVector mode,
                             IntegerVector L, IntegerVector nch, IntegerVector band) {
  const int N = A.nrow(), D = A.ncol();
  if (B.nrow() != N || B.ncol() != D) stop("A and B must have the same shape");
  NumericVector w(off.size(), 1.0), clip(off.size());
  std::vector<Layer> ly = make_layers(off, len, mode, L, nch, band, w, clip, D);
  const int nl = ly.size();
  std::vector<std::vector<double>> tab(nl);
  for (int l = 0; l < nl; l++) if (ly[l].mode == DTW) tab[l] = new_table(ly[l]);
  std::vector<double> a(D), b(D);
  NumericMatrix out(N, nl);
  for (int n = 0; n < N; n++) {
    for (int d = 0; d < D; d++) { a[d] = A(n, d); b[d] = B(n, d); }
    for (int l = 0; l < nl; l++) out(n, l) = layer_dist(&a[ly[l].off], &b[ly[l].off], ly[l], tab[l]);
  }
  return out;
}

// ----------------------------------------------------------------------------------- topographic product
// Bauer & Pawelzik (1992), as ns_core.cpp of novel_SOM: dV distances between prototypes in the data space,
// dA distances on the map; ties on the map are broken by the data distance
// [[Rcpp::export]]
double topographic_product_cpp(NumericMatrix dV, NumericMatrix dA) {
  const int K = dV.nrow();
  if (K < 3) return 0;
  double total = 0;
  std::vector<int> idx(K - 1);
  for (int j = 0; j < K; j++) {
    int c = 0;
    for (int k = 0; k < K; k++) if (k != j) idx[c++] = k;
    std::vector<int> oV = idx, oA = idx;
    std::sort(oV.begin(), oV.end(), [&](int a, int b) {
      if (dV(j, a) != dV(j, b)) return dV(j, a) < dV(j, b);
      if (dA(j, a) != dA(j, b)) return dA(j, a) < dA(j, b);
      return a < b; });
    std::sort(oA.begin(), oA.end(), [&](int a, int b) {
      if (dA(j, a) != dA(j, b)) return dA(j, a) < dA(j, b);
      if (dV(j, a) != dV(j, b)) return dV(j, a) < dV(j, b);
      return a < b; });
    double logprod = 0;
    for (int k = 1; k <= K - 1; k++) {
      double q1 = dV(j, oA[k - 1]) / dV(j, oV[k - 1]);
      double q2 = dA(j, oA[k - 1]) / dA(j, oV[k - 1]);
      if (!(q1 > 0) || !(q2 > 0) || !std::isfinite(q1) || !std::isfinite(q2)) { q1 = q2 = 1; }
      logprod += std::log(q1) + std::log(q2);
      total += logprod / (2.0 * k);
    }
  }
  return total / ((double) K * (K - 1));
}
