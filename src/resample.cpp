// Observed and resampled pointwise statistics (the "front end").
//
// One-sample sign flips and paired swaps act on "units": single observations
// (fdatest's scheme) or blocks of observations whose signs move together.
// Unit-level 0/1 indicators are drawn as stats::rbinom(q, 1, 0.5) per resample,
// so with units = observations the RNG stream, and every resampled statistic,
// match fdatest's os_single_perm() / ts_single_perm(). Unpaired permutations
// replicate sample(n) via R_unif_index().
//
// Ported from fdatest (Pini, Stamm, Vantini; GPL >= 3), development version
// commit 67bb41c: R/utils-test-os.R and R/utils-test-ts.R.

#include "common.h"
#include <R_ext/Random.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

namespace {

// fdatest's ts_pointwise_stat(): alt 0 = two.sided, 1 = greater, 2 = less.
//   sign_diff <- sign(meandiff); sign_diff[sign_diff == -1] <- 0
//   greater: (meandiff * sign_diff)^2,  less: (meandiff * (sign_diff - 1))^2
inline double alt_stat(double d, int alt) {
  if (alt == 0) return d * d;
  const double sd = ISNAN(d) ? d : (d > 0 ? 1.0 : 0.0);
  const double a = (alt == 1) ? d * sd : d * (sd - 1.0);
  return a * a;
}

// diag(stats::cov(x)) as R computes it (use = "everything", C function
// cov_na_1): two-pass long double mean rounded to double, then the long double
// sum of squared deviations divided by n - 1. The steps follow R's
// src/library/stats/src/cov.c (Copyright The R Core Team, GPL >= 2).
inline double r_cov_var(const double* v, int n) {
  if (n <= 1) return NA_REAL;
  ldouble sum = 0.0L;
  for (int i = 0; i < n; ++i) sum += v[i];
  ldouble tmp = sum / n;
  if (R_FINITE((double)tmp)) {
    sum = 0.0L;
    for (int i = 0; i < n; ++i) sum += (v[i] - tmp);
    tmp = tmp + sum / n;
  }
  const ldouble xm = (double)tmp;
  sum = 0.0L;
  for (int i = 0; i < n; ++i) sum += (v[i] - xm) * (v[i] - xm);
  return (double)(sum / (R_xlen_t)(n - 1));
}

// 0/1 indicators of resample b, one per unit: drawn as rbinom(q, 1, 0.5)
// (one uniform per unit, in unit order) or read from a user-supplied matrix.
// For size 1 and prob 0.5, R's rbinom() (inversion branch of nmath/rbinom.c)
// consumes exactly one unif_rand() u and returns u >= 0.5; calling unif_rand()
// directly gives the same draws and stream at a fifth of the cost.
inline void unit_indicators(int b, int q, bool draw, const IntegerMatrix& ind,
                            std::vector<int>& f) {
  if (draw) {
    for (int g = 0; g < q; ++g) f[g] = unif_rand() >= 0.5;
  } else {
    for (int g = 0; g < q; ++g) f[g] = ind(b, g);
  }
}

// Unit sums of the columns of x (n x p): U[k * q + g] = sum over i in unit g of
// x[i, k], accumulated in long double in observation order. For singleton units
// this is x itself, so the statistics below reproduce colMeans() of the
// sign-flipped or swapped data exactly.
std::vector<double> unit_sums(const NumericMatrix& x, const IntegerVector& unit, int q) {
  const int n = x.nrow(), p = x.ncol();
  std::vector<double> U((size_t)q * p);
  std::vector<ldouble> acc(q);
  for (int k = 0; k < p; ++k) {
    std::fill(acc.begin(), acc.end(), 0.0L);
    const double* col = x.begin() + (R_xlen_t)k * n;
    for (int i = 0; i < n; ++i) acc[unit[i]] += col[i];
    for (int g = 0; g < q; ++g) U[(size_t)k * q + g] = (double)acc[g];
  }
  return U;
}

}  // namespace

// One-sample sign flips. x = data - mu (n x p). Indicator 1 keeps the sign
// (fdatest: signs <- rbinom(n, 1, 0.5) * 2 - 1). Statistic: alt_stat of the mean
// of the flipped curves y. With standardize = TRUE it is divided by var() of y,
// computed curve by curve on every resample exactly as R's
// colMeans(y)^2 / apply(y, 2, var) (a two-pass variance: a one-pass formula
// cancels catastrophically where the curves are nearly equal).
// [[Rcpp::export]]
List cpp_stats_onesample(NumericMatrix x, IntegerVector unit, int q, int B,
                         IntegerMatrix ind, int alt, bool standardize) {
  const int n = x.nrow(), p = x.ncol();
  const bool draw = ind.nrow() == 0;
  NumericVector t0(p);
  NumericMatrix tp(B, p);
  std::vector<int> f(q, 1);

  if (standardize) {
    std::vector<double> sg(n), v(n);
    auto stats_std = [&](double* out, R_xlen_t stride) {
      for (int i = 0; i < n; ++i) sg[i] = f[unit[i]] ? 1.0 : -1.0;
      for (int k = 0; k < p; ++k) {
        const double* col = x.begin() + (R_xlen_t)k * n;
        ldouble s = 0.0L;
        for (int i = 0; i < n; ++i) {
          v[i] = col[i] * sg[i];
          s += v[i];
        }
        out[k * stride] = alt_stat((double)(s / n), alt) / r_cov_var(v.data(), n);
      }
    };
    stats_std(t0.begin(), 1);
    for (int b = 0; b < B; ++b) {
      if ((b & 255) == 0) Rcpp::checkUserInterrupt();
      unit_indicators(b, q, draw, ind, f);
      stats_std(tp.begin() + b, B);
    }
    return List::create(_["t0"] = t0, _["t_perm"] = tp);
  }

  const std::vector<double> U = unit_sums(x, unit, q);
  // interleaved: W[2 (k q + g)] = -U, W[2 (k q + g) + 1] = +U, so that the sign
  // is an index offset and the long double add reads its operand from memory
  // (negation is exact: these are the values of fdatest's coeff * signs)
  std::vector<double> W(2 * (size_t)q * p);
  for (size_t c = 0; c < (size_t)q * p; ++c) {
    W[2 * c] = -U[c];
    W[2 * c + 1] = U[c];
  }
  std::vector<int> idx(q);
  auto finish = [&](ldouble s, int) { return alt_stat((double)(s / n), alt); };
  // four columns in lockstep (independent long double chains, same arithmetic)
  auto stats = [&](double* out, R_xlen_t stride) {
    for (int g = 0; g < q; ++g) idx[g] = 2 * g + f[g];
    int k = 0;
    for (; k + 3 < p; k += 4) {
      const double* w0 = W.data() + 2 * (size_t)k * q;
      const double *w1 = w0 + 2 * q, *w2 = w1 + 2 * q, *w3 = w2 + 2 * q;
      ldouble s0 = 0.0L, s1 = 0.0L, s2 = 0.0L, s3 = 0.0L;
      for (int g = 0; g < q; ++g) {
        const int i = idx[g];
        s0 += w0[i];
        s1 += w1[i];
        s2 += w2[i];
        s3 += w3[i];
      }
      out[k * stride] = finish(s0, k);
      out[(k + 1) * stride] = finish(s1, k + 1);
      out[(k + 2) * stride] = finish(s2, k + 2);
      out[(k + 3) * stride] = finish(s3, k + 3);
    }
    for (; k < p; ++k) {
      const double* w = W.data() + 2 * (size_t)k * q;
      ldouble s = 0.0L;
      for (int g = 0; g < q; ++g) s += w[idx[g]];
      out[k * stride] = finish(s, k);
    }
  };

  stats(t0.begin(), 1);
  for (int b = 0; b < B; ++b) {
    if ((b & 255) == 0) Rcpp::checkUserInterrupt();
    unit_indicators(b, q, draw, ind, f);
    stats(tp.begin() + b, B);
  }
  return List::create(_["t0"] = t0, _["t_perm"] = tp);
}

// Paired two-sample swaps without standardization. x1 = data1 - mu, x2 = data2
// (n x p each). Indicator 1 swaps the unit's curves between the two samples
// (fdatest: if_perm <- rbinom(n1, 1, 0.5)). Statistic: alt_stat of the
// difference of the two sample means.
// [[Rcpp::export]]
List cpp_stats_paired(NumericMatrix x1, NumericMatrix x2, IntegerVector unit, int q,
                      int B, IntegerMatrix ind, int alt) {
  const int n = x1.nrow(), p = x1.ncol();
  const bool draw = ind.nrow() == 0;
  const std::vector<double> U1 = unit_sums(x1, unit, q), U2 = unit_sums(x2, unit, q);
  // interleaved: W[2 (k q + g)] = U1, W[2 (k q + g) + 1] = U2, so that the swap
  // is an index offset rather than a branch
  std::vector<double> W(2 * (size_t)q * p);
  for (size_t c = 0; c < (size_t)q * p; ++c) {
    W[2 * c] = U1[c];
    W[2 * c + 1] = U2[c];
  }

  NumericVector t0(p);
  NumericMatrix tp(B, p);
  std::vector<int> f(q, 0);

  auto finish = [&](ldouble s1, ldouble s2) {
    const double m1 = (double)(s1 / n), m2 = (double)(s2 / n);
    return alt_stat(m1 - m2, alt);
  };
  // two columns in lockstep: four independent long double chains
  auto stats = [&](double* out, R_xlen_t stride) {
    int k = 0;
    for (; k + 1 < p; k += 2) {
      const double* w = W.data() + 2 * (size_t)k * q;
      const double* v = w + 2 * (size_t)q;
      ldouble a1 = 0.0L, a2 = 0.0L, b1 = 0.0L, b2 = 0.0L;
      for (int g = 0; g < q; ++g) {
        const int i1 = 2 * g + f[g], i2 = 2 * g + 1 - f[g];
        a1 += w[i1];
        a2 += w[i2];
        b1 += v[i1];
        b2 += v[i2];
      }
      out[k * stride] = finish(a1, a2);
      out[(k + 1) * stride] = finish(b1, b2);
    }
    for (; k < p; ++k) {
      const double* w = W.data() + 2 * (size_t)k * q;
      ldouble s1 = 0.0L, s2 = 0.0L;
      for (int g = 0; g < q; ++g) {
        s1 += w[2 * g + f[g]];
        s2 += w[2 * g + 1 - f[g]];
      }
      out[k * stride] = finish(s1, s2);
    }
  };

  stats(t0.begin(), 1);
  for (int b = 0; b < B; ++b) {
    if ((b & 255) == 0) Rcpp::checkUserInterrupt();
    unit_indicators(b, q, draw, ind, f);
    stats(tp.begin() + b, B);
  }
  return List::create(_["t0"] = t0, _["t_perm"] = tp);
}

// Two-sample statistics at the observation level, as fdatest's
// ts_pointwise_stat() on the permuted coefficient matrix: used for paired swaps
// with standardization and for unpaired label permutations.
// coeff = rbind(data1 - mu, data2) (N x p), first n1 rows = sample 1.
//   paired:   ind is B x q of 0/1 swap indicators per unit (or 0 rows: draw);
//   unpaired: ind is B x N of 1-based permutations (or 0 rows: draw as sample(N)).
// [[Rcpp::export]]
List cpp_stats_twosample_obs(NumericMatrix coeff, int n1, bool paired, IntegerVector unit,
                             int q, int B, IntegerMatrix ind, int alt, bool standardize) {
  const int N = coeff.nrow(), p = coeff.ncol(), n2 = N - n1;
  const bool draw = ind.nrow() == 0;
  std::vector<int> g1(n1), g2(n2), f(paired ? q : 0), perm(paired ? 0 : N), work;
  std::vector<double> v1(n1), v2(n2);
  for (int i = 0; i < n1; ++i) g1[i] = i;
  for (int i = 0; i < n2; ++i) g2[i] = n1 + i;

  NumericVector t0(p);
  NumericMatrix tp(B, p);

  auto stats = [&](double* out, R_xlen_t stride) {
    for (int k = 0; k < p; ++k) {
      const double* col = coeff.begin() + (R_xlen_t)k * N;
      ldouble s1 = 0.0L, s2 = 0.0L;
      for (int i = 0; i < n1; ++i) {
        v1[i] = col[g1[i]];
        s1 += v1[i];
      }
      for (int i = 0; i < n2; ++i) {
        v2[i] = col[g2[i]];
        s2 += v2[i];
      }
      const double m1 = (double)(s1 / n1), m2 = (double)(s2 / n2);
      double t = alt_stat(m1 - m2, alt);
      if (standardize) {
        // sp <- ((n1 - 1) * s1 + (n - n1 - 1) * s2) / (n - 2): R rounds each
        // product; volatile keeps the compiler from fusing them into an FMA
        volatile double a1 = (double)(n1 - 1) * r_cov_var(v1.data(), n1);
        volatile double a2 = (double)(n2 - 1) * r_cov_var(v2.data(), n2);
        const double sp = (a1 + a2) / (double)(N - 2);
        t = t / sp;
      }
      out[k * stride] = t;
    }
  };

  stats(t0.begin(), 1);
  for (int b = 0; b < B; ++b) {
    if ((b & 255) == 0) Rcpp::checkUserInterrupt();
    if (paired) {
      unit_indicators(b, q, draw, ind, f);
      for (int i = 0; i < n1; ++i) {
        const int w = f[unit[i]];
        g1[i] = w ? n1 + i : i;
        g2[i] = w ? i : n1 + i;
      }
    } else {
      if (draw) {
        // sample(N) as R's do_sample() in src/main/random.c (Copyright The R
        // Core Team, GPL >= 2): x[i] = i; iy[i] = x[j] + 1, x[j] = x[--n],
        // j = R_unif_index(n)
        work.resize(N);
        for (int i = 0; i < N; ++i) work[i] = i;
        int dn = N;
        for (int i = 0; i < N; ++i) {
          const int j = (int)R_unif_index((double)dn);
          perm[i] = work[j];
          work[j] = work[--dn];
        }
      } else {
        for (int i = 0; i < N; ++i) perm[i] = ind(b, i) - 1;
      }
      for (int i = 0; i < n1; ++i) g1[i] = perm[i];
      for (int i = 0; i < n2; ++i) g2[i] = perm[n1 + i];
    }
    stats(tp.begin() + b, B);
  }
  return List::create(_["t0"] = t0, _["t_perm"] = tp);
}
