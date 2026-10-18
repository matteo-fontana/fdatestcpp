// P-value engine on precomputed statistics: unadjusted, IWT, TWT and Global.
//
// Inputs everywhere: t0 (length p, observed pointwise statistics) and t_perm
// (B x p, R column-major, resampled pointwise statistics). A resample counts
// against T0 when T* >= T0, as in fdatest.
//
// Ported from fdatest (Pini, Stamm, Vantini; GPL >= 3), development version
// commit 67bb41c: R/utils-adjust-iwt.R, R/utils-adjust-twt.R,
// R/utils-adjust-global.R and src/pval_correct.cpp.

#include "common.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

bool all_finite(const double* x, R_xlen_t n) {
  for (R_xlen_t i = 0; i < n; ++i)
    if (!R_FINITE(x[i])) return false;
  return true;
}

// Pointwise counts #{b : T*_b(k) >= T0(k)}. na[k] flags columns where some
// comparison involves NaN: fdatest's colSums() of the logical matrix is NA.
void pointwise_counts(const double* t0, const double* tp, int B, int p,
                      std::vector<int>& cnt, std::vector<char>& na) {
  cnt.assign(p, 0);
  na.assign(p, 0);
  for (int k = 0; k < p; ++k) {
    const double th = t0[k];
    const double* col = tp + (R_xlen_t)k * B;
    int c = 0;
    bool bad = ISNAN(th);
    for (int b = 0; b < B; ++b) {
      c += (col[b] >= th);
      if (ISNAN(col[b])) bad = true;
    }
    cnt[k] = c;
    na[k] = bad;
  }
}

// ---------------------------------------------------------------------------
// IWT
//
// Interval (j, L) = grid points j, j+1, ..., j+L-1 (0-based; mod p when
// recycling). fdatest stores its p-value in pvalue_matrix[p - L + 1, j + 1].
//
// The interval statistic of a row is accumulated from each start j, extending
// the interval one point at a time: a long double running sum rounded to
// double (exactly what R's sum(x[inf:sup], na.rm = TRUE) does for T0, with its
// overflow to Inf, and rowSums() does for T*) or a running max. That is O(p^2)
// per resample, the same as prefix-sum differences, but reproduces fdatest's
// rounding exactly.
//
// cnt[j * p + L - 1] = #{b : T*_b(j, L) >= T0(j, L)}.
// GENERAL = some input is not finite: NaN is dropped as na.rm = TRUE does
// (0 in sums, -Inf in max) and na flags comparisons that are NaN (Inf - Inf).

// Starts j..j+3 of one resample in lockstep: four independent accumulation
// chains hide the latency of the long double add (or of max). Each chain does
// exactly the single-start arithmetic, so the results are unchanged.
template <bool MAX>
inline void iwt_four_starts(const double* x, const double* th, int* cnt, int p, int j,
                           bool recycle) {
  const double* x0 = x + j;
  const double *x1 = x0 + 1, *x2 = x0 + 2, *x3 = x0 + 3;
  const double* h0 = th + (size_t)j * p;
  const double *h1 = h0 + p, *h2 = h1 + p, *h3 = h2 + p;
  int* c0 = cnt + (size_t)j * p;
  int *c1 = c0 + p, *c2 = c1 + p, *c3 = c2 + p;
  const int Lc = recycle ? p : p - j - 3;  // length range common to the four starts
  if (MAX) {
    double m0 = R_NegInf, m1 = R_NegInf, m2 = R_NegInf, m3 = R_NegInf;
    for (int L = 0; L < Lc; ++L) {
      if (x0[L] > m0) m0 = x0[L];
      if (x1[L] > m1) m1 = x1[L];
      if (x2[L] > m2) m2 = x2[L];
      if (x3[L] > m3) m3 = x3[L];
      c0[L] += (m0 >= h0[L]);
      c1[L] += (m1 >= h1[L]);
      c2[L] += (m2 >= h2[L]);
      c3[L] += (m3 >= h3[L]);
    }
    if (!recycle) {  // starts j+2, j+1, j are 1, 2, 3 points longer
      for (int L = Lc; L < Lc + 1; ++L) { if (x2[L] > m2) m2 = x2[L]; c2[L] += (m2 >= h2[L]); }
      for (int L = Lc; L < Lc + 2; ++L) { if (x1[L] > m1) m1 = x1[L]; c1[L] += (m1 >= h1[L]); }
      for (int L = Lc; L < Lc + 3; ++L) { if (x0[L] > m0) m0 = x0[L]; c0[L] += (m0 >= h0[L]); }
    }
  } else {
    ldouble s0 = 0.0L, s1 = 0.0L, s2 = 0.0L, s3 = 0.0L;
    for (int L = 0; L < Lc; ++L) {
      s0 += x0[L];
      s1 += x1[L];
      s2 += x2[L];
      s3 += x3[L];
      c0[L] += ((double)s0 >= h0[L]);
      c1[L] += ((double)s1 >= h1[L]);
      c2[L] += ((double)s2 >= h2[L]);
      c3[L] += ((double)s3 >= h3[L]);
    }
    if (!recycle) {
      for (int L = Lc; L < Lc + 1; ++L) { s2 += x2[L]; c2[L] += ((double)s2 >= h2[L]); }
      for (int L = Lc; L < Lc + 2; ++L) { s1 += x1[L]; c1[L] += ((double)s1 >= h1[L]); }
      for (int L = Lc; L < Lc + 3; ++L) { s0 += x0[L]; c0[L] += ((double)s0 >= h0[L]); }
    }
  }
}

template <bool MAX, bool GENERAL>
void iwt_counts(const double* t0, const double* tp, int B, int p, bool recycle,
                std::vector<int>& cnt, std::vector<char>& na) {
  const double fill = MAX ? R_NegInf : 0.0;
  const size_t pp = (size_t)p * p;
  std::vector<double> x(2 * (size_t)p), th(pp, 0.0);
  cnt.assign(pp, 0);
  na.assign(GENERAL ? pp : 0, 0);

  for (int k = 0; k < p; ++k) {
    double v = t0[k];
    if (GENERAL && ISNAN(v)) v = fill;
    x[k] = v;
    x[k + p] = v;
  }
  for (int j = 0; j < p; ++j) {
    const int Lmax = recycle ? p : p - j;
    const double* xj = x.data() + j;
    double* thj = th.data() + (size_t)j * p;
    if (MAX) {
      double m = R_NegInf;
      for (int L = 0; L < Lmax; ++L) {
        if (xj[L] > m) m = xj[L];
        thj[L] = m;
      }
    } else {
      ldouble s = 0.0L;
      for (int L = 0; L < Lmax; ++L) {
        s += xj[L];
        thj[L] = r_sum_round(s);  // sum(t0_2x[inf:sup], na.rm = TRUE)
      }
    }
  }

  for (int b = 0; b < B; ++b) {
    if ((b & 255) == 0) Rcpp::checkUserInterrupt();
    for (int k = 0; k < p; ++k) {
      double v = tp[b + (R_xlen_t)k * B];
      if (GENERAL && ISNAN(v)) v = fill;
      x[k] = v;
      x[k + p] = v;
    }
    int j0 = 0;
    if (!GENERAL) {
      for (; j0 + 3 < p; j0 += 4) iwt_four_starts<MAX>(x.data(), th.data(), cnt.data(), p, j0, recycle);
    }
    for (int j = j0; j < p; ++j) {
      const int Lmax = recycle ? p : p - j;
      const double* xj = x.data() + j;
      const double* thj = th.data() + (size_t)j * p;
      int* cj = cnt.data() + (size_t)j * p;
      if (MAX) {
        double m = R_NegInf;
        for (int L = 0; L < Lmax; ++L) {
          if (xj[L] > m) m = xj[L];
          cj[L] += (m >= thj[L]);
        }
      } else {
        ldouble s = 0.0L;
        for (int L = 0; L < Lmax; ++L) {
          s += xj[L];
          const double v = (double)s;
          cj[L] += (v >= thj[L]);
          if (GENERAL && (ISNAN(v) || ISNAN(thj[L]))) na[(size_t)j * p + L] = 1;
        }
      }
    }
  }
}

// Adjusted p-values: adj[k] = max over all intervals containing k of their
// p-values, ignoring NA, which is what fdatest's pval_correct_cpp() returns in
// its first row. Dynamic programme in O(p^2):
//   Q(j, L) = max(P(j, L), Q(j - 1, L + 1), Q(j, L + 1)),  adj[k] = Q(k, 1).
// pm is the p x p column-major p-value matrix in fdatest's layout.
std::vector<double> iwt_correct(const std::vector<double>& pm, int p, bool recycle) {
  auto P = [&](int j, int L) { return pm[(size_t)(p - L) + (size_t)j * p]; };
  std::vector<double> q(p, NA_REAL), qn(p, NA_REAL);
  if (recycle) {
    for (int j = 0; j < p; ++j) q[j] = P(j, p);
    for (int L = p - 1; L >= 1; --L) {
      for (int j = 0; j < p; ++j) {
        const double sup = nanmax_rm(q[(j - 1 + p) % p], q[j]);
        qn[j] = nanmax_rm(P(j, L), sup);
      }
      std::swap(q, qn);
    }
  } else {
    q[0] = P(0, p);
    for (int L = p - 1; L >= 1; --L) {
      const int nj = p - L + 1;  // starts 0..p-L
      for (int j = 0; j < nj; ++j) {
        double sup = NA_REAL;
        if (j >= 1) sup = nanmax_rm(sup, q[j - 1]);
        if (j <= p - L - 1) sup = nanmax_rm(sup, q[j]);
        qn[j] = nanmax_rm(P(j, L), sup);
      }
      std::swap(q, qn);
    }
  }
  return q;
}

}  // namespace

// Binary digits of the accumulator (64 for x86 long double). Results equal
// fdatest's (R's LDOUBLE arithmetic) when this matches R's long double
// precision; tests use it to decide between exact and tolerant comparisons.
// [[Rcpp::export(rng = false)]]
int cpp_acc_digits() {
  int digits = 0;
  ldouble_unit_roundoff(&digits);
  return digits;
}

// [[Rcpp::export(rng = false)]]
NumericVector cpp_pointwise(NumericVector t0, NumericMatrix t_perm, bool plus_one) {
  const int p = t0.size(), B = t_perm.nrow();
  std::vector<int> cnt;
  std::vector<char> na;
  pointwise_counts(t0.begin(), t_perm.begin(), B, p, cnt, na);
  NumericVector out(p);
  for (int k = 0; k < p; ++k) out[k] = na[k] ? NA_REAL : pv_div(cnt[k], B, plus_one);
  return out;
}

// [[Rcpp::export(rng = false)]]
List cpp_iwt(NumericVector t0, NumericMatrix t_perm, bool recycle, bool use_max,
             bool plus_one, bool want_matrix) {
  const int p = t0.size(), B = t_perm.nrow();
  const double* pt0 = t0.begin();
  const double* ptp = t_perm.begin();
  const bool general = !(all_finite(pt0, p) && all_finite(ptp, (R_xlen_t)B * p));

  std::vector<int> pcnt, cnt;
  std::vector<char> pna, na;
  pointwise_counts(pt0, ptp, B, p, pcnt, pna);
  if (use_max) {
    if (general) iwt_counts<true, true>(pt0, ptp, B, p, recycle, cnt, na);
    else iwt_counts<true, false>(pt0, ptp, B, p, recycle, cnt, na);
  } else {
    if (general) iwt_counts<false, true>(pt0, ptp, B, p, recycle, cnt, na);
    else iwt_counts<false, false>(pt0, ptp, B, p, recycle, cnt, na);
  }

  // p-value matrix in fdatest's layout; last row = unadjusted p-values
  std::vector<double> pm((size_t)p * p, NA_REAL);
  NumericVector unadj(p);
  for (int k = 0; k < p; ++k) {
    unadj[k] = pna[k] ? NA_REAL : pv_div(pcnt[k], B, plus_one);
    pm[(size_t)(p - 1) + (size_t)k * p] = unadj[k];
  }
  for (int j = 0; j < p; ++j) {
    const int Lmax = recycle ? p : p - j;
    for (int L = 2; L <= Lmax; ++L) {
      const size_t c = (size_t)j * p + (L - 1);
      pm[(size_t)(p - L) + (size_t)j * p] =
          (general && na[c]) ? NA_REAL : pv_div(cnt[c], B, plus_one);
    }
  }
  std::vector<double> adj = iwt_correct(pm, p, recycle);
  NumericVector adjusted(adj.begin(), adj.end());

  if (!want_matrix) {
    return List::create(_["unadjusted_pvalues"] = unadj,
                        _["adjusted_pvalues"] = adjusted,
                        _["global_pvalue"] = pm[0]);
  }
  NumericMatrix M(p, p);
  std::copy(pm.begin(), pm.end(), M.begin());
  return List::create(_["unadjusted_pvalues"] = unadj,
                      _["adjusted_pvalues"] = adjusted,
                      _["pvalue_matrix"] = M,
                      _["global_pvalue"] = pm[0]);
}

// IWT correction of a given p-value matrix (for tests against pval_correct_cpp)
// [[Rcpp::export(rng = false)]]
NumericVector cpp_iwt_correct(NumericMatrix pvalue_matrix, bool recycle) {
  const int p = pvalue_matrix.nrow();
  std::vector<double> pm(pvalue_matrix.begin(), pvalue_matrix.end());
  std::vector<double> adj = iwt_correct(pm, p, recycle);
  return NumericVector(adj.begin(), adj.end());
}

// ---------------------------------------------------------------------------
// TWT
//
// fdatest's family: for each threshold th in c(0, sort(unique(pval)), 1) the
// sub-level set {pval <= th} and the super-level set {pval > th}. With grid
// points sorted by p-value (ties by index) these are prefixes and suffixes of
// the sorted order: prefixes [0, cb[h]) for h = 1..U and suffixes [cb[h], m)
// for h = 1..U-1, where cb[h] ends the h-th group of equal p-values. The
// thresholds 0 and 1 only add empty sets or the full set (= prefix U).
// adjusted[k] = max(pval[k], p-values of all sets containing k), NA
// propagating as R's max() does.
//
// Set statistics. max: running max along the sorted order (exact). integral:
// fdatest sums each set in grid-index order in long double (sum() for T0,
// rowSums() for T*); the sorted-order running sum differs from that only by
// rounding, bounded by
//   |F - s| <= (2^-53 + 2 (m - 1) u) sum|t|,
// with u the unit roundoff of the accumulator (2^-64 for x86 long double),
// so the comparison with T0 is decided by s unless |s - F0| is within twice
// that bound; then the set is re-summed in index order (exact). The bound
// assumes no overflow, so sets whose sum of |t| exceeds DBL_MAX / 4, or whose
// observed sum is infinite, also take the exact path. Ties between T* and T0
// are the only ordinary cases that reach it, so the cost stays O(B p). With
// non-finite inputs every set is summed in index order.

// [[Rcpp::export(rng = false)]]
List cpp_twt(NumericVector t0, NumericMatrix t_perm, bool use_max, bool plus_one,
             bool brute) {
  const int p = t0.size(), B = t_perm.nrow();
  const double* pt0 = t0.begin();
  const double* ptp = t_perm.begin();

  std::vector<int> pcnt;
  std::vector<char> pna;
  pointwise_counts(pt0, ptp, B, p, pcnt, pna);
  NumericVector unadj(p);
  for (int k = 0; k < p; ++k) unadj[k] = pna[k] ? NA_REAL : pv_div(pcnt[k], B, plus_one);

  // valid grid points sorted by p-value (equivalently by count), ties by index
  std::vector<int> ord;
  for (int k = 0; k < p; ++k)
    if (!pna[k]) ord.push_back(k);
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return pcnt[a] < pcnt[b]; });
  const int m = (int)ord.size();
  std::vector<int> rank(p, -1);
  for (int pos = 0; pos < m; ++pos) rank[ord[pos]] = pos;
  std::vector<int> cb(1, 0);
  for (int pos = 1; pos < m; ++pos)
    if (pcnt[ord[pos]] != pcnt[ord[pos - 1]]) cb.push_back(pos);
  if (m > 0) cb.push_back(m);
  const int U = (int)cb.size() - 1;

  // sum over the set of sorted positions [lo, hi), in grid-index order,
  // rounded as sum() (observed) or rowSums() (resampled)
  auto exact_sum = [&](const double* v, int lo, int hi, bool observed) -> double {
    ldouble s = 0.0L;
    for (int k = 0; k < p; ++k) {
      const int r = rank[k];
      if (r >= lo && r < hi) s += v[k];
    }
    return observed ? r_sum_round(s) : (double)s;
  };

  // observed set statistics
  std::vector<double> f0P(U + 1, 0.0), f0S(U + 1, 0.0);
  if (use_max) {
    double mx = R_NegInf;
    int pos = 0;
    for (int h = 1; h <= U; ++h) {
      for (; pos < cb[h]; ++pos) if (pt0[ord[pos]] > mx) mx = pt0[ord[pos]];
      f0P[h] = mx;
    }
    mx = R_NegInf;
    pos = m;
    for (int h = U - 1; h >= 1; --h) {
      while (pos > cb[h]) {
        --pos;
        if (pt0[ord[pos]] > mx) mx = pt0[ord[pos]];
      }
      f0S[h] = mx;
    }
  } else {
    for (int h = 1; h <= U; ++h) f0P[h] = exact_sum(pt0, 0, cb[h], true);
    for (int h = 1; h < U; ++h) f0S[h] = exact_sum(pt0, cb[h], m, true);
  }

  const bool finite = all_finite(pt0, p) && all_finite(ptp, (R_xlen_t)B * p);
  const bool fast_sum = !use_max && finite && !brute;
  // twice the bound above: 2^-52 + 4 p u
  const ldouble tol_k = (ldouble)std::ldexp(1.0, -52) + 4 * (ldouble)p * ldouble_unit_roundoff();

  std::vector<int> cP(U + 1, 0), cS(U + 1, 0);
  std::vector<char> naP(U + 1, 0), naS(U + 1, 0);
  std::vector<double> r(p);

  // comparison of a set's sorted-order sum s (sum of |t| = a) with F0
  const ldouble big = (ldouble)DBL_MAX / 4;
  auto guarded = [&](ldouble s, ldouble a, double F0, int lo, int hi) -> int {
    if (R_FINITE(F0) && a <= big) {
      const ldouble d = s - (ldouble)F0;
      const ldouble tol = tol_k * a;
      if (d > tol) return 1;
      if (-d > tol) return 0;
      if (a == 0.0L) return 0.0 >= F0;  // every term is zero
    }
    return exact_sum(r.data(), lo, hi, false) >= F0;
  };

  for (int b = 0; b < B; ++b) {
    if ((b & 255) == 0) Rcpp::checkUserInterrupt();
    for (int k = 0; k < p; ++k) r[k] = ptp[b + (R_xlen_t)k * B];

    if (use_max) {
      double mx = R_NegInf;
      int pos = 0;
      for (int h = 1; h <= U; ++h) {
        for (; pos < cb[h]; ++pos) if (r[ord[pos]] > mx) mx = r[ord[pos]];
        cP[h] += (mx >= f0P[h]);
      }
      mx = R_NegInf;
      pos = m;
      for (int h = U - 1; h >= 1; --h) {
        while (pos > cb[h]) {
          --pos;
          if (r[ord[pos]] > mx) mx = r[ord[pos]];
        }
        cS[h] += (mx >= f0S[h]);
      }
    } else if (fast_sum) {
      ldouble s = 0.0L, a = 0.0L;
      int pos = 0;
      for (int h = 1; h <= U; ++h) {
        for (; pos < cb[h]; ++pos) {
          const double v = r[ord[pos]];
          s += v;
          a += std::fabs(v);
        }
        cP[h] += guarded(s, a, f0P[h], 0, cb[h]);
      }
      s = 0.0L;
      a = 0.0L;
      pos = m;
      for (int h = U - 1; h >= 1; --h) {
        while (pos > cb[h]) {
          --pos;
          const double v = r[ord[pos]];
          s += v;
          a += std::fabs(v);
        }
        cS[h] += guarded(s, a, f0S[h], cb[h], m);
      }
    } else {
      for (int h = 1; h <= U; ++h) {
        const double F = exact_sum(r.data(), 0, cb[h], false);
        if (ISNAN(F) || ISNAN(f0P[h])) naP[h] = 1;
        else cP[h] += (F >= f0P[h]);
      }
      for (int h = 1; h < U; ++h) {
        const double F = exact_sum(r.data(), cb[h], m, false);
        if (ISNAN(F) || ISNAN(f0S[h])) naS[h] = 1;
        else cS[h] += (F >= f0S[h]);
      }
    }
  }

  std::vector<double> pvP(U + 1, NA_REAL), pvS(U + 1, NA_REAL);
  for (int h = 1; h <= U; ++h) pvP[h] = naP[h] ? NA_REAL : pv_mean(cP[h], B, plus_one);
  for (int h = 1; h < U; ++h) pvS[h] = naS[h] ? NA_REAL : pv_mean(cS[h], B, plus_one);

  // group g (1..U) lies in prefixes h >= g and in suffixes h <= g - 1
  std::vector<double> preMax(U + 2, R_NegInf), sufMax(U + 1, R_NegInf);
  for (int h = U; h >= 1; --h) preMax[h] = nanmax_prop(preMax[h + 1], pvP[h]);
  for (int g = 2; g <= U; ++g) sufMax[g] = nanmax_prop(sufMax[g - 1], pvS[g - 1]);
  NumericVector adjusted(p, NA_REAL);
  for (int g = 1; g <= U; ++g) {
    const double pg = unadj[ord[cb[g - 1]]];
    const double ag = nanmax_prop(nanmax_prop(pg, preMax[g]), sufMax[g]);
    for (int pos = cb[g - 1]; pos < cb[g]; ++pos) adjusted[ord[pos]] = ag;
  }

  return List::create(_["unadjusted_pvalues"] = unadj,
                      _["adjusted_pvalues"] = adjusted,
                      _["global_pvalue"] = (U >= 1) ? pvP[U] : NA_REAL);
}

// ---------------------------------------------------------------------------
// Global (fdatest's p_adjust_global): sum(t0) vs rowSums(t_coeff), or max,
// without na.rm, so any NaN makes the p-value NA.

// [[Rcpp::export(rng = false)]]
List cpp_global(NumericVector t0, NumericMatrix t_perm, bool use_max, bool plus_one) {
  const int p = t0.size(), B = t_perm.nrow();
  const double* pt0 = t0.begin();
  const double* ptp = t_perm.begin();

  std::vector<int> pcnt;
  std::vector<char> pna;
  pointwise_counts(pt0, ptp, B, p, pcnt, pna);
  NumericVector unadj(p);
  for (int k = 0; k < p; ++k) unadj[k] = pna[k] ? NA_REAL : pv_div(pcnt[k], B, plus_one);

  bool bad = false;
  int c = 0;
  if (use_max) {
    double F0 = R_NegInf;
    for (int k = 0; k < p; ++k) {
      if (ISNAN(pt0[k])) bad = true;
      else if (pt0[k] > F0) F0 = pt0[k];
    }
    std::vector<double> mx(B, R_NegInf);
    for (int k = 0; k < p; ++k) {
      const double* col = ptp + (R_xlen_t)k * B;
      for (int b = 0; b < B; ++b) {
        if (ISNAN(col[b])) bad = true;
        else if (col[b] > mx[b]) mx[b] = col[b];
      }
    }
    if (!bad)
      for (int b = 0; b < B; ++b) c += (mx[b] >= F0);
  } else {
    ldouble s0 = 0.0L;
    for (int k = 0; k < p; ++k) s0 += pt0[k];
    const double F0 = r_sum_round(s0);  // sum(t0)
    std::vector<ldouble> acc(B, 0.0L);
    for (int k = 0; k < p; ++k) {
      const double* col = ptp + (R_xlen_t)k * B;
      for (int b = 0; b < B; ++b) acc[b] += col[b];
    }
    bad = ISNAN(F0);
    for (int b = 0; b < B; ++b) {
      const double F = (double)acc[b];
      if (ISNAN(F)) bad = true;
      c += (F >= F0);
    }
  }
  const double glob = bad ? NA_REAL : pv_mean(c, B, plus_one);
  return List::create(_["unadjusted_pvalues"] = unadj,
                      _["adjusted_pvalues"] = NumericVector(p, glob),
                      _["global_pvalue"] = glob);
}
