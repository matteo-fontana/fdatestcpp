// Shared helpers for fdatestcpp.
//
// Bit-exactness with fdatest: fdatest computes everything with R primitives
// (colMeans, rowSums, sum, cov, mean) that accumulate in LDOUBLE (80-bit
// extended precision with gcc/mingw) and round to double at the end. The code
// here accumulates in `long double` in the same order and rounds at the same
// points, so with identical inputs the doubles it produces are identical.

#ifndef FDATESTCPP_COMMON_H
#define FDATESTCPP_COMMON_H

#include <Rcpp.h>
#include <cfloat>
#include <cmath>
#include <limits>

// Accumulator type: long double, as R's LDOUBLE. Its precision is whatever
// the platform provides (80-bit extended on x86_64, double on arm64 macOS,
// quadruple on aarch64 Linux); the code only assumes it is at least double.
// FDATESTCPP_DOUBLE_ACC forces double, to test that case on x86_64.
#ifdef FDATESTCPP_DOUBLE_ACC
typedef double ldouble;
#else
typedef long double ldouble;
#endif

// Unit roundoff of the accumulator, measured at run time (as R measures
// .Machine$longdouble.eps) rather than taken from std::numeric_limits: under
// valgrind, which has no 80-bit arithmetic, or with a reduced x87 precision
// control word, long double arithmetic runs below its nominal precision.
// digits receives the number of binary digits (64 for x86 long double).
inline ldouble ldouble_unit_roundoff(int* digits = nullptr) {
  volatile ldouble one = 1, u = 1, s;
  int d = 0;
  do {
    u = u / 2;
    ++d;
    s = one + u;
  } while (s != one);
  if (digits) *digits = d;
  return u;
}

// Rounding of a long double sum. fdatest sums observed statistics with sum()
// and resampled ones with rowSums(). R's sum() (rsum() in summary.c) maps a
// total beyond the double range to +/-Inf; rowSums(), like a plain cast,
// rounds to nearest, so totals up to DBL_MAX + 2^970 become DBL_MAX.
inline double r_sum_round(ldouble s) {
  if (s > DBL_MAX) return R_PosInf;
  if (s < -DBL_MAX) return R_NegInf;
  return (double)s;
}

// p-value from a count of resamples with T* >= T0.
// fdatest computes it in two ways, which can differ in the last bit:
//   IWT matrix and unadjusted p-values: integer count / n_perm (double division)
//   TWT and Global: mean(logical), i.e. (double)((long double)count / B)
// plus_one = TRUE gives (1 + count) / (B + 1) instead.
inline double pv_div(int cnt, int B, bool plus_one) {
  if (plus_one) return ((double)cnt + 1.0) / ((double)B + 1.0);
  return (double)cnt / (double)B;
}

inline double pv_mean(int cnt, int B, bool plus_one) {
  if (plus_one) return ((double)cnt + 1.0) / ((double)B + 1.0);
  return (double)((ldouble)cnt / (ldouble)B);
}

// max() ignoring NaN (fdatest's pval_correct_cpp and na.rm = TRUE)
inline double nanmax_rm(double a, double b) {
  if (ISNAN(a)) return b;
  if (ISNAN(b)) return a;
  return a > b ? a : b;
}

// max() propagating NaN (R's max() without na.rm)
inline double nanmax_prop(double a, double b) {
  if (ISNAN(a) || ISNAN(b)) return NA_REAL;
  return a > b ? a : b;
}

#endif
