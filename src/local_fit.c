#include <R.h>
#include <Rinternals.h>
/* No fused multiply-add contraction (clang contracts a*b + c by default on
 * arm64), so that results are identical across compilers and CPUs. */
#ifdef __clang__
#pragma clang fp contract(off)
#endif
#include <math.h>

/*
 * Direct local polynomial regression, evaluated exactly at each requested
 * point (locfit evaluates on an adaptive tree and interpolates instead).
 *
 * At x0 the bandwidth is max(h, d_k), where d_k is the distance to the k-th
 * nearest fitting point and k = n * nn (truncated, at least 1), as in
 * locfit. Points with
 * |x - x0| < bandwidth get weight w * (1 - |u|^3)^3, u = (x - x0) / bandwidth,
 * and a weighted polynomial of degree `degree` in u is fitted; its intercept
 * is the fitted value. If the local design is singular the degree is
 * lowered (quadratic, then linear, then a weighted mean).
 *
 * x must be sorted in increasing order; weights must be positive.
 */

/* Solve the (p x p) system a * beta = b in place by Gaussian elimination with
 * partial pivoting. Returns 0 if a pivot is negligible relative to scale. */
static int solve_small(double a[3][3], double b[3], int p, double scale)
{
    int i, j, k, pivot;
    for (k = 0; k < p; k++) {
        double best = fabs(a[k][k]), tmp;
        pivot = k;
        for (i = k + 1; i < p; i++) {
            if (fabs(a[i][k]) > best) { best = fabs(a[i][k]); pivot = i; }
        }
        if (!(best > 1e-10 * scale)) return 0;
        if (pivot != k) {
            for (j = 0; j < p; j++) { tmp = a[k][j]; a[k][j] = a[pivot][j]; a[pivot][j] = tmp; }
            tmp = b[k]; b[k] = b[pivot]; b[pivot] = tmp;
        }
        for (i = k + 1; i < p; i++) {
            double factor = a[i][k] / a[k][k];
            for (j = k; j < p; j++) a[i][j] -= factor * a[k][j];
            b[i] -= factor * b[k];
        }
    }
    for (k = p - 1; k >= 0; k--) {
        for (j = k + 1; j < p; j++) b[k] -= a[k][j] * b[j];
        b[k] /= a[k][k];
    }
    return 1;
}

/* First index with x[index] >= value (x sorted). */
static R_xlen_t lower_bound(const double *x, R_xlen_t n, double value)
{
    R_xlen_t left = 0, right = n, mid;
    while (left < right) {
        mid = left + (right - left) / 2;
        if (x[mid] < value) left = mid + 1; else right = mid;
    }
    return left;
}

SEXP C_local_fit(SEXP x_sexp, SEXP y_sexp, SEXP w_sexp, SEXP eval_sexp,
                 SEXP nn_sexp, SEXP h_sexp, SEXP degree_sexp)
{
    R_xlen_t n, n_eval, e, i, k, left, right, first, last;
    const double *x, *y, *w, *x_eval;
    double nn, h_fixed, *out;
    int degree;
    SEXP result;

    if (TYPEOF(x_sexp) != REALSXP || TYPEOF(y_sexp) != REALSXP ||
        TYPEOF(w_sexp) != REALSXP || TYPEOF(eval_sexp) != REALSXP)
        error("Local fit inputs must be double vectors.");
    n = XLENGTH(x_sexp);
    if (XLENGTH(y_sexp) != n || XLENGTH(w_sexp) != n)
        error("Local fit inputs must have equal lengths.");
    if (n < 1) error("Local fit needs at least one point.");
    n_eval = XLENGTH(eval_sexp);
    x = REAL(x_sexp); y = REAL(y_sexp); w = REAL(w_sexp); x_eval = REAL(eval_sexp);
    nn = asReal(nn_sexp); h_fixed = asReal(h_sexp); degree = asInteger(degree_sexp);
    if (!R_FINITE(nn) || nn <= 0 || !R_FINITE(h_fixed) || h_fixed < 0 ||
        degree < 0 || degree > 2)
        error("Invalid local fit bandwidth or degree.");
    for (i = 0; i < n; i++) {
        if (!R_FINITE(x[i]) || !R_FINITE(y[i]) || !R_FINITE(w[i]) || w[i] <= 0)
            error("Local fit inputs must be finite with positive weights.");
        if (i > 0 && x[i] < x[i - 1]) error("Local fit positions must be sorted.");
    }
    /* Tolerate rounding so that, e.g., nn = 30 / n gives k = 30 as in locfit. */
    k = (R_xlen_t) (n * nn + 1e-8);
    if (k < 1) k = 1;
    if (k > n) k = n;

    PROTECT(result = allocVector(REALSXP, n_eval));
    out = REAL(result);
    for (e = 0; e < n_eval; e++) {
        double x0 = x_eval[e], bandwidth, d_k = 0.0;
        double moments[5] = {0, 0, 0, 0, 0}, rhs[3] = {0, 0, 0};
        int p, solved = 0;

        /* Distance to the k-th nearest fitting point. */
        right = lower_bound(x, n, x0);
        left = right - 1;
        for (i = 0; i < k; i++) {
            if (left >= 0 && (right >= n || x0 - x[left] <= x[right] - x0)) {
                d_k = x0 - x[left]; left--;
            } else {
                d_k = x[right] - x0; right++;
            }
        }
        bandwidth = d_k > h_fixed ? d_k : h_fixed;
        if (!(bandwidth > 0)) bandwidth = 1.0;

        first = lower_bound(x, n, x0 - bandwidth);
        last = lower_bound(x, n, x0 + bandwidth);
        for (i = first; i < last; i++) {
            double u = (x[i] - x0) / bandwidth, au = fabs(u), kernel, wu, u2;
            if (au >= 1.0) continue;
            kernel = 1.0 - au * au * au;
            wu = w[i] * kernel * kernel * kernel;
            u2 = u * u;
            moments[0] += wu; moments[1] += wu * u; moments[2] += wu * u2;
            moments[3] += wu * u2 * u; moments[4] += wu * u2 * u2;
            rhs[0] += wu * y[i]; rhs[1] += wu * y[i] * u; rhs[2] += wu * y[i] * u2;
        }
        if (!(moments[0] > 0)) {
            out[e] = NA_REAL;
            continue;
        }
        for (p = degree + 1; p >= 1 && !solved; p--) {
            double a[3][3], b[3];
            int r, c;
            for (r = 0; r < p; r++) {
                for (c = 0; c < p; c++) a[r][c] = moments[r + c];
                b[r] = rhs[r];
            }
            if (solve_small(a, b, p, moments[0])) {
                out[e] = b[0];
                solved = 1;
            }
        }
        if (!solved) out[e] = NA_REAL;
        if ((e & 65535) == 0) R_CheckUserInterrupt();
    }
    UNPROTECT(1);
    return result;
}
