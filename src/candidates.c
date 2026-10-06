#include <R.h>
#include <Rinternals.h>
/* No fused multiply-add contraction (clang contracts a*b + c by default on
 * arm64), so that results are identical across compilers and CPUs. */
#ifdef __clang__
#pragma clang fp contract(off)
#endif
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

/*
 * Compiled versions of the candidate-discovery steps. They follow the R
 * reference implementations (.compute_candidate_smoothing_stats_r,
 * .refine_candidate_segments_r and .trim_candidate_shape) step by step.
 */

/* ---------------------------------------------------------------- helpers */

static void sort_doubles(double *x, int n)
{
    int i, j;
    for (i = 1; i < n; i++) {
        double value = x[i];
        for (j = i - 1; j >= 0 && x[j] > value; j--) x[j + 1] = x[j];
        x[j + 1] = value;
    }
}

/* Median as in matrixStats: the mean of the two middle values for even n. */
static double median_sorted(const double *x, int n)
{
    int half = n / 2;
    if (n % 2) return x[half];
    return (x[half - 1] + x[half]) / 2.0;
}

/* matrixStats::rowMads(constant = 1) of the n values in x (overwritten). */
static double mad_of(double *x, int n)
{
    double center;
    int i;
    if (n == 0) return NA_REAL;
    sort_doubles(x, n);
    center = median_sorted(x, n);
    for (i = 0; i < n; i++) x[i] = fabs(x[i] - center);
    sort_doubles(x, n);
    return median_sorted(x, n);
}

/* mean() as in base R: long double sum, then one refinement pass. */
static double r_mean(const double *x, R_xlen_t n)
{
    long double s = 0.0, t = 0.0;
    R_xlen_t i;
    for (i = 0; i < n; i++) s += x[i];
    s /= n;
    if (R_FINITE((double) s)) {
        for (i = 0; i < n; i++) t += (x[i] - s);
        s += t / n;
    }
    return (double) s;
}

/* ------------------------------------------------------- smoothing stats */

/*
 * Locus statistics of one chromosome from dense matrices: `pos` holds the
 * sorted locus positions and S/T are (locus x sample) column-major matrices
 * with zeros for uncovered samples. `sample_group` gives 1 or 2 for each
 * sample column. Returns the columns of the R table plus `invalid_weight`,
 * the number of loci whose effect is finite but whose smoothing weight is
 * not positive and finite.
 */
SEXP C_smoothing_stats(SEXP pos_sexp, SEXP s_sexp, SEXP t_sexp, SEXP group_sexp,
                       SEXP cov_q75_sexp)
{
    R_xlen_t n_loci = XLENGTH(pos_sexp), locus;
    int n_samples = (int) XLENGTH(group_sexp), j;
    const double *pos = REAL(pos_sexp), *S_in = REAL(s_sexp), *T_in = REAL(t_sexp);
    const int *group = INTEGER(group_sexp);
    double cov_q75 = asReal(cov_q75_sexp), invalid = 0.0;
    double *S_row, *T_row, *values;
    const char *names[] = {"pos", "group1_pooled", "group2_pooled", "beta",
        "beta_raw_std", "mad_group1", "mad_group2", "sd_raw", "cov_mean",
        "weight", "n_covered_group1", "n_covered_group2", "invalid_weight"};
    int n_columns = 13, k;
    double *column[12];
    SEXP result, result_names;

    if (XLENGTH(s_sexp) != n_loci * n_samples || XLENGTH(t_sexp) != n_loci * n_samples)
        error("S and T must be (locus x sample) matrices.");
    for (locus = 1; locus < n_loci; locus++)
        if (!(pos[locus] > pos[locus - 1])) error("Loci must be unique and sorted.");

    PROTECT(result = allocVector(VECSXP, n_columns));
    PROTECT(result_names = allocVector(STRSXP, n_columns));
    for (k = 0; k < 12; k++) {
        SET_VECTOR_ELT(result, k, allocVector(REALSXP, n_loci));
        column[k] = REAL(VECTOR_ELT(result, k));
    }
    for (k = 0; k < n_columns; k++) SET_STRING_ELT(result_names, k, mkChar(names[k]));
    setAttrib(result, R_NamesSymbol, result_names);

    S_row = (double *) R_alloc(n_samples, sizeof(double));
    T_row = (double *) R_alloc(n_samples, sizeof(double));
    values = (double *) R_alloc(n_samples, sizeof(double));
    for (locus = 0; locus < n_loci; locus++) {
        double s1 = 0, t1 = 0, s2 = 0, t2 = 0, total = 0;
        double pooled1, pooled2, mad1, mad2, beta, sd_raw, cov_mean, weight;
        double covered1 = 0, covered2 = 0, floor_value;
        int n1 = 0, n2 = 0;
        for (j = 0; j < n_samples; j++) {
            S_row[j] = S_in[locus + (R_xlen_t) j * n_loci];
            T_row[j] = T_in[locus + (R_xlen_t) j * n_loci];
        }
        for (j = 0; j < n_samples; j++) {
            total += T_row[j];
            if (group[j] == 1) {
                s1 += S_row[j]; t1 += T_row[j];
                if (T_row[j] > 0) covered1 += 1;
            } else if (group[j] == 2) {
                s2 += S_row[j]; t2 += T_row[j];
                if (T_row[j] > 0) covered2 += 1;
            }
        }
        pooled1 = t1 <= 0 ? NA_REAL : s1 / t1;
        pooled2 = t2 <= 0 ? NA_REAL : s2 / t2;
        for (j = 0; j < n_samples; j++)
            if (group[j] == 1 && T_row[j] > 0) values[n1++] = S_row[j] / T_row[j];
        mad1 = mad_of(values, n1);
        if (!ISNA(mad1)) mad1 = mad1 * 1.4826;
        for (j = 0; j < n_samples; j++)
            if (group[j] == 2 && T_row[j] > 0) values[n2++] = S_row[j] / T_row[j];
        mad2 = mad_of(values, n2);
        if (!ISNA(mad2)) mad2 = mad2 * 1.4826;
        if (covered1 == 1 && covered2 == 1) {
            mad1 = 1; mad2 = 1;
        } else if (covered1 == 1 && covered2 >= 2) {
            mad1 = mad2;
        } else if (covered2 == 1 && covered1 >= 2) {
            mad2 = mad1;
        }
        beta = (ISNA(pooled1) || ISNA(pooled2)) ? NA_REAL : pooled1 - pooled2;
        /* Intentional source parity with dmrseq: both 1.4826 factors. */
        if (ISNA(mad1) || ISNA(mad2)) {
            sd_raw = NA_REAL;
        } else {
            /* Separate statements: no fused multiply-add (clang on arm64
             * contracts within an expression), so results match R exactly. */
            double square1 = mad1 * mad1, square2 = mad2 * mad2;
            sd_raw = 1.4826 * sqrt(square1 + square2);
        }
        if (R_FINITE(sd_raw) && sd_raw < 1e-5) sd_raw = 1e-5;
        cov_mean = total / n_samples;
        floor_value = 1 / (cov_mean > 5 ? cov_mean : 5);
        weight = ISNA(sd_raw) ? NA_REAL :
            (cov_mean < cov_q75 ? cov_mean : cov_q75) /
            (sd_raw > floor_value ? sd_raw : floor_value);
        if (R_FINITE(beta) && (!R_FINITE(weight) || weight <= 0)) invalid += 1;
        column[0][locus] = pos[locus];
        column[1][locus] = pooled1;
        column[2][locus] = pooled2;
        column[3][locus] = beta;
        column[4][locus] = (ISNA(beta) || ISNA(sd_raw)) ? NA_REAL :
            beta / (sd_raw * 2 / sqrt((double) n_samples));
        column[5][locus] = mad1;
        column[6][locus] = mad2;
        column[7][locus] = sd_raw;
        column[8][locus] = cov_mean;
        column[9][locus] = weight;
        column[10][locus] = covered1;
        column[11][locus] = covered2;
    }
    SET_VECTOR_ELT(result, 12, ScalarReal(invalid));
    UNPROTECT(2);
    return result;
}

/*
 * Region sums of S and T over the loci with start - 1 < pos <= end (as
 * findInterval() on prefix sums), for (locus x sample) matrices.
 */
SEXP C_region_sums(SEXP pos_sexp, SEXP s_sexp, SEXP t_sexp, SEXP start_sexp, SEXP end_sexp)
{
    R_xlen_t n_loci = XLENGTH(pos_sexp), n_regions = XLENGTH(start_sexp), r, i;
    int n_samples, j;
    const double *pos = REAL(pos_sexp), *S_in = REAL(s_sexp), *T_in = REAL(t_sexp),
        *start = REAL(start_sexp), *end = REAL(end_sexp);
    double *S_out, *T_out;
    SEXP result, S_sexp, T_sexp, names;

    n_samples = n_loci ? (int) (XLENGTH(s_sexp) / n_loci) : 0;
    if (XLENGTH(t_sexp) != XLENGTH(s_sexp) || XLENGTH(end_sexp) != n_regions)
        error("Invalid region-sum inputs.");
    PROTECT(S_sexp = allocMatrix(REALSXP, (int) n_regions, n_samples));
    PROTECT(T_sexp = allocMatrix(REALSXP, (int) n_regions, n_samples));
    S_out = REAL(S_sexp);
    T_out = REAL(T_sexp);
    for (r = 0; r < n_regions; r++) {
        R_xlen_t left = 0, right = n_loci, mid, first;
        double low = start[r] - 1;
        while (left < right) {
            mid = left + (right - left) / 2;
            if (pos[mid] <= low) left = mid + 1; else right = mid;
        }
        first = left;
        for (j = 0; j < n_samples; j++) {
            long double s = 0, t = 0;
            for (i = first; i < n_loci && pos[i] <= end[r]; i++) {
                s += S_in[i + (R_xlen_t) j * n_loci];
                t += T_in[i + (R_xlen_t) j * n_loci];
            }
            S_out[r + (R_xlen_t) j * n_regions] = (double) s;
            T_out[r + (R_xlen_t) j * n_regions] = (double) t;
        }
    }
    PROTECT(result = allocVector(VECSXP, 2));
    PROTECT(names = allocVector(STRSXP, 2));
    SET_VECTOR_ELT(result, 0, S_sexp);
    SET_VECTOR_ELT(result, 1, T_sexp);
    SET_STRING_ELT(names, 0, mkChar("S"));
    SET_STRING_ELT(names, 1, mkChar("T"));
    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(4);
    return result;
}

/* ------------------------------------------------- candidate refinement */

static int direction_of(double x, double threshold)
{
    int result;
    result = (x >= threshold || fabs(x - threshold) <= sqrt(DBL_EPSILON)) ? 1 : 0;
    if (x <= -threshold) result = -1;
    return result;
}

/* Slope test of lm(x[index] ~ index) for index = first..last (1-based):
 * TRUE when slope * sign > 0 and its two-sided p-value < 0.01. */
static int significant_slope(const double *x, int first, int last, int sign)
{
    int n = last - first + 1, i;
    long double mean_i = 0, mean_x = 0, sxx = 0, sxy = 0, rss = 0;
    double slope, intercept, se, statistic, p_value;
    if (n <= 4) return 0;
    for (i = first; i <= last; i++) { mean_i += i; mean_x += x[i - 1]; }
    mean_i /= n; mean_x /= n;
    for (i = first; i <= last; i++) {
        sxx += (i - mean_i) * (i - mean_i);
        sxy += (i - mean_i) * (x[i - 1] - mean_x);
    }
    slope = (double) (sxy / sxx);
    intercept = (double) (mean_x - slope * mean_i);
    for (i = first; i <= last; i++) {
        long double residual = x[i - 1] - (intercept + slope * i);
        rss += residual * residual;
    }
    se = sqrt((double) (rss / (n - 2)) / (double) sxx);
    statistic = slope / se;
    p_value = 2 * pt(-fabs(statistic), n - 2, 1, 0);
    return R_FINITE(p_value) && slope * sign > 0 && p_value < 0.01;
}

/* .trim_candidate_shape(): returns the retained 1-based range [*first, *last]. */
static void trim_shape(const double *x, int n, int min_cpgs, int *first_out, int *last_out)
{
    int mid = 1, first = 1, last = n, i;
    double maximum = x[0], minimum = x[0], ratio, cut;
    *first_out = 1;
    *last_out = n;
    if (n <= min_cpgs) return;
    for (i = 1; i < n; i++) {
        if (x[i] > maximum) { maximum = x[i]; mid = i + 1; }
        if (x[i] < minimum) minimum = x[i];
    }
    ratio = maximum / minimum;
    if (ISNAN(ratio) || ratio <= 4.0 / 3.0) return;
    cut = (0.5 * (maximum - minimum) + minimum + 0.75 * r_mean(x, n)) / 2;
    if (significant_slope(x, 1, mid, 1)) {
        double a = nearbyint(mid - 0.125 * n), b = mid - 2, c = R_PosInf, value;
        if (a < 1) a = 1;
        if (b < 1) b = 1;
        for (i = 1; i <= mid; i++) if (x[i - 1] >= cut) { c = i; break; }
        value = a < b ? a : b;
        if (c < value) value = c;
        first = (int) value;
    }
    if (significant_slope(x, mid, n, -1)) {
        double a = nearbyint(mid + 0.125 * n), b = mid + 2, c = R_NegInf, value;
        if (a > n) a = n;
        if (b > n) b = n;
        for (i = n; i >= mid; i--) if (x[i - 1] >= cut) { c = i; break; }
        value = a > b ? a : b;
        if (c > value) value = c;
        last = (int) value;
    }
    if (last - first + 1 >= min_cpgs) {
        *first_out = first;
        *last_out = last;
    }
}

SEXP C_trim_shape(SEXP x_sexp, SEXP min_cpgs_sexp)
{
    int first, last, i;
    SEXP result;
    trim_shape(REAL(x_sexp), (int) XLENGTH(x_sexp), asInteger(min_cpgs_sexp), &first, &last);
    PROTECT(result = allocVector(INTSXP, last - first + 1));
    for (i = first; i <= last; i++) INTEGER(result)[i - first] = i;
    UNPROTECT(1);
    return result;
}

/*
 * Candidate segments of loci sorted by (chromosome code, position). Region
 * clusters use all loci; loci without a finite scan value are then dropped,
 * consecutive loci with equal (chromosome, cluster, direction) form a
 * segment, and non-zero segments are refined by the standardized raw signal
 * and trimmed. Returns one row per segment with at least one retained locus.
 */
SEXP C_refine_segments(SEXP chr_sexp, SEXP pos_sexp, SEXP smooth_sexp, SEXP raw_sexp,
                       SEXP block_sexp, SEXP smoothed_sexp, SEXP threshold_sexp,
                       SEXP max_gap_sexp, SEXP min_cpgs_sexp)
{
    R_xlen_t n = XLENGTH(pos_sexp), i, m = 0, start, end, n_out = 0, capacity = 64;
    const int *chr = INTEGER(chr_sexp), *block = INTEGER(block_sexp),
        *smoothed = LOGICAL(smoothed_sexp);
    const double *pos = REAL(pos_sexp), *smooth = REAL(smooth_sexp), *raw = REAL(raw_sexp);
    double threshold = asReal(threshold_sexp), max_gap = asReal(max_gap_sexp);
    int min_cpgs = asInteger(min_cpgs_sexp), cluster = 0, segment = 0;
    R_xlen_t *keep;
    int *cluster_of, *direction, *segment_of, *dir_raw;
    double *x;
    /* outputs */
    int *o_chr, *o_segment, *o_direction, *o_n, *o_block, *o_smoothed;
    double *o_start, *o_end, *o_max, *o_mean;
    const char *names[] = {"chr_code", "segment", "direction", "region_start",
        "region_end", "max_abs_smooth", "mean_smooth", "n_loci", "block", "n_smoothed"};
    SEXP result, result_names;
    int k;

    for (i = 1; i < n; i++) {
        if (chr[i] < chr[i - 1] || (chr[i] == chr[i - 1] && pos[i] < pos[i - 1]))
            error("Loci must be sorted by chromosome and position.");
    }
    keep = (R_xlen_t *) R_alloc(n > 0 ? n : 1, sizeof(R_xlen_t));
    cluster_of = (int *) R_alloc(n > 0 ? n : 1, sizeof(int));
    for (i = 0; i < n; i++) {
        if (i == 0 || chr[i] != chr[i - 1]) cluster = 1;
        else if (pos[i] - pos[i - 1] > max_gap) cluster++;
        cluster_of[i] = cluster;
        if (R_FINITE(smooth[i])) keep[m++] = i;
    }
    direction = (int *) R_alloc(m > 0 ? m : 1, sizeof(int));
    segment_of = (int *) R_alloc(m > 0 ? m : 1, sizeof(int));
    dir_raw = (int *) R_alloc(m > 0 ? m : 1, sizeof(int));
    x = (double *) R_alloc(m > 0 ? m : 1, sizeof(double));
    for (i = 0; i < m; i++) {
        R_xlen_t r = keep[i];
        direction[i] = direction_of(smooth[r], threshold);
        /* NA raw values never match a direction. */
        dir_raw[i] = R_FINITE(raw[r]) ? direction_of(raw[r], threshold) : 99;
        if (i == 0 || chr[r] != chr[keep[i - 1]] || cluster_of[r] != cluster_of[keep[i - 1]] ||
            direction[i] != direction[i - 1]) segment++;
        segment_of[i] = segment;
    }

    o_chr = (int *) R_alloc(capacity, sizeof(int));
    o_segment = (int *) R_alloc(capacity, sizeof(int));
    o_direction = (int *) R_alloc(capacity, sizeof(int));
    o_n = (int *) R_alloc(capacity, sizeof(int));
    o_block = (int *) R_alloc(capacity, sizeof(int));
    o_smoothed = (int *) R_alloc(capacity, sizeof(int));
    o_start = (double *) R_alloc(capacity, sizeof(double));
    o_end = (double *) R_alloc(capacity, sizeof(double));
    o_max = (double *) R_alloc(capacity, sizeof(double));
    o_mean = (double *) R_alloc(capacity, sizeof(double));

    for (start = 0; start < m; start = end) {
        int dir, retained_first, retained_last, length, j;
        double maximum = 0;
        int n_smooth = 0;
        for (end = start; end < m && segment_of[end] == segment_of[start]; end++) ;
        dir = direction[start];
        if (dir == 0) continue;
        length = (int) (end - start);
        retained_first = 1;
        retained_last = length;
        if (length > min_cpgs) {
            int lo = 0, hi = 0;
            for (j = 1; j <= length; j++) {
                if (dir_raw[start + j - 1] == dir) { if (!lo) lo = j; hi = j; }
            }
            if (lo && hi - lo + 1 >= min_cpgs) { retained_first = lo; retained_last = hi; }
        }
        if (retained_last - retained_first + 1 > min_cpgs) {
            int count = retained_last - retained_first + 1, tf, tl;
            for (j = 0; j < count; j++)
                x[j] = smooth[keep[start + retained_first - 1 + j]] * dir;
            trim_shape(x, count, min_cpgs, &tf, &tl);
            retained_last = retained_first + tl - 1;
            retained_first = retained_first + tf - 1;
        }
        length = retained_last - retained_first + 1;
        for (j = 0; j < length; j++) {
            R_xlen_t r = keep[start + retained_first - 1 + j];
            double value = fabs(smooth[r]);
            if (j == 0 || value > maximum) maximum = value;
            x[j] = smooth[r];
            if (smoothed[r] == TRUE) n_smooth++;
        }
        if (n_out == capacity) {
            R_xlen_t next = capacity * 2;
#define GROW(field, type) do { type *g = (type *) R_alloc(next, sizeof(type)); \
            memcpy(g, field, capacity * sizeof(type)); field = g; } while (0)
            GROW(o_chr, int); GROW(o_segment, int); GROW(o_direction, int);
            GROW(o_n, int); GROW(o_block, int); GROW(o_smoothed, int);
            GROW(o_start, double); GROW(o_end, double); GROW(o_max, double);
            GROW(o_mean, double);
#undef GROW
            capacity = next;
        }
        o_chr[n_out] = chr[keep[start]];
        o_segment[n_out] = segment_of[start];
        o_direction[n_out] = dir;
        o_start[n_out] = pos[keep[start + retained_first - 1]];
        o_end[n_out] = pos[keep[start + retained_last - 1]];
        o_max[n_out] = maximum;
        o_mean[n_out] = r_mean(x, length);
        o_n[n_out] = length;
        o_block[n_out] = block[keep[start + retained_first - 1]];
        o_smoothed[n_out] = n_smooth;
        n_out++;
    }

    PROTECT(result = allocVector(VECSXP, 10));
    PROTECT(result_names = allocVector(STRSXP, 10));
    {
        SEXP v;
        v = allocVector(INTSXP, n_out); SET_VECTOR_ELT(result, 0, v);
        memcpy(INTEGER(v), o_chr, n_out * sizeof(int));
        v = allocVector(INTSXP, n_out); SET_VECTOR_ELT(result, 1, v);
        memcpy(INTEGER(v), o_segment, n_out * sizeof(int));
        v = allocVector(INTSXP, n_out); SET_VECTOR_ELT(result, 2, v);
        memcpy(INTEGER(v), o_direction, n_out * sizeof(int));
        v = allocVector(REALSXP, n_out); SET_VECTOR_ELT(result, 3, v);
        memcpy(REAL(v), o_start, n_out * sizeof(double));
        v = allocVector(REALSXP, n_out); SET_VECTOR_ELT(result, 4, v);
        memcpy(REAL(v), o_end, n_out * sizeof(double));
        v = allocVector(REALSXP, n_out); SET_VECTOR_ELT(result, 5, v);
        memcpy(REAL(v), o_max, n_out * sizeof(double));
        v = allocVector(REALSXP, n_out); SET_VECTOR_ELT(result, 6, v);
        memcpy(REAL(v), o_mean, n_out * sizeof(double));
        v = allocVector(INTSXP, n_out); SET_VECTOR_ELT(result, 7, v);
        memcpy(INTEGER(v), o_n, n_out * sizeof(int));
        v = allocVector(INTSXP, n_out); SET_VECTOR_ELT(result, 8, v);
        memcpy(INTEGER(v), o_block, n_out * sizeof(int));
        v = allocVector(INTSXP, n_out); SET_VECTOR_ELT(result, 9, v);
        memcpy(INTEGER(v), o_smoothed, n_out * sizeof(int));
    }
    for (k = 0; k < 10; k++) SET_STRING_ELT(result_names, k, mkChar(names[k]));
    setAttrib(result, R_NamesSymbol, result_names);
    UNPROTECT(2);
    return result;
}
