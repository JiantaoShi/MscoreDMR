#include <R.h>
#include <Rinternals.h>
/* No fused multiply-add contraction (clang contracts a*b + c by default on
 * arm64), so that results are identical across compilers and CPUs. */
#ifdef __clang__
#pragma clang fp contract(off)
#endif

/*
 * Sum of count * n^2 over reads overlapping each region, where n is the number
 * of CpGs shared by the read and the region (CpG index ranges [a, b] and
 * [region_a, region_b]).
 *
 * Reads must be sorted by physical start `lo`. Reads that physically overlap
 * a region but share no CpG with it are counted in `empty`, weighted by
 * `n_distinct` (the number of distinct mHap records collapsed into the row).
 *
 * Returns list(N = numeric(length(region_start)), empty = numeric(1)).
 */
SEXP C_region_nsq(SEXP lo, SEXP hi, SEXP a, SEXP b, SEXP count,
                  SEXP n_distinct, SEXP region_start, SEXP region_end,
                  SEXP region_a, SEXP region_b)
{
    R_xlen_t n_reads, n_regions, i, r, left, right, mid;
    const double *p_lo, *p_hi, *p_count, *p_nd, *p_rs, *p_re;
    const int *p_a, *p_b, *p_ra, *p_rb;
    double *run_max_hi, *p_out, empty = 0.0;
    SEXP out, N, empty_sexp, names;

    if (TYPEOF(lo) != REALSXP || TYPEOF(hi) != REALSXP ||
        TYPEOF(count) != REALSXP || TYPEOF(n_distinct) != REALSXP ||
        TYPEOF(region_start) != REALSXP || TYPEOF(region_end) != REALSXP)
        error("Coordinates, counts and record multiplicities must be double.");
    if (TYPEOF(a) != INTSXP || TYPEOF(b) != INTSXP ||
        TYPEOF(region_a) != INTSXP || TYPEOF(region_b) != INTSXP)
        error("CpG indexes must be integer.");
    n_reads = XLENGTH(lo);
    n_regions = XLENGTH(region_start);
    if (XLENGTH(hi) != n_reads || XLENGTH(a) != n_reads ||
        XLENGTH(b) != n_reads || XLENGTH(count) != n_reads ||
        XLENGTH(n_distinct) != n_reads)
        error("Read vectors must have equal lengths.");
    if (XLENGTH(region_end) != n_regions || XLENGTH(region_a) != n_regions ||
        XLENGTH(region_b) != n_regions)
        error("Region vectors must have equal lengths.");

    p_lo = REAL(lo); p_hi = REAL(hi); p_count = REAL(count);
    p_nd = REAL(n_distinct); p_a = INTEGER(a); p_b = INTEGER(b);
    p_rs = REAL(region_start); p_re = REAL(region_end);
    p_ra = INTEGER(region_a); p_rb = INTEGER(region_b);

    /* Running maximum of read ends: the first read that can reach a region
     * start is found by binary search on this non-decreasing array. */
    run_max_hi = (double *) R_alloc(n_reads > 0 ? n_reads : 1, sizeof(double));
    for (i = 0; i < n_reads; i++) {
        if (i > 0 && p_lo[i] < p_lo[i - 1])
            error("Reads must be sorted by start coordinate.");
        run_max_hi[i] = (i > 0 && run_max_hi[i - 1] > p_hi[i]) ?
            run_max_hi[i - 1] : p_hi[i];
    }

    PROTECT(N = allocVector(REALSXP, n_regions));
    p_out = REAL(N);
    for (r = 0; r < n_regions; r++) {
        double start = p_rs[r], end = p_re[r], total = 0.0;
        int ra = p_ra[r], rb = p_rb[r];
        left = 0;
        right = n_reads;
        while (left < right) {
            mid = left + (right - left) / 2;
            if (run_max_hi[mid] < start) left = mid + 1; else right = mid;
        }
        for (i = left; i < n_reads && p_lo[i] <= end; i++) {
            int first, last;
            double shared;
            if (p_hi[i] < start) continue;
            first = p_a[i] > ra ? p_a[i] : ra;
            last = p_b[i] < rb ? p_b[i] : rb;
            shared = (double) last - (double) first + 1.0;
            if (shared > 0) total += p_count[i] * shared * shared;
            else empty += p_nd[i];
        }
        p_out[r] = total;
        if ((r & 1023) == 0) R_CheckUserInterrupt();
    }

    PROTECT(empty_sexp = ScalarReal(empty));
    PROTECT(out = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, N);
    SET_VECTOR_ELT(out, 1, empty_sexp);
    PROTECT(names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("N"));
    SET_STRING_ELT(names, 1, mkChar("empty"));
    setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(4);
    return out;
}
