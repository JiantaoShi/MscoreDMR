#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_region_nsq(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                         SEXP, SEXP);
extern SEXP C_scan_mhap(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_local_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_smoothing_stats(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_region_sums(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_refine_segments(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_trim_shape(SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"C_region_nsq", (DL_FUNC) &C_region_nsq, 10},
    {"C_scan_mhap", (DL_FUNC) &C_scan_mhap, 4},
    {"C_local_fit", (DL_FUNC) &C_local_fit, 7},
    {"C_smoothing_stats", (DL_FUNC) &C_smoothing_stats, 5},
    {"C_region_sums", (DL_FUNC) &C_region_sums, 5},
    {"C_refine_segments", (DL_FUNC) &C_refine_segments, 9},
    {"C_trim_shape", (DL_FUNC) &C_trim_shape, 2},
    {NULL, NULL, 0}
};

void R_init_MscoreDMR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
