#' MscoreDMR: differentially methylated regions from mHap data
#'
#' Detects differentially methylated regions from read-level methylation
#' haplotype (mHap) files. The main entry point is [dmr_mscore()], which runs
#' candidate discovery ([find_mscore_candidates()]), two-stage GLS testing
#' ([fit_mscore_gls()]) and permutation-based FDR estimation in one call.
#'
#' @section Example data:
#' `system.file("extdata", package = "MscoreDMR")` contains a toy data set
#' used by the examples: `toy_cpg.tsv`, a headerless CpG annotation (chromosome
#' and 1-based position) with 60 CpGs on `chr1`, and six gzipped mHap files
#' (`control1.mhap.gz` to `control3.mhap.gz` and `case1.mhap.gz` to
#' `case3.mhap.gz`). Each mHap file has 20 reads covering three CpGs; case reads
#' from 800 to 1440 bp are methylated at the first two CpGs; all other reads
#' are fully unmethylated, so the toy data contain one hypermethylated region.
#' `ESCC_DMR_subset.txt` lists 150 hg19 regions (BED-style 0-based starts) that
#' are hypermethylated, hypomethylated or unchanged (`NC`) in esophageal
#' squamous-cell carcinoma, for checking analyses of the public ESCC mHap
#' data described in the README.
#'
#' @keywords internal
#' @import data.table
#' @useDynLib MscoreDMR, .registration = TRUE
"_PACKAGE"
