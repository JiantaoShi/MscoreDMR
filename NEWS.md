# MscoreDMR 0.1.0

* First packaged release: `dmr_mscore()`, `find_mscore_candidates()` and
  `fit_mscore_gls()`, with streaming mHap ingestion and tabix index reuse.
* Added DESCRIPTION, NAMESPACE, LICENSE, package documentation, toy data in
  `inst/extdata` and runnable examples.

## Performance (results unchanged)

* mHap files are read by a C reader (zlib) that validates, CpG-maps and
  counts records in one pass and writes a compact per-chromosome read cache.
  Unsorted inputs fall back to the R reader.
* Region N_sq is computed in C from the read cache instead of re-querying
  reads for the observed pass and every permutation.
* Candidate smoothing calls `locfit.raw()` and `predict()` directly,
  skipping the unused variance band of `preplot()`.
* No RDS serialisation inside the pipeline: per-chromosome locus counts are
  kept in memory as dense (CpG x sample) S/T matrices, shared by forked
  workers, and region S/T sums are computed from them in C. Only the
  compact read cache is written to disk (binary, in `cache_dir`).
* Smoothing statistics (pooled methylation, MADs, weights), candidate
  segmentation, refinement and edge trimming run in C; the R versions are
  kept as test references and give identical results. Smoothing blocks are
  sliced as row ranges, and short blocks are reported once per chromosome.
* New `workers` argument of `dmr_mscore()` runs the observed and permutation
  passes in forked worker processes.

## Candidate smoothing

* New `smoother` argument of `dmr_mscore()` and `find_mscore_candidates()`.
  The default `"direct"` fits the weighted local quadratic regression exactly
  at every CpG in C (same tricube kernel and bandwidth as locfit; equal to
  locfit evaluated at the data to ~1e-8). `"locfit"` keeps the previous
  tree-interpolated locfit fit. Results differ only by locfit's
  interpolation error.
