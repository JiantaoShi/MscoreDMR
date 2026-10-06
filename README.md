# MscoreDMR

MscoreDMR finds differentially methylated regions (DMRs) genome-wide from
read-level methylation haplotype ([mHap](https://jiantaoshi.github.io/mHap/))
files. It needs only the mHap files and a CpG position index:

1. **CpG alignment.** Every character of each mHap methylation string is
   mapped to an absolute CpG coordinate from the CpG index. Each read
   contributes to every CpG it covers with its global state (methylated if
   any of its CpGs is methylated), giving per-sample M-scores.
2. **Candidate discovery.** The between-group M-score difference is smoothed
   by weighted local quadratic regression and segmented into candidate
   regions, following the two-group candidate discovery of
   [dmrseq](https://bioconductor.org/packages/dmrseq/).
3. **Testing.** Each candidate is tested with a two-stage Beta-Bernoulli
   generalized least-squares (GLS) model for any design formula and
   coefficient.
4. **Significance.** Wald statistics are calibrated against a pooled null
   built from group-label permutations, giving empirical p-values and
   Benjamini-Hochberg FDR.

The mHap files are streamed once by a compiled reader, the heavy steps run in
C, and permutations can run in parallel, so a six-sample whole-genome analysis
takes about 10 minutes on one server node (see [Performance](#performance)).

## Installation

MscoreDMR contains C code, so a C compiler and zlib are required (Linux: gcc
and zlib headers; macOS: Xcode Command Line Tools, `xcode-select --install`).

```r
install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors", "Rsamtools"))
install.packages(c("data.table", "locfit", "matrixStats"))

# From a source checkout or a built tarball
install.packages("path/to/MscoreDMR", repos = NULL, type = "source")
```

Tested with R 4.4 (Linux, gcc) and R 4.5 (macOS arm64, clang).

## Input

Two kinds of input are needed: one mHap file per sample and one CpG index
that matches the genome build of the mHap files. Files may be plain text,
gzip or bgzip compressed. Tabix indexes (`.tbi`) are not needed by
`dmr_mscore()`; they are only used if you query the files yourself.

### mHap files

An mHap file has six tab-separated, headerless columns: chromosome, start,
end, methylation string (`0`/`1` per CpG), read count and strand. `start` and
`end` are the 1-based positions of the first and last CpG of the read.

The example uses a public esophageal squamous-cell carcinoma (ESCC) dataset
with three tumor and three normal samples (hg19):

| Sample | Group | mHap | Index |
|---|---|---|---|
| SRX8208812 | Tumor | [SRX8208812.mhap.gz](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208812.mhap.gz) | [.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208812.mhap.gz.tbi) |
| SRX8208813 | Tumor | [SRX8208813.mhap.gz](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208813.mhap.gz) | [.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208813.mhap.gz.tbi) |
| SRX8208814 | Tumor | [SRX8208814.mhap.gz](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208814.mhap.gz) | [.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208814.mhap.gz.tbi) |
| SRX8208802 | Normal | [SRX8208802.mhap.gz](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208802.mhap.gz) | [.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208802.mhap.gz.tbi) |
| SRX8208803 | Normal | [SRX8208803.mhap.gz](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208803.mhap.gz) | [.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208803.mhap.gz.tbi) |
| SRX8208804 | Normal | [SRX8208804.mhap.gz](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208804.mhap.gz) | [.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208804.mhap.gz.tbi) |

Each file is about 400 MB.

### CpG index

A headerless file listing every CpG of the genome; the first two columns are
the chromosome and the 1-based CpG position (further columns are ignored):

- [hg19_CpG.gz](http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz)
  ([hg19_CpG.gz.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz.tbi))

For other genomes, use the
[annotation file](https://jiantaoshi.github.io/mHap/AnnotationFiles.html) that
matches the coordinates of your mHap files.

### Downloading the example data

```bash
mkdir -p ESCC && cd ESCC
base=http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public
for s in SRX8208812 SRX8208813 SRX8208814; do wget -c "$base/tumor/$s.mhap.gz"; done
for s in SRX8208802 SRX8208803 SRX8208804; do wget -c "$base/normal/$s.mhap.gz"; done
wget -c http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz
```

## Quick start on the bundled toy data

```r
library(MscoreDMR)

extdata <- system.file("extdata", package = "MscoreDMR")
samples <- c(paste0("control", 1:3), paste0("case", 1:3))
mhap_files <- setNames(file.path(extdata, paste0(samples, ".mhap.gz")), samples)
design_data <- data.frame(
  sample = samples,
  group = factor(rep(c("control", "case"), each = 3), levels = c("control", "case"))
)
set.seed(1)
dmrs <- dmr_mscore(mhap_files, file.path(extdata, "toy_cpg.tsv"), design_data,
                   formula = ~ group, coef = "groupcase", B = 2,
                   threshold = 0.05, bp_span = 1200, min_in_span = 20)
dmrs
```

## Whole-genome analysis of the ESCC example

```r
library(MscoreDMR)

samples <- c("SRX8208812", "SRX8208813", "SRX8208814",   # tumor
             "SRX8208802", "SRX8208803", "SRX8208804")   # normal
mhap_files <- setNames(file.path("ESCC", paste0(samples, ".mhap.gz")), samples)

design_data <- data.frame(
  sample = samples,
  group  = factor(rep(c("Tumor", "Normal"), each = 3), levels = c("Normal", "Tumor"))
)

set.seed(20260902)
dmrs <- dmr_mscore(
  mhap_files  = mhap_files,
  cpg_files   = "ESCC/hg19_CpG.gz",
  design_data = design_data,
  formula     = ~ group,
  coef        = "groupTumor",
  B           = 9,            # group-label permutations
  workers     = 5,            # parallel passes (forked processes; 1 on Windows)
  cache_dir   = "tmp"         # scratch space for the compact read cache
)
saveRDS(dmrs, "ESCC_dmrs.rds")
```

### Choosing `coef`

The design matrix comes from `model.matrix(formula, design_data)`. With
`group` a factor whose levels are `c("Normal", "Tumor")`, Normal is the
reference and the Tumor-versus-Normal coefficient is `groupTumor`. A positive
Wald statistic then means a higher M-score in Tumor. `coef` can also be a
coefficient index or a numeric contrast vector, and the formula may contain
covariates, provided it has one two-level grouping variable (named `group`
or `groups`, or uniquely identifiable from `coef`) used for the permutations.

### Main parameters

| Argument | Default | Meaning |
|---|---|---|
| `B` | 10 | Number of group-label permutations for the null (reduced to the number of distinct non-redundant relabellings, e.g. 9 for 3 vs 3). |
| `threshold` | 0.10 | Minimum absolute smoothed M-score difference of a candidate CpG. |
| `min_cpgs` | 5 | Minimum number of CpGs per candidate region. |
| `max_gap_smooth` / `max_gap_region` | 2500 / 1000 | Maximum CpG gap (bp) within a smoothing block / a candidate region. |
| `bp_span`, `min_in_span` | 1000, 30 | Bandwidth of the local regression: at least `bp_span` bp and about `min_in_span` CpGs. |
| `smoother` | `"direct"` | `"direct"` fits the local regression exactly at every CpG (compiled); `"locfit"` uses locfit's interpolated fit as in dmrseq. |
| `workers` | 1 | Number of forked worker processes for the observed and permutation passes. Results do not depend on it. |
| `chunk_nrows` | 100000 | Lines per batch when reading the CpG index (and unsorted mHap files, which use the R reader). |
| `cache_dir` | `tempdir()` | Parent directory of the private scratch cache, removed when the call returns. |

## Output

`dmr_mscore()` returns a `GRanges` with one row per tested candidate region
(1-based, inclusive coordinates) and these metadata columns:

| Column | Meaning |
|---|---|
| `direction` | Sign of the smoothed difference, first minus second level of the grouping factor (`"+"`: higher in the first level). |
| `max_abs_smooth`, `mean_smooth` | Maximum absolute and mean smoothed M-score difference in the region. |
| `n_loci` | Number of CpGs in the region. |
| `block`, `n_smoothed` | Smoothing block and number of smoothed CpGs (diagnostics). |
| `wald_statistic` | GLS Wald statistic of `coef`. |
| `empirical_p_value` | Permutation p-value, (number of pooled null statistics at least as extreme + 1) / (null size + 1). |
| `fdr` | Benjamini-Hochberg adjusted `empirical_p_value`. |
| `null_pool_size` | Number of null statistics pooled from the permutations. |

`S4Vectors::metadata(dmrs)$qc` records input and filtering counts (records
read, invalid records by cause, CpGs screened, regions dropped, permutations
used), and `metadata(dmrs)$candidate_discovery` the discovery parameters.

```r
library(GenomicRanges)
sig <- dmrs[dmrs$fdr <= 0.05]
table(ifelse(sig$wald_statistic > 0, "Hyper in tumor", "Hypo in tumor"))
write.csv(as.data.frame(dmrs), "ESCC_dmrs.csv", row.names = FALSE)
```

## Checking the ESCC result against reference regions

`ESCC_DMR_subset.txt` (shipped in `inst/extdata`) lists 150 regions known to
be hypermethylated, hypomethylated or unchanged (`NC`) in tumor. The
`Category` column is used only for this check.

```r
library(data.table)
ref <- fread(system.file("extdata", "ESCC_DMR_subset.txt", package = "MscoreDMR"))
ref_gr <- GRanges(ref$Chr, IRanges(ref$Start + 1, ref$End))   # BED-style starts
hits <- findOverlaps(ref_gr, dmrs)
best <- data.table(ref = queryHits(hits), fdr = dmrs$fdr[subjectHits(hits)],
                   wald = dmrs$wald_statistic[subjectHits(hits)])[order(fdr), .SD[1], by = ref]
call <- rep("NC", nrow(ref))
significant <- best[fdr <= 0.05]
call[significant$ref] <- ifelse(significant$wald > 0, "Hyper", "Hypo")
table(reference = ref$Category, MscoreDMR = call)
```

## Performance

Whole-genome analysis of the six ESCC samples (530,722,975 mHap records,
27.6 million covered CpGs, hg19), B = 9 permutations, on one Linux server
node (Intel Xeon Gold 6126, R 4.4, gcc 4.8):

| Implementation | Workers | Wall time | CPU time | Peak memory |
|---|---|---|---|---|
| Original R implementation (`run02`) | 1 | 82,614 s (22.9 h) | – | – |
| MscoreDMR 0.1.0 | 5 | **573 s (9.6 min)** | 1,116 s | 10.1 GB |

Peak memory is the largest resident set of a single process; the
per-chromosome CpG count matrices are kept in memory and shared by the forked
workers. The compact read cache in `cache_dir` needs about 2 GB of disk per
sample and is removed when the call returns. A single chromosome (chr21, 1/73
of the genome) takes 18 s with one worker.

### Results compared with the original implementation

The current version and the original R implementation read identical input
(530,722,975 valid records; 27,522,427 CpGs pass the coverage screen) and fit
the same GLS model: for the 108,045 candidates with identical boundaries the
Wald statistics agree to 1.6e-14. Candidate boundaries differ slightly
because the default `smoother = "direct"` fits the local regression exactly at
every CpG, whereas locfit (used by the original) interpolates its fit.

| | Original (`run02`) | MscoreDMR 0.1.0 |
|---|---|---|
| Tested candidate regions | 282,695 | 306,653 |
| Regions with FDR ≤ 0.05 | 122,737 (595.7 Mb) | 131,233 (584.5 Mb) |
| Median candidate width | 1,280 bp | 1,179 bp |
| Pooled null statistics | 1,845,876 | 2,013,323 |
| ESCC reference regions called correctly | 147 / 150 | 148 / 150 |

Agreement between the two runs: 93.1% of the significant base pairs of the
original are significant in MscoreDMR and 94.8% vice versa (base-pair Jaccard
index 0.886); overlapping candidates have Wald statistics correlated at 0.973
with the same sign in 99.97% of cases.

## Other functions

- `find_mscore_candidates()`: candidate discovery only, for two groups.
- `fit_mscore_gls()`: the two-stage GLS model for summary matrices.

## License

MIT. Portions of the candidate-discovery code are adapted from dmrseq (MIT
License).
