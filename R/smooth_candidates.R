#' Construct candidate regions from CpG-aligned mHap files
#'
#' Maps every methylation-string character to an absolute CpG coordinate from
#' `cpg_files`, aggregates read-level global methylation states at the true CpG
#' loci, and smooths the between-group M-score difference in independent
#' genomic blocks. mHap `start` and `end` locate covered CpGs only; they are
#' never used as proxy CpG loci.
#'
#' @param mhap_files Named character vector or list of mHap file paths. Each
#'   file has six headerless columns: chromosome, start, end, methylation
#'   string, count, and strand.
#' @param cpg_files One shared CpG annotation path, or a named vector/list with
#'   one path per sample. Each headerless annotation has chromosome and 1-based
#'   CpG position in its first two columns. Bgzip/tabix and plain TSV files are
#'   supported.
#' @param groups Vector assigning each sample to one of exactly two groups.
#' @param group_levels Optional length-two vector specifying the order used for
#'   `group_levels[1] - group_levels[2]`.
#' @param threshold Per-locus absolute candidate cutoff, also used on the
#'   standardized raw signal for edge refinement. Defaults to `0.10`.
#' @param max_gap_smooth Maximum distance between adjacent true CpG loci in one
#'   smoothing block. Defaults to `2500` bp.
#' @param max_gap_region Maximum distance between adjacent true CpG loci in one
#'   candidate region. Defaults to `1000` bp, independently of smoothing gaps.
#' @param min_cpgs Minimum number of valid true-CpG loci retained in a
#'   candidate region after two edge refinements. Defaults to `5`, at least `3`.
#' @param bp_span Local regression distance bandwidth in bp. Defaults to `1000`.
#' @param min_in_span Adaptive nearest-neighbor span parameter. Defaults to `30`.
#' @param chunk_nrows Maximum mHap records per input batch. Defaults to `100000`.
#' @param cache_dir Parent directory for a private temporary disk cache. Defaults
#'   to [tempdir()]. The private cache is removed when the call returns.
#' @param smoother Local regression used for candidate smoothing. `"direct"`
#'   (default) fits the weighted local quadratic regression exactly at every
#'   CpG in compiled code. `"locfit"` uses [locfit::locfit.raw()], which fits
#'   at the vertices of an adaptive tree and interpolates between them (the
#'   dmrseq behaviour); both use the same tricube kernel and bandwidth
#'   `max(bp_span, distance to the k-th nearest CpG)`, so results differ only
#'   by locfit's interpolation error.
#' @details Candidate discovery follows the ordinary two-group, local-DMR path
#'   of dmrseq commit bd5187a759c362927344a6e049b4924fa0d5368a, using M-scores.
#'   Loci without positive pooled coverage in both observed groups are screened
#'   out before candidate discovery, with counts recorded in QC. Missing
#'   sample coverage is not a zero M-score. Fitting failures stop with context.
#'   No maximum candidate width or CpG count is imposed.
#' @return A [GenomicRanges::GRanges] object with candidate-level metadata.
#'   Global ingestion and alignment QC is stored in `metadata(result)$qc`.
#' @export
#'
#' @examples
#' extdata <- system.file("extdata", package = "MscoreDMR")
#' samples <- c(paste0("control", 1:3), paste0("case", 1:3))
#' mhap_files <- stats::setNames(
#'   file.path(extdata, paste0(samples, ".mhap.gz")), samples
#' )
#' candidates <- find_mscore_candidates(
#'   mhap_files, file.path(extdata, "toy_cpg.tsv"),
#'   groups = rep(c("control", "case"), each = 3),
#'   group_levels = c("case", "control"),
#'   threshold = 0.05, bp_span = 1200, min_in_span = 20
#' )
#' candidates
find_mscore_candidates <- function(mhap_files, cpg_files, groups,
                                   group_levels = NULL, threshold = 0.10,
                                   max_gap_smooth = 2500,
                                   max_gap_region = 1000,
                                   min_cpgs = 5, bp_span = 1000,
                                   min_in_span = 30, chunk_nrows = 100000L,
                                   cache_dir = NULL,
                                   smoother = c("direct", "locfit")) {
  smoother <- match.arg(smoother)
  prepared <- .prepare_mscore_store(mhap_files, cpg_files, chunk_nrows,
                                     cache_dir, keep_reads = FALSE)
  on.exit(unlink(prepared$root, recursive = TRUE), add = TRUE)
  sample_names <- names(prepared$mhap_files)
  group_info <- .normalise_groups(groups, sample_names, group_levels)
  .validate_smoothing_parameters(
    threshold, max_gap_smooth, max_gap_region, min_cpgs, bp_span, min_in_span
  )
  prepared <- .screen_mscore_store(prepared, group_info$groups,
                                   group_info$levels, min_in_span)
  candidates <- .smooth_mscore_store(
    prepared, group_info$groups, group_info$levels,
    threshold, max_gap_smooth, max_gap_region, min_cpgs, bp_span, min_in_span,
    smoother = smoother
  )
  S4Vectors::metadata(candidates)$qc <- prepared$qc
  S4Vectors::metadata(candidates)$cpg_annotation_files <-
    unique(unname(prepared$cpg_files))
  S4Vectors::metadata(candidates)$coordinate_system <- "1-based inclusive"
  candidates
}

.screen_observed_candidate_loci <- function(loci, groups, group_levels) {
  # Screen once under observed labels. Never subset the read table or CpG cache.
  coverage <- loci[, .(eligible =
    sum(T[sample %in% names(groups)[groups == group_levels[1L]]], na.rm = TRUE) > 0 &
    sum(T[sample %in% names(groups)[groups == group_levels[2L]]], na.rm = TRUE) > 0
  ), by = c("chr", "pos")]
  keep <- coverage[which(coverage$eligible), c("chr", "pos"), with = FALSE]
  dropped <- nrow(coverage) - nrow(keep)
  message("Candidate coverage screen: retained ", nrow(keep), " of ",
          nrow(coverage), " CpGs; excluded ", dropped,
          " without positive pooled coverage in both observed groups.")
  list(loci = loci[keep, on = c("chr", "pos"), nomatch = 0L],
       qc = list(candidate_loci_before_screen = nrow(coverage),
                 candidate_loci_retained = nrow(keep),
                 candidate_loci_dropped_group_coverage = dropped))
}

.normalise_groups <- function(groups, sample_names, group_levels = NULL) {
  if (length(groups) != length(sample_names) || anyNA(groups)) {
    stop("`groups` must contain one non-missing value per sample.",
         call. = FALSE)
  }
  if (!is.null(names(groups)) && all(nzchar(names(groups)))) {
    if (anyDuplicated(names(groups)) ||
        !setequal(names(groups), sample_names)) {
      stop("Named `groups` must match the mHap sample names.",
           call. = FALSE)
    }
    groups <- groups[match(sample_names, names(groups))]
  }
  groups <- as.character(groups)
  observed <- unique(groups)
  if (length(observed) != 2L) {
    stop("`groups` must identify exactly two groups.", call. = FALSE)
  }
  if (is.null(group_levels)) {
    group_levels <- observed
  } else {
    group_levels <- as.character(group_levels)
    if (length(group_levels) != 2L || anyNA(group_levels) ||
        !setequal(group_levels, observed)) {
      stop("`group_levels` must contain the two values present in `groups`.",
           call. = FALSE)
    }
  }
  list(groups = stats::setNames(groups, sample_names), levels = group_levels)
}

.validate_smoothing_parameters <- function(threshold, max_gap_smooth,
                                           max_gap_region, min_cpgs,
                                           bp_span = 1000, min_in_span = 30) {
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("`threshold` must be one finite, non-negative number.",
         call. = FALSE)
  }
  if (!is.numeric(max_gap_smooth) || length(max_gap_smooth) != 1L ||
      !is.finite(max_gap_smooth) || max_gap_smooth < 0) {
    stop("`max_gap_smooth` must be one finite, non-negative number.",
         call. = FALSE)
  }
  if (!is.numeric(max_gap_region) || length(max_gap_region) != 1L ||
      !is.finite(max_gap_region) || max_gap_region < 0) {
    stop("`max_gap_region` must be one finite, non-negative number.",
         call. = FALSE)
  }
  if (!is.numeric(min_cpgs) || length(min_cpgs) != 1L ||
      !is.finite(min_cpgs) || min_cpgs < 3 ||
      min_cpgs > .Machine$integer.max || min_cpgs != floor(min_cpgs)) {
    stop("`min_cpgs` must be one integer at least 3.", call. = FALSE)
  }
  for (parameter in c("bp_span", "min_in_span")) {
    value <- get(parameter)
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop("`", parameter, "` must be one finite positive number.", call. = FALSE)
    }
  }
}

.compute_candidate_smoothing_stats <- function(master_locus_dt, groups,
                                                group_levels,
                                                cov_q75 = NULL,
                                                min_in_span = 30L) {
  sample_names <- names(groups)
  if (is.null(sample_names) || anyDuplicated(sample_names) ||
      anyNA(sample_names) || any(!nzchar(sample_names))) {
    stop("`groups` must be uniquely named by sample.", call. = FALSE)
  }
  if (length(group_levels) != 2L ||
      !setequal(unname(groups), group_levels)) {
    stop("`group_levels` must contain the two values present in `groups`.",
         call. = FALSE)
  }
  required_columns <- c("chr", "pos", "sample", "S", "T")
  missing_columns <- setdiff(required_columns, names(master_locus_dt))
  if (length(missing_columns)) {
    stop("The master locus table is missing: ",
         paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(master_locus_dt) == 0L) return(data.table::data.table())
  sample_index <- match(master_locus_dt$sample, sample_names)
  if (anyNA(sample_index)) {
    stop("Locus-table samples must be present in `groups`.", call. = FALSE)
  }
  S_values <- as.numeric(master_locus_dt$S)
  T_values <- as.numeric(master_locus_dt$T)
  if (any(!is.finite(S_values)) || any(!is.finite(T_values)) ||
      any(S_values < 0) || any(T_values < 0) || any(S_values > T_values)) {
    stop("The master locus table must satisfy 0 <= S <= T with finite values.",
         call. = FALSE)
  }
  n_samples <- length(sample_names)
  if (is.null(cov_q75)) {
    coverage <- master_locus_dt[, .(
      cov_mean = sum(as.numeric(T)) / n_samples
    ), by = .(chr, pos)]
    positive_coverage <- coverage$cov_mean[
      is.finite(coverage$cov_mean) & coverage$cov_mean > 0
    ]
    if (!length(positive_coverage)) {
      stop("No finite positive CpG coverage was available for smoothing.",
           call. = FALSE)
    }
    cov_q75 <- as.numeric(stats::quantile(
      positive_coverage, probs = 0.75, na.rm = TRUE,
      names = FALSE, type = 7
    ))
  }
  if (length(cov_q75) != 1L || !is.finite(cov_q75) || cov_q75 <= 0) {
    stop("`cov_q75` must be one finite, positive number.", call. = FALSE)
  }
  sample_group <- match(unname(groups[sample_names]), group_levels)
  dense <- .dense_loci_from_long(master_locus_dt, sample_names)
  results <- lapply(names(dense), function(chromosome) {
    .smoothing_stats_table(chromosome, dense[[chromosome]], sample_group, cov_q75)
  })
  result <- data.table::rbindlist(results, use.names = TRUE)
  data.table::setorder(result, chr, pos)
  chromosome_density <- result[, .(
    density_span = as.numeric(min_in_span) *
      (max(pos) - min(pos) + 1) / .N
  ), by = chr]
  data.table::setattr(result, "cov_q75", cov_q75)
  data.table::setattr(result, "density_span", mean(chromosome_density$density_span))
  result
}

# Locus statistics of one chromosome from dense S/T matrices (compiled).
.smoothing_stats_table <- function(chromosome, loci, sample_group, cov_q75) {
  columns <- .Call(C_smoothing_stats, as.numeric(loci$pos), loci$S, loci$T,
                   as.integer(sample_group), as.numeric(cov_q75))
  if (columns$invalid_weight > 0) {
    stop("Internal error: finite beta produced an invalid smoothing weight.",
         call. = FALSE)
  }
  columns$invalid_weight <- NULL
  data.table::as.data.table(c(list(chr = rep(chromosome, length(columns$pos))), columns))
}

# R reference implementation of .compute_candidate_smoothing_stats(), kept
# for regression tests of the compiled version.
.compute_candidate_smoothing_stats_r <- function(master_locus_dt, groups,
                                                  group_levels,
                                                  cov_q75 = NULL,
                                                  min_in_span = 30L) {
  sample_names <- names(groups)
  if (is.null(sample_names) || anyDuplicated(sample_names) ||
      anyNA(sample_names) || any(!nzchar(sample_names))) {
    stop("`groups` must be uniquely named by sample.", call. = FALSE)
  }
  if (length(group_levels) != 2L ||
      !setequal(unname(groups), group_levels)) {
    stop("`group_levels` must contain the two values present in `groups`.",
         call. = FALSE)
  }
  required_columns <- c("chr", "pos", "sample", "S", "T")
  missing_columns <- setdiff(required_columns, names(master_locus_dt))
  if (length(missing_columns)) {
    stop("The master locus table is missing: ",
         paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(master_locus_dt) == 0L) return(data.table::data.table())
  if (any(!master_locus_dt$sample %in% sample_names)) {
    stop("Locus-table samples must be present in `groups`.", call. = FALSE)
  }
  S_values <- master_locus_dt$S
  T_values <- master_locus_dt$T
  if (any(!is.finite(S_values)) || any(!is.finite(T_values)) ||
      any(S_values < 0) || any(T_values < 0) || any(S_values > T_values)) {
    stop("The master locus table must satisfy 0 <= S <= T with finite values.",
         call. = FALSE)
  }

  n_samples <- length(sample_names)
  coverage <- master_locus_dt[, .(
    cov_mean = sum(as.numeric(T)) / n_samples
  ), by = .(chr, pos)]
  if (is.null(cov_q75)) {
    positive_coverage <- coverage$cov_mean[
      is.finite(coverage$cov_mean) & coverage$cov_mean > 0
    ]
    if (!length(positive_coverage)) {
      stop("No finite positive CpG coverage was available for smoothing.",
           call. = FALSE)
    }
    cov_q75 <- as.numeric(stats::quantile(
      positive_coverage, probs = 0.75, na.rm = TRUE,
      names = FALSE, type = 7
    ))
  }
  if (length(cov_q75) != 1L || !is.finite(cov_q75) || cov_q75 <= 0) {
    stop("`cov_q75` must be one finite, positive number.", call. = FALSE)
  }

  group1_columns <- which(unname(groups[sample_names]) == group_levels[1L])
  group2_columns <- which(unname(groups[sample_names]) == group_levels[2L])
  chromosome_results <- vector("list", length(unique(coverage$chr)))
  chromosomes <- unique(coverage$chr)
  for (chromosome_index in seq_along(chromosomes)) {
    chromosome <- chromosomes[chromosome_index]
    chromosome_loci <- coverage[chr == chromosome]
    data.table::setorder(chromosome_loci, pos)
    chromosome_values <- master_locus_dt[chr == chromosome]
    row_index <- match(chromosome_values$pos, chromosome_loci$pos)
    column_index <- match(chromosome_values$sample, sample_names)
    if (anyNA(row_index) || anyNA(column_index)) {
      stop("Internal error while aligning chromosome matrices.",
           call. = FALSE)
    }
    matrix_names <- list(as.character(chromosome_loci$pos), sample_names)
    matrix_dimensions <- c(nrow(chromosome_loci), n_samples)
    S_matrix <- matrix(0, matrix_dimensions[1L], matrix_dimensions[2L],
                       dimnames = matrix_names)
    T_matrix <- matrix(0, matrix_dimensions[1L], matrix_dimensions[2L],
                       dimnames = matrix_names)
    matrix_index <- cbind(row_index, column_index)
    S_matrix[matrix_index] <- as.numeric(chromosome_values$S)
    T_matrix[matrix_index] <- as.numeric(chromosome_values$T)

    prop_matrix <- S_matrix / T_matrix
    prop_matrix[T_matrix <= 0] <- NA_real_
    S_group1 <- matrixStats::rowSums2(
      S_matrix, cols = group1_columns, na.rm = TRUE
    )
    T_group1 <- matrixStats::rowSums2(
      T_matrix, cols = group1_columns, na.rm = TRUE
    )
    S_group2 <- matrixStats::rowSums2(
      S_matrix, cols = group2_columns, na.rm = TRUE
    )
    T_group2 <- matrixStats::rowSums2(
      T_matrix, cols = group2_columns, na.rm = TRUE
    )
    group1_pooled <- S_group1 / T_group1
    group2_pooled <- S_group2 / T_group2
    group1_pooled[T_group1 <= 0] <- NA_real_
    group2_pooled[T_group2 <= 0] <- NA_real_

    mad_group1 <- matrixStats::rowMads(
      prop_matrix, cols = group1_columns, na.rm = TRUE, constant = 1.4826
    )
    mad_group2 <- matrixStats::rowMads(
      prop_matrix, cols = group2_columns, na.rm = TRUE, constant = 1.4826
    )
    n_covered_group1 <- matrixStats::rowSums2(
      T_matrix[, group1_columns, drop = FALSE] > 0
    )
    n_covered_group2 <- matrixStats::rowSums2(
      T_matrix[, group2_columns, drop = FALSE] > 0
    )

    both_single <- n_covered_group1 == 1L & n_covered_group2 == 1L
    borrow_group2 <- n_covered_group1 == 1L & n_covered_group2 >= 2L
    borrow_group1 <- n_covered_group2 == 1L & n_covered_group1 >= 2L
    mad_group1[both_single] <- 1
    mad_group2[both_single] <- 1
    mad_group1[borrow_group2] <- mad_group2[borrow_group2]
    mad_group2[borrow_group1] <- mad_group1[borrow_group1]

    beta <- group1_pooled - group2_pooled
    # Intentional source parity: rowMads scaling AND the outer factor occur
    # in the pinned dmrseq implementation. Do not silently remove either.
    sd_raw <- 1.4826 * sqrt(mad_group1^2 + mad_group2^2)
    sd_raw[is.finite(sd_raw) & sd_raw < 1e-5] <- 1e-5
    cov_mean <- matrixStats::rowSums2(T_matrix, na.rm = TRUE) / n_samples
    weight <- pmin(cov_mean, cov_q75) /
      pmax(sd_raw, 1 / pmax(cov_mean, 5))
    keep <- is.finite(beta)
    if (any(!is.finite(weight[keep]) | weight[keep] <= 0)) {
      stop("Internal error: finite beta produced an invalid smoothing weight.",
           call. = FALSE)
    }
    chromosome_results[[chromosome_index]] <- data.table::data.table(
      chr = chromosome,
      pos = chromosome_loci$pos,
      group1_pooled = group1_pooled,
      group2_pooled = group2_pooled,
      beta = beta,
      beta_raw_std = beta / (sd_raw * 2 / sqrt(n_samples)),
      mad_group1 = mad_group1,
      mad_group2 = mad_group2,
      sd_raw = sd_raw,
      cov_mean = cov_mean,
      weight = weight,
      n_covered_group1 = n_covered_group1,
      n_covered_group2 = n_covered_group2
    )
  }
  result <- data.table::rbindlist(chromosome_results, use.names = TRUE)
  data.table::setorder(result, chr, pos)
  chromosome_density <- result[, .(
    density_span = as.numeric(min_in_span) *
      (max(pos) - min(pos) + 1) / .N
  ), by = chr]
  data.table::setattr(result, "cov_q75", cov_q75)
  data.table::setattr(result, "density_span", mean(chromosome_density$density_span))
  result
}

.adaptive_locfit_nn <- function(positions, density_span,
                                min_in_span = 30L, max_nn = 0.75) {
  positions <- as.numeric(positions)
  if (!length(positions) || !is.finite(density_span) || density_span <= 0) {
    stop("Cannot calculate an adaptive locfit neighborhood.", call. = FALSE)
  }
  physical_range <- max(positions) - min(positions) + 1
  nn <- min(
    density_span / physical_range,
    as.numeric(min_in_span) / length(positions),
    max_nn
  )
  if (!is.finite(nn) || nn <= 0 || nn > max_nn) {
    stop("Adaptive locfit `nn` must be in (0, 0.75].", call. = FALSE)
  }
  nn
}

.run_smoothing_pipeline <- function(master_locus_dt, groups, group_levels,
                                    threshold, max_gap_smooth = 2500,
                                    max_gap_region = 1000, min_cpgs = 5,
                                    bp_span = 1000, min_in_span = 30,
                                    check_group_coverage = TRUE,
                                    cov_q75 = NULL, density_span = NULL,
                                    smoother = "direct") {
  smoother <- match.arg(smoother, c("direct", "locfit"))
  .validate_smoothing_parameters(
    threshold, max_gap_smooth, max_gap_region, min_cpgs, bp_span, min_in_span
  )
  sample_names <- names(groups)
  if (is.null(sample_names) || anyDuplicated(sample_names)) {
    stop("`groups` must be named by sample.", call. = FALSE)
  }
  if (nrow(master_locus_dt) == 0L) {
    message("No valid CpG-aligned loci were available for smoothing.")
    return(GenomicRanges::GRanges())
  }
  if (any(!unique(master_locus_dt$sample) %in% sample_names)) {
    stop("Locus-table samples must be present in `groups`.", call. = FALSE)
  }
  locus_stats <- .compute_candidate_smoothing_stats(
    master_locus_dt, groups, group_levels, cov_q75 = cov_q75,
    min_in_span = min_in_span
  )
  if (is.null(density_span)) density_span <- attr(locus_stats, "density_span")
  .smooth_locus_stats(locus_stats, threshold, max_gap_smooth, max_gap_region,
                      min_cpgs, bp_span, min_in_span, check_group_coverage,
                      density_span, smoother)
}

# Smoothing, segmentation and refinement of a locus statistics table.
.smooth_locus_stats <- function(locus_stats, threshold, max_gap_smooth,
                                max_gap_region, min_cpgs, bp_span, min_in_span,
                                check_group_coverage, density_span, smoother) {
  if (check_group_coverage && any(!is.finite(locus_stats$beta))) {
    stop("Observed candidate input has ", sum(!is.finite(locus_stats$beta)),
         " locus/loci without pooled coverage in both groups. Filter these ",
         "CpGs before analysis; do not replace missing M-scores with zero.",
         call. = FALSE)
  }
  # Keep the master coordinate universe for density and gap calculations,
  # including loci with a missing pooled effect in a permutation.
  data.table::setorder(locus_stats, chr, pos)
  locus_stats[, gap := pos - data.table::shift(pos), by = chr]
  locus_stats[, block := cumsum(
    is.na(gap) | gap > max_gap_smooth
  ), by = chr]
  # Blocks are contiguous row ranges of the sorted table.
  n_loci <- nrow(locus_stats)
  chromosome_values <- locus_stats$chr
  block_values <- locus_stats$block
  positions <- as.numeric(locus_stats$pos)
  beta_values <- locus_stats$beta
  weight_values <- locus_stats$weight
  smooth_beta <- locus_stats$beta_raw_std
  smoothed_flag <- logical(n_loci)
  starts <- c(1L, which(chromosome_values[-1L] != chromosome_values[-n_loci] |
                          block_values[-1L] != block_values[-n_loci]) + 1L)
  ends <- c(starts[-1L] - 1L, n_loci)
  short_blocks <- 0L
  for (b in seq_along(starts)) {
    rows <- seq.int(starts[b], ends[b])
    chromosome <- chromosome_values[starts[b]]
    if (b == 1L || chromosome != chromosome_values[starts[b - 1L]]) {
      message("Processing ", chromosome, "...")
    }
    finite <- is.finite(beta_values[rows])
    if (length(rows) < min_cpgs || sum(finite) < min_cpgs) {
      short_blocks <- short_blocks + 1L
    } else {
      block_positions <- positions[rows]
      adaptive_nn <- .adaptive_locfit_nn(
        block_positions, density_span = density_span,
        min_in_span = min_in_span, max_nn = 0.75
      )
      smooth_beta[rows] <- tryCatch({
        # Omit missing responses only for fitting, retaining the complete
        # block for nn and prediction.
        fit_rows <- rows[finite]
        if (any(!is.finite(positions[fit_rows])) ||
            any(!is.finite(weight_values[fit_rows]) | weight_values[fit_rows] <= 0)) {
          stop("Non-finite coordinates or invalid fitting weights.")
        }
        values <- if (smoother == "direct") {
          # Exact local quadratic fit at every CpG (positions are sorted).
          .Call(C_local_fit, positions[fit_rows], as.numeric(beta_values[fit_rows]),
                as.numeric(weight_values[fit_rows]), block_positions,
                as.numeric(adaptive_nn), as.numeric(bp_span), 2L)
        } else {
          # locfit.raw() with alpha = c(nn, h) is the fit behind
          # locfit(y ~ lp(x, nn, h, deg = 2)); predict() gives the values of
          # preplot() without the unused variance band.
          fit <- locfit::locfit.raw(
            positions[fit_rows], beta_values[fit_rows],
            weights = weight_values[fit_rows],
            alpha = c(adaptive_nn, bp_span), deg = 2, kern = "tricube",
            family = "gaussian", maxk = 10000
          )
          as.numeric(stats::predict(fit, newdata = block_positions))
        }
        if (length(values) != length(rows) || any(!is.finite(values))) {
          stop("Prediction length mismatch or non-finite predictions.")
        }
        values
      },
        error = function(error) {
          stop("Candidate smoothing failed in ", chromosome, " / block ",
               block_values[starts[b]], " [", min(block_positions), ", ",
               max(block_positions), "]: ", conditionMessage(error),
               call. = FALSE)
        }
      )
      smoothed_flag[rows] <- TRUE
    }
    if (b == length(starts) || chromosome != chromosome_values[starts[b + 1L]]) {
      if (short_blocks) {
        message("  ", short_blocks, " block(s) on ", chromosome,
                " have fewer than ", min_cpgs,
                " usable CpGs; using standardized raw differences.")
      }
      short_blocks <- 0L
    }
  }
  locus_stats[, `:=`(smooth_beta = smooth_beta, smoothed = smoothed_flag)]
  smoothed <- locus_stats
  regions <- .refine_candidate_segments(
    smoothed, threshold, max_gap_region, min_cpgs
  )
  if (!nrow(regions)) {
    message("No candidate regions passed refinement and filtering.")
    return(GenomicRanges::GRanges())
  }
  data.table::setorder(regions, chr, region_start, region_end)
  result <- GenomicRanges::GRanges(
    seqnames = regions$chr,
    ranges = IRanges::IRanges(regions$region_start, regions$region_end)
  )
  S4Vectors::mcols(result)$direction <- ifelse(regions$direction > 0, "+", "-")
  S4Vectors::mcols(result)$max_abs_smooth <- regions$max_abs_smooth
  S4Vectors::mcols(result)$mean_smooth <- regions$mean_smooth
  S4Vectors::mcols(result)$n_loci <- regions$n_loci
  S4Vectors::mcols(result)$block <- regions$block
  S4Vectors::mcols(result)$n_smoothed <- regions$n_smoothed
  S4Vectors::metadata(result)$candidate_discovery <- list(
    reference_commit = "bd5187a759c362927344a6e049b4924fa0d5368a",
    smoother = smoother,
    threshold = threshold, bp_span = bp_span, min_in_span = min_in_span,
    min_cpgs = min_cpgs, max_gap_smooth = max_gap_smooth,
    max_gap_region = max_gap_region
  )
  result
}

.refine_candidate_segments <- function(smoothed, threshold,
                                       max_gap_region, min_cpgs) {
  if (!nrow(smoothed)) return(data.table::data.table())
  if (!"beta_raw_std" %in% names(smoothed)) {
    stop("Candidate refinement requires `beta_raw_std`.", call. = FALSE)
  }
  loci <- data.table::data.table(
    chr = smoothed$chr, pos = as.numeric(smoothed$pos),
    smooth_beta = as.numeric(smoothed$smooth_beta),
    beta_raw_std = as.numeric(smoothed$beta_raw_std),
    block = if ("block" %in% names(smoothed)) as.integer(smoothed$block) else NA_integer_,
    smoothed = if ("smoothed" %in% names(smoothed)) as.logical(smoothed$smoothed) else TRUE
  )
  data.table::setorder(loci, chr, pos)
  chromosomes <- unique(loci$chr)
  regions <- .Call(C_refine_segments, match(loci$chr, chromosomes), loci$pos,
                   loci$smooth_beta, loci$beta_raw_std, loci$block, loci$smoothed,
                   as.numeric(threshold), as.numeric(max_gap_region),
                   as.integer(min_cpgs))
  regions <- data.table::as.data.table(regions)
  if (!nrow(regions)) return(data.table::data.table())
  regions[, chr := chromosomes[chr_code]]
  regions[, chr_code := NULL]
  data.table::setcolorder(regions, "chr")
  regions <- regions[n_loci >= as.integer(min_cpgs)]
  if (nrow(regions) && max(regions$n_loci) > 1000L) {
    message("Large candidates with more than 1000 CpGs detected; no size cap applied.")
  }
  data.table::setorder(regions, chr, region_start, region_end)
  regions
}

# R reference implementation of .refine_candidate_segments(), kept for
# regression tests of the compiled version.
.refine_candidate_segments_r <- function(smoothed, threshold,
                                         max_gap_region, min_cpgs) {
  if (!nrow(smoothed)) return(data.table::data.table())
  smoothed <- data.table::copy(smoothed)
  data.table::setorder(smoothed, chr, pos)
  if (!"beta_raw_std" %in% names(smoothed)) {
    stop("Candidate refinement requires `beta_raw_std`.", call. = FALSE)
  }
  if (!"smoothed" %in% names(smoothed)) smoothed[, smoothed := TRUE]
  # Region clusters use ALL master coordinates, independently of smoothing
  # clusters. Missing scan values are omitted only after these clusters exist.
  smoothed[, region_cluster := cumsum(
    is.na(data.table::shift(pos)) | pos - data.table::shift(pos) > max_gap_region
  ), by = chr]
  smoothed <- smoothed[is.finite(smooth_beta)]
  if (!nrow(smoothed)) return(data.table::data.table())
  smoothed[, direction := .candidate_direction(smooth_beta, threshold)]

  # Background values remain in the run-length encoding so two same-sign
  # peaks cannot bridge subthreshold CpGs.
  smoothed[, segment := data.table::rleid(
    chr, region_cluster, direction
  )]
  candidate_loci <- smoothed[
    is.finite(smooth_beta) & direction != 0
  ]
  if (!nrow(candidate_loci)) return(data.table::data.table())

  candidate_loci[, keep_locus := {
    retained <- seq_len(.N)
    if (.N > min_cpgs) {
      matching <- which(.candidate_direction(beta_raw_std, threshold) ==
                          direction[1L])
      if (length(matching) && max(matching) - min(matching) + 1L >= min_cpgs) {
        retained <- seq.int(min(matching), max(matching))
      }
    }
    if (length(retained) > min_cpgs) {
      retained <- retained[.trim_candidate_shape(
        smooth_beta[retained] * direction[1L], min_cpgs
      )]
    }
    seq_len(.N) %in% retained
  }, by = .(chr, segment, direction)]
  candidate_loci <- candidate_loci[keep_locus == TRUE]
  if (!nrow(candidate_loci)) return(data.table::data.table())

  regions <- candidate_loci[, .(
    region_start = min(pos),
    region_end = max(pos),
    max_abs_smooth = max(abs(smooth_beta)),
    mean_smooth = mean(smooth_beta),
    n_loci = .N,
    block = block[1L],
    n_smoothed = sum(smoothed)
  ), by = .(chr, segment, direction)]
  regions <- regions[
    n_loci >= as.integer(min_cpgs)
  ]
  if (nrow(regions) && max(regions$n_loci) > 1000L) {
    message("Large candidates with more than 1000 CpGs detected; no size cap applied.")
  }
  data.table::setorder(regions, chr, region_start, region_end)
  regions
}

# Match bumphunter's tolerant positive comparison and ordinary negative bound.
.candidate_direction <- function(x, threshold) {
  result <- rep(NA_integer_, length(x))
  valid <- is.finite(x)
  result[valid] <- as.integer(x[valid] >= threshold |
    abs(x[valid] - threshold) <= .Machine$double.eps^0.5)
  result[valid & x <= -threshold] <- -1L
  result
}

# Adapted from dmrseq trimEdges(), commit bd5187a759c362927344a6e049b4924fa0d5368a.
# x is direction-adjusted scan signal, not raw M-score or genomic distance.
.trim_candidate_shape <- function(x, min_cpgs) {
  n <- length(x)
  original <- seq_len(n)
  if (n <= min_cpgs) return(original)
  mid <- which.max(x)
  first <- 1L
  last <- n
  ratio <- max(x) / min(x)
  if (is.na(ratio) || ratio <= 4 / 3) return(original)
  cut <- (0.5 * (max(x) - min(x)) + min(x) + 0.75 * mean(x)) / 2
  significant_slope <- function(index, direction) {
    if (length(index) <= 4L) return(FALSE)
    coefficients <- summary(stats::lm(x[index] ~ index))$coefficients
    nrow(coefficients) == 2L && is.finite(coefficients[2L, 4L]) &&
      coefficients[2L, 1L] * direction > 0 && coefficients[2L, 4L] < 0.01
  }
  left <- seq_len(mid)
  right <- seq.int(mid, n)
  if (significant_slope(left, 1)) {
    first <- min(max(1, round(mid - 0.125 * n)), max(1, mid - 2),
                 min(left[x[left] >= cut]))
  }
  if (significant_slope(right, -1)) {
    last <- max(min(round(mid + 0.125 * n), n), min(mid + 2, n),
                max(right[x[right] >= cut]))
  }
  if (last - first + 1L >= min_cpgs) seq.int(first, last) else original
}
