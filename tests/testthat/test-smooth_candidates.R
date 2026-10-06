test_that("smoothing uses true CpGs and does not bridge a 2500 bp gap", {
  positions <- c(seq(100, 1000, by = 100), seq(4000, 4900, by = 100))
  cpg_file <- tempfile(fileext = ".cpg.tsv")
  data.table::fwrite(
    data.frame(chr = "chr1", pos = positions), cpg_file,
    sep = "\t", col.names = FALSE
  )
  make_mhap <- function(methylated) {
    data.frame(
      chr = "chr1", start = positions, end = positions,
      string = if (methylated) "1" else "0",
      count = 10, strand = "+"
    )
  }
  mhap_files <- vapply(seq_len(4), function(index) {
    path <- tempfile(fileext = ".mhap")
    data.table::fwrite(make_mhap(index <= 2), path,
                       sep = "\t", col.names = FALSE)
    path
  }, character(1))
  names(mhap_files) <- paste0("sample", seq_along(mhap_files))
  on.exit(unlink(c(cpg_file, mhap_files)), add = TRUE)

  candidates <- find_mscore_candidates(
    mhap_files, cpg_file,
    groups = c("group1", "group1", "group2", "group2"),
    threshold = 0.10
  )

  expect_s4_class(candidates, "GRanges")
  expect_length(candidates, 2L)
  expect_equal(IRanges::start(candidates), c(100L, 4000L))
  expect_equal(IRanges::end(candidates), c(1000L, 4900L))
  expect_equal(S4Vectors::mcols(candidates)$direction, c("+", "+"))
  expect_equal(S4Vectors::metadata(candidates)$coordinate_system,
               "1-based inclusive")
  expect_equal(S4Vectors::metadata(candidates)$qc$invalid_mhap_records, 0L)
})

test_that("region gap splitting is independent of the smoothing gap", {
  smoothed <- data.table::data.table(
    chr = "chr1", pos = c(100, 200, 1700, 1800), block = 1L,
    beta_raw_std = rep(0.20, 4), smooth_beta = rep(0.20, 4)
  )

  regions <- MscoreDMR:::.refine_candidate_segments(
    smoothed, threshold = 0.10, max_gap_region = 1000, min_cpgs = 2
  )

  expect_equal(regions$region_start, c(100, 1700))
  expect_equal(regions$region_end, c(200, 1800))
  expect_equal(regions$n_loci, c(2L, 2L))
})

test_that("raw-effect edge trimming recalculates candidate summaries", {
  smoothed <- data.table::data.table(
    chr = "chr1", pos = seq(100, 500, by = 100), block = 1L,
    beta_raw_std = c(-0.05, 0.20, -0.01, 0.30, -0.02),
    smooth_beta = c(0.11, 0.12, 0.13, 0.30, 0.14)
  )

  regions <- MscoreDMR:::.refine_candidate_segments(
    smoothed, threshold = 0.10, max_gap_region = 1000, min_cpgs = 3
  )

  expect_equal(nrow(regions), 1L)
  expect_equal(regions$region_start, 200)
  expect_equal(regions$region_end, 400)
  expect_equal(regions$n_loci, 3L)
  expect_equal(regions$max_abs_smooth, 0.30)
  expect_equal(regions$mean_smooth, mean(c(0.12, 0.13, 0.30)))

  retained <- MscoreDMR:::.refine_candidate_segments(
    smoothed, threshold = 0.10, max_gap_region = 1000, min_cpgs = 4
  )
  expect_equal(retained$n_loci, 5L)
})

test_that("candidate parameter validation separates smoothing and region gaps", {
  expect_no_error(MscoreDMR:::.validate_smoothing_parameters(0.1, 1000, 1500, 5))
  expect_error(
    MscoreDMR:::.validate_smoothing_parameters(0.1, 2500, 1000, 0),
    "at least 3"
  )
})

test_that("candidate inputs use pooled effects and dmrseq-style weights", {
  samples <- c("g1a", "g1b", "g2a", "g2b")
  groups <- stats::setNames(c("g1", "g1", "g2", "g2"), samples)
  values <- data.table::data.table(
    chr = "chr1",
    pos = rep(c(100, 200, 300, 400, 500), each = 4L),
    sample = rep(samples, times = 5L),
    S = c(
      1, 0, 0, 0,
      10, 0, 0, 0,
      0, 0, 5, 0,
      1, 0, 2, 8,
      10, 10, 0, 0
    ),
    T = c(
      1, 9, 1, 9,
      10, 0, 10, 0,
      0, 0, 10, 0,
      1, 0, 10, 10,
      10, 10, 10, 10
    )
  )

  stats <- MscoreDMR:::.compute_candidate_smoothing_stats(
    values, groups, c("g1", "g2")
  )

  first <- stats[pos == 100]
  expect_equal(first$group1_pooled, 0.1)
  expect_equal(first$group2_pooled, 0)
  expect_equal(first$beta, 0.1)
  expect_false(isTRUE(all.equal(first$beta, mean(c(1, 0)))))
  expect_equal(first$mad_group1, 1.4826 * 0.5)
  expect_equal(first$mad_group2, 0)
  expect_equal(first$sd_raw, 1.4826^2 * 0.5)
  expect_equal(first$beta_raw_std, first$beta / first$sd_raw)

  both_single <- stats[pos == 200]
  expect_equal(both_single$cov_mean, 5)
  expect_equal(both_single$n_covered_group1, 1)
  expect_equal(both_single$n_covered_group2, 1)
  expect_equal(both_single$mad_group1, 1)
  expect_equal(both_single$mad_group2, 1)

  no_group1_depth <- stats[pos == 300]
  expect_true(is.na(no_group1_depth$beta))

  borrowed <- stats[pos == 400]
  expect_equal(borrowed$mad_group2, 1.4826 * 0.3)
  expect_equal(borrowed$mad_group1, borrowed$mad_group2)

  bounded <- stats[pos == 500]
  expect_equal(bounded$sd_raw, 1e-5)
  expect_equal(attr(stats, "cov_q75"), 5.25)
  expected_weights <- pmin(stats$cov_mean, attr(stats, "cov_q75")) /
    pmax(stats$sd_raw, 1 / pmax(stats$cov_mean, 5))
  finite_beta <- is.finite(stats$beta)
  expect_equal(stats$weight[finite_beta], expected_weights[finite_beta])
  expect_true(all(is.finite(stats$weight[finite_beta])))
  expect_true(all(stats$weight[finite_beta] > 0))

  reversed <- MscoreDMR:::.compute_candidate_smoothing_stats(
    values, groups, c("g2", "g1")
  )
  expect_equal(reversed$beta, -stats$beta)
})

test_that("adaptive locfit neighborhood is a bounded fraction", {
  expect_equal(
    MscoreDMR:::.adaptive_locfit_nn(
      positions = seq(100, 1000, by = 100),
      density_span = 3000, min_in_span = 30, max_nn = 0.75
    ),
    0.75
  )
  expect_equal(
    MscoreDMR:::.adaptive_locfit_nn(
      positions = seq(100, 10000, length.out = 100),
      density_span = 1000, min_in_span = 30, max_nn = 0.75
    ),
    1000 / 9901
  )
})

test_that("subthreshold same-sign loci split peaks instead of bridging them", {
  input <- data.table::data.table(
    chr = "chr1", pos = seq_len(12) * 100, block = 1L,
    smooth_beta = c(rep(0.2, 5), 0.02, 0.03, rep(0.2, 5)),
    beta_raw_std = 1
  )
  regions <- MscoreDMR:::.refine_candidate_segments(input, 0.1, 1000, 5)
  expect_equal(regions$n_loci, c(5L, 5L))
  expect_equal(regions$region_start, c(100, 800))
  expect_equal(regions$region_end, c(500, 1200))
  input[, `:=`(smooth_beta = -smooth_beta, beta_raw_std = -beta_raw_std)]
  negative <- MscoreDMR:::.refine_candidate_segments(input, 0.1, 1000, 5)
  expect_equal(negative$region_start, regions$region_start)
  expect_equal(negative$region_end, regions$region_end)
  expect_true(all(negative$direction == -1))
})

test_that("cutoff equality, absent raw support and minimum-sized regions are retained", {
  input <- data.table::data.table(
    chr = "chr1", pos = seq_len(5) * 100, block = 1L,
    smooth_beta = 0.1, beta_raw_std = -1
  )
  regions <- MscoreDMR:::.refine_candidate_segments(input, 0.1, 1000, 5)
  expect_equal(regions$n_loci, 5L)
  input <- data.table::rbindlist(list(input, input[5][, pos := 600]))
  expect_equal(MscoreDMR:::.refine_candidate_segments(input, 0.1, 1000, 5)$n_loci, 6L)
  expect_equal(MscoreDMR:::.candidate_direction(c(0.1, -0.1, 0, NA), 0.1),
               c(1L, -1L, 0L, NA_integer_))
})

test_that("shape refinement follows the pinned peak and shoulder rule", {
  x <- c(seq(0.11, 0.9, length.out = 10), seq(0.8, 0.11, length.out = 9))
  x <- x + rep(c(0, 0.001, -0.001), length.out = length(x))
  # Hand-calculated from the pinned source's cut and peak-relative guards.
  expect_equal(MscoreDMR:::.trim_candidate_shape(x, 5), 5:15)
  expect_equal(MscoreDMR:::.trim_candidate_shape(rep(0.2, 10), 5), 1:10)
  expect_equal(MscoreDMR:::.trim_candidate_shape(x, 18), seq_along(x))
})

test_that("small-block fallback uses standardized signal and min_cpgs", {
  groups <- stats::setNames(c("a", "a", "b", "b"), paste0("s", 1:4))
  values <- data.table::CJ(pos = c(100, 200, 400, 500), sample = names(groups))
  values[, `:=`(chr = "chr1", T = 10, S = ifelse(sample %in% c("s1", "s2"), 8, 2))]
  result <- MscoreDMR:::.run_smoothing_pipeline(
    values, groups, c("a", "b"), 0.1, max_gap_smooth = 100,
    max_gap_region = 1000, min_cpgs = 3
  )
  expect_length(result, 1)
  expect_equal(S4Vectors::mcols(result)$n_smoothed, 0L)
  expect_equal(S4Vectors::mcols(result)$mean_smooth, 0.6 / 1e-5)
})

test_that("observed coverage screening keeps a fixed locus universe without mutating input", {
  groups <- c(a1 = "a", a2 = "a", b1 = "b", b2 = "b")
  input <- data.table::data.table(chr = "chr1", pos = c(100, 100, 200),
                                  sample = c("a1", "b1", "a2"), S = 1, T = 2)
  original <- data.table::copy(input)
  screened <- MscoreDMR:::.screen_observed_candidate_loci(input, groups, c("a", "b"))
  expect_equal(unique(screened$loci$pos), 100)
  expect_equal(screened$qc$candidate_loci_dropped_group_coverage, 1L)
  expect_equal(input, original)
  # a1 and b1 can land in the same permuted group; the selected coordinates
  # remain fixed even though the new effect is missing.
  permuted <- c(a1 = "a", a2 = "b", b1 = "a", b2 = "b")
  expect_length(MscoreDMR:::.run_smoothing_pipeline(
    screened$loci, permuted, c("a", "b"), 0.1, check_group_coverage = FALSE
  ), 0L)
  empty <- MscoreDMR:::.screen_observed_candidate_loci(input[3], groups, c("a", "b"))
  expect_equal(nrow(empty$loci), 0L)
  expect_equal(names(empty$loci), names(input))
})

test_that("observed missing group coverage fails explicitly", {
  groups <- c(a1 = "a", a2 = "a", b1 = "b", b2 = "b")
  values <- data.table::data.table(chr = "chr1", pos = 100, sample = "a1", S = 1, T = 2)
  expect_error(MscoreDMR:::.run_smoothing_pipeline(values, groups, c("a", "b"), 0.1),
               "without pooled coverage")
  expect_length(MscoreDMR:::.run_smoothing_pipeline(
    values, groups, c("a", "b"), 0.1, check_group_coverage = FALSE
  ), 0)
})

test_that("permutation smoothing fits finite rows and predicts all block coordinates", {
  groups <- c(a1 = "a", a2 = "a", b1 = "b", b2 = "b")
  positions <- seq_len(20) * 100
  values <- data.table::CJ(pos = positions, sample = names(groups))
  values[, `:=`(chr = "chr1", T = 10,
               S = ifelse(sample %in% c("a1", "a2"), 8, 2))]
  # Missing coverage in one group at an internal locus of a permutation.
  values <- values[!(pos == 1000 & sample %in% c("b1", "b2"))]
  original_locfit <- locfit::locfit.raw
  fitted_positions <- NULL
  testthat::local_mocked_bindings(
    locfit.raw = function(x, y, weights, ...) {
      expect_true(all(is.finite(y)))
      expect_true(all(is.finite(weights) & weights > 0))
      fitted_positions <<- x
      original_locfit(x, y, weights = weights, ...)
    }, .package = "locfit"
  )
  expect_warning(result <- MscoreDMR:::.run_smoothing_pipeline(
    values, groups, c("a", "b"), 0.1, check_group_coverage = FALSE,
    smoother = "locfit"
  ), NA)
  expect_equal(fitted_positions, positions[positions != 1000])
  expect_length(result, 1L)
  expect_equal(S4Vectors::mcols(result)$n_loci, 20L)
  expect_equal(S4Vectors::mcols(result)$n_smoothed, 20L)
  expect_equal(as.numeric(GenomicRanges::start(result)), 100)
  expect_equal(as.numeric(GenomicRanges::end(result)), 2000)
})

test_that("locfit failures stop with genomic context rather than raw fallback", {
  groups <- c(a1 = "a", a2 = "a", b1 = "b", b2 = "b")
  values <- data.table::CJ(pos = seq_len(10) * 100, sample = names(groups))
  values[, `:=`(chr = "chr1", T = 10, S = ifelse(sample %in% c("a1", "a2"), 8, 2))]
  testthat::local_mocked_bindings(locfit.raw = function(...) stop("mock fitting failure"),
                                 .package = "locfit")
  expect_error(MscoreDMR:::.run_smoothing_pipeline(values, groups, c("a", "b"), 0.1,
                                                   smoother = "locfit"),
               "chr1 / block 1.*mock fitting failure")
})

test_that("locfit.raw with predict equals the formula fit with preplot", {
  set.seed(3)
  posi <- sort(sample(1:50000, 400))
  yi <- sin(posi / 3000) + rnorm(400, sd = 0.2)
  wi <- runif(400, 0.5, 5)
  lp <- locfit::lp
  for (nn in c(0.05, 0.3, 0.75)) {
    formula_fit <- locfit::locfit(yi ~ lp(posi, nn = nn, h = 1000, deg = 2),
      data = data.frame(posi, yi, wi), weights = wi, family = "gaussian",
      kern = "tricube", maxk = 10000)
    pp <- stats::preplot(formula_fit, where = "data", band = "local",
                         newdata = data.frame(posi = posi))
    raw_fit <- locfit::locfit.raw(posi, yi, weights = wi, alpha = c(nn, 1000),
      deg = 2, kern = "tricube", family = "gaussian", maxk = 10000)
    expect_identical(as.numeric(stats::predict(raw_fit, newdata = posi)),
                     as.numeric(pp$trans(pp$fit)))
  }
})

test_that("the direct C fit equals locfit evaluated exactly at the data", {
  set.seed(11)
  for (case in 1:12) {
    n <- sample(c(5, 40, 400), 1)
    x <- sort(sample(1:30000, n)) + 0
    y <- sin(x / 2500) + rnorm(n, sd = 0.2)
    w <- runif(n, 0.2, 5)
    # nn = 30 / n reproduces the rounding case k = n * nn = 30.
    nn <- min(0.75, if (case %% 2) 30 / n else runif(1, 0.05, 0.75))
    h <- sample(c(0, 500, 1000), 1)
    lp <- locfit::lp
    data <- data.frame(x = x, y = y, w = w)
    exact <- suppressWarnings(locfit::locfit(y ~ lp(x, nn = nn, h = h, deg = 2),
      data = data, weights = w, kern = "tricube", maxk = 10000,
      ev = locfit::dat()))
    direct <- .Call(MscoreDMR:::C_local_fit, x, y, w, x, nn, h, 2L)
    expect_equal(direct, as.numeric(stats::fitted(exact)), tolerance = 1e-6)
  }
  expect_error(.Call(MscoreDMR:::C_local_fit, c(2, 1), c(0, 0), c(1, 1), 1, 0.5, 0, 2L),
               "sorted")
  expect_error(.Call(MscoreDMR:::C_local_fit, 1, 0, 0, 1, 0.5, 0, 2L), "positive")
})

test_that("direct and locfit smoothers give close candidate scans", {
  groups <- c(a1 = "a", a2 = "a", b1 = "b", b2 = "b")
  set.seed(5)
  values <- data.table::CJ(pos = sort(sample(1:40000, 300)), sample = names(groups))
  values[, `:=`(chr = "chr1", T = 20,
    S = round(20 * plogis(ifelse(sample %in% c("a1", "a2") & pos > 15000 & pos < 25000,
                                 2, 0) + rnorm(.N, sd = 0.3))))]
  direct <- MscoreDMR:::.run_smoothing_pipeline(values, groups, c("a", "b"), 0.1)
  locfit_result <- MscoreDMR:::.run_smoothing_pipeline(values, groups, c("a", "b"), 0.1,
                                                       smoother = "locfit")
  expect_identical(S4Vectors::metadata(direct)$candidate_discovery$smoother, "direct")
  expect_gt(length(direct), 0L)
  expect_equal(length(direct), length(locfit_result))
  # Differences stem from locfit's tree interpolation only.
  expect_lt(max(abs(S4Vectors::mcols(direct)$max_abs_smooth -
                    S4Vectors::mcols(locfit_result)$max_abs_smooth)), 0.1)
})

test_that("compiled smoothing statistics equal the R reference", {
  set.seed(21)
  for (case in 1:15) {
    n_samples <- sample(c(4, 5, 6, 8), 1)
    samples <- paste0("s", seq_len(n_samples))
    groups <- stats::setNames(sample(rep(c("A", "B"), length.out = n_samples)), samples)
    loci <- data.table::CJ(chr = c("chr1", "chr2"), pos = sort(sample(1:5000, 60)),
                           sample = samples)
    # Missing samples, zero-coverage rows and loci uncovered in one group.
    loci <- loci[runif(.N) > 0.25]
    loci[, T := sample(c(0, 1, 3, 10, 40), .N, replace = TRUE)]
    loci[, S := round(T * runif(.N))]
    loci <- loci[sample(.N)]
    for (q75 in list(NULL, 7.5)) {
      reference <- MscoreDMR:::.compute_candidate_smoothing_stats_r(
        loci, groups, c("A", "B"), cov_q75 = q75, min_in_span = 20)
      compiled <- MscoreDMR:::.compute_candidate_smoothing_stats(
        loci, groups, c("A", "B"), cov_q75 = q75, min_in_span = 20)
      expect_equal(compiled, reference, tolerance = 0)
    }
  }
})

test_that("compiled candidate refinement and trimming equal the R reference", {
  set.seed(22)
  for (case in 1:40) {
    n <- sample(c(10, 80, 400), 1)
    pos <- sort(sample(1:(n * 150), n))
    chr <- rep(c("chr1", "chr10", "chr2"), length.out = n)[order(runif(n))]
    wave <- sin(seq_len(n) / sample(3:12, 1)) * runif(1, 0.05, 0.5)
    smoothed <- data.table::data.table(chr = chr, pos = as.numeric(pos),
      smooth_beta = wave + rnorm(n, sd = 0.02),
      beta_raw_std = wave * 4 + rnorm(n, sd = 0.3),
      block = sample(1:4, n, replace = TRUE), smoothed = runif(n) > 0.2)
    smoothed[sample(n, n %/% 10), smooth_beta := NA_real_]
    smoothed[sample(n, n %/% 10), beta_raw_std := NA_real_]
    threshold <- sample(c(0, 0.05, 0.1, 0.2), 1)
    min_cpgs <- sample(3:6, 1)
    expect_equal(
      MscoreDMR:::.refine_candidate_segments(smoothed, threshold, 1000, min_cpgs),
      MscoreDMR:::.refine_candidate_segments_r(smoothed, threshold, 1000, min_cpgs),
      tolerance = 0)
  }
  for (case in 1:300) {
    n <- sample(4:60, 1)
    peak <- sample(n, 1)
    x <- 0.1 + exp(-((seq_len(n) - peak) / runif(1, 1, n / 2))^2) * runif(1, 0, 2) +
      rnorm(n, sd = runif(1, 0, 0.1))
    min_cpgs <- sample(3:8, 1)
    expect_identical(.Call(MscoreDMR:::C_trim_shape, x, min_cpgs),
                     MscoreDMR:::.trim_candidate_shape(x, min_cpgs))
  }
})
