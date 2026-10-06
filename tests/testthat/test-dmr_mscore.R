test_that("dmr_mscore reuses input IO and calculates permutation FDR", {
  block_starts <- seq(100, 2000, by = 100)
  cpg_positions <- as.vector(rbind(
    block_starts, block_starts + 20, block_starts + 40
  ))
  cpg_file <- tempfile(fileext = ".cpg.tsv")
  data.table::fwrite(
    data.frame(chr = "chr1", pos = cpg_positions), cpg_file,
    sep = "\t", col.names = FALSE
  )
  sample_names <- paste0("sample", seq_len(6))
  mhap_files <- vapply(seq_along(sample_names), function(sample_index) {
    mock_mhap <- data.frame(
      chr = "chr1", start = block_starts, end = block_starts + 40,
      string = if (sample_index <= 3L) "000" else "100",
      count = 10, strand = "+"
    )
    path <- tempfile(fileext = ".mhap")
    data.table::fwrite(mock_mhap, path, sep = "\t", col.names = FALSE)
    path
  }, character(1))
  names(mhap_files) <- sample_names
  on.exit(unlink(c(cpg_file, mhap_files)), add = TRUE)

  design_data <- data.frame(
    sample = sample_names,
    group = factor(rep(c("control", "case"), each = 3),
                   levels = c("control", "case"))
  )
  original_stream_reader <- getFromNamespace(".stream_text_batches", "MscoreDMR")
  original_scanner <- getFromNamespace(".scan_mhap_file", "MscoreDMR")
  mhap_read_count <- 0L
  cpg_read_count <- 0L
  smoothing_group_signatures <- character()
  smoothing_spans <- numeric()
  original_stats_table <- getFromNamespace(".smoothing_stats_table", "MscoreDMR")
  original_smooth <- getFromNamespace(".smooth_locus_stats", "MscoreDMR")
  testthat::local_mocked_bindings(
    .stream_text_batches = function(path, ...) {
      if (path %in% mhap_files) mhap_read_count <<- mhap_read_count + 1L
      else cpg_read_count <<- cpg_read_count + 1L
      original_stream_reader(path, ...)
    },
    .scan_mhap_file = function(path, ...) {
      mhap_read_count <<- mhap_read_count + 1L
      original_scanner(path, ...)
    },
    .smoothing_stats_table = function(chromosome, loci, sample_group, cov_q75) {
      smoothing_group_signatures <<- c(
        smoothing_group_signatures,
        paste(c("control", "case")[sample_group], collapse = "|")
      )
      original_stats_table(chromosome, loci, sample_group, cov_q75)
    },
    .smooth_locus_stats = function(locus_stats, threshold, max_gap_smooth,
                                   max_gap_region, min_cpgs, bp_span,
                                   min_in_span, ...) {
      smoothing_spans <<- c(smoothing_spans, min_in_span)
      original_smooth(locus_stats, threshold, max_gap_smooth, max_gap_region,
                      min_cpgs, bp_span, min_in_span, ...)
    },
    .package = "MscoreDMR"
  )

  set.seed(2026)
  expect_message(
    result <- dmr_mscore(
      mhap_files, cpg_file, design_data, ~ group, "groupcase",
      B = 10, threshold = 0.05, bp_span = 1200, min_in_span = 20
    ),
    "reducing B to 9"
  )

  expect_s4_class(result, "GRanges")
  expect_gt(length(result), 0L)
  expect_equal(mhap_read_count, 6L)
  expect_equal(cpg_read_count, 1L)
  expect_equal(length(smoothing_group_signatures), 10L)
  expect_equal(length(unique(smoothing_group_signatures)), 10L)
  expect_equal(smoothing_spans, rep(20, 10))
  expect_equal(S4Vectors::metadata(result)$candidate_discovery$bp_span, 1200)
  expect_equal(
    smoothing_group_signatures[1L],
    paste(as.character(design_data$group), collapse = "|")
  )
  expect_true(all(c(
    "wald_statistic", "empirical_p_value", "fdr", "null_pool_size"
  ) %in% names(S4Vectors::mcols(result))))
  expect_true(all(is.finite(S4Vectors::mcols(result)$wald_statistic)))
  expect_true(all(S4Vectors::mcols(result)$empirical_p_value > 0 &
                  S4Vectors::mcols(result)$empirical_p_value <= 1))
  expect_true(all(S4Vectors::mcols(result)$fdr > 0 &
                  S4Vectors::mcols(result)$fdr <= 1))
  expect_gte(S4Vectors::mcols(result)$null_pool_size[1], 1L)
  qc <- S4Vectors::metadata(result)$qc
  expect_equal(qc$total_mhap_records, 120L)
  expect_equal(qc$invalid_mhap_records, 0L)
  expect_equal(qc$nonredundant_permutations_available, 9)
  expect_equal(qc$permutations_requested, 10L)
  expect_equal(qc$permutations_generated, 9L)
  expect_lte(qc$permutations_successful, 9L)

  # Same public workflow on bgzipped, indexed inputs; the compact read cache
  # is used either way and must give identical results.
  indexed_files <- vapply(mhap_files, Rsamtools::bgzip, character(1))
  on.exit(unlink(c(indexed_files, paste0(indexed_files, ".tbi"))), add = TRUE)
  for (path in indexed_files) Rsamtools::indexTabix(path, seq = 1L,
    start = 2L, end = 3L, zeroBased = FALSE)
  original_prepare <- getFromNamespace(".prepare_mscore_store", "MscoreDMR")
  testthat::local_mocked_bindings(.prepare_mscore_store = function(...) {
    store <- original_prepare(...)
    expect_identical(store$read_backend, "compact")
    expect_length(store$read_paths, 6L)
    expect_length(store$read_sources, 0L)
    expect_length(store$read_chunks, 0L)
    store
  }, .package = "MscoreDMR")
  set.seed(2026)
  indexed_result <- dmr_mscore(indexed_files, cpg_file, design_data, ~ group,
    "groupcase", B = 10, threshold = 0.05, bp_span = 1200, min_in_span = 20)
  expect_equal(as.data.frame(indexed_result), as.data.frame(result), tolerance = 0)
  expect_equal(S4Vectors::metadata(indexed_result)$qc, qc)
})

test_that("parallel passes give results identical to a sequential run", {
  skip_on_os("windows")
  extdata <- system.file("extdata", package = "MscoreDMR")
  samples <- c(paste0("control", 1:3), paste0("case", 1:3))
  mhap_files <- stats::setNames(
    file.path(extdata, paste0(samples, ".mhap.gz")), samples
  )
  design_data <- data.frame(sample = samples,
    group = factor(rep(c("control", "case"), each = 3),
                   levels = c("control", "case")))
  run <- function(workers) {
    set.seed(7)
    messages <- character()
    result <- withCallingHandlers(
      dmr_mscore(mhap_files, file.path(extdata, "toy_cpg.tsv"), design_data,
        ~ group, "groupcase", B = 5, threshold = 0.05, bp_span = 1200,
        min_in_span = 20, workers = workers),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    list(result = result, messages = messages)
  }
  sequential <- run(1L)
  parallel <- run(2L)
  expect_identical(parallel$result, sequential$result)
  expect_gt(length(sequential$result), 0L)
  # Worker messages are replayed in the parent, one per permutation.
  expect_equal(sum(grepl("^Permutation [0-9]+ of 5", parallel$messages)), 5L)
  expect_error(dmr_mscore(mhap_files, file.path(extdata, "toy_cpg.tsv"),
    design_data, ~ group, "groupcase", workers = 0), "workers")
})

test_that("permutations are unique, non-redundant, and reduced when necessary", {
  observed <- c("control", "control", "case", "case")

  expect_message(
    plan <- MscoreDMR:::.generate_unique_permutations(observed, B = 10),
    NA
  )
  expect_equal(plan$available, 2)
  expect_equal(nrow(plan$assignments), 2L)
  expect_equal(length(unique(apply(plan$assignments, 1, paste,
                                   collapse = "|"))), 2L)
  expect_false(any(apply(plan$assignments, 1, function(assignment) {
    all(unname(assignment) == unname(observed))
  })))
  reversed <- ifelse(observed == "control", "case", "control")
  expect_false(any(apply(plan$assignments, 1, identical, reversed)))
  expect_true(all(apply(plan$assignments, 1, function(x) {
    all(as.integer(table(x)) == c(2L, 2L))
  })))

})

test_that("unbalanced permutations exclude original-condition-pure small groups", {
  observed <- c(rep("small", 3), rep("large", 5))
  names(observed) <- paste0("sample", seq_along(observed))
  set.seed(2026)
  plan <- MscoreDMR:::.generate_unique_permutations(observed, B = 100)

  expect_equal(plan$available, 45)
  expect_equal(nrow(plan$assignments), 45L)
  expect_equal(colnames(plan$assignments), names(observed))
  expect_equal(length(unique(apply(plan$assignments, 1, paste,
                                   collapse = "|"))), 45L)
  expect_false(any(apply(plan$assignments, 1, function(assignment) {
    all(unname(assignment) == unname(observed))
  })))
  expect_true(all(apply(plan$assignments, 1, function(assignment) {
    sum(assignment == "small") == 3L &&
      sum(assignment == "large") == 5L &&
      length(unique(observed[assignment == "small"])) == 2L
  })))
})

test_that("permutation selection prioritises the lowest dmrseq similarity tier", {
  observed <- c(rep("small", 4), rep("large", 6))
  set.seed(2026)
  plan <- MscoreDMR:::.generate_unique_permutations(observed, B = 10)

  similarity <- apply(plan$assignments, 1, function(assignment) {
    max(table(observed[assignment == "small"]))
  })
  expect_equal(nrow(plan$assignments), 10L)
  expect_true(all(similarity == 2L))
})

test_that("unbalanced permutation generation is independent of level order", {
  small_first <- c(rep("small", 3), rep("large", 5))
  large_first <- c(rep("large", 5), rep("small", 3))

  set.seed(2026)
  plan_small_first <- MscoreDMR:::.generate_unique_permutations(
    small_first, B = 10
  )
  set.seed(2026)
  plan_large_first <- MscoreDMR:::.generate_unique_permutations(
    large_first, B = 10
  )

  expect_equal(plan_small_first$available, 45)
  expect_equal(plan_large_first$available, 45)
  expect_true(all(rowSums(plan_small_first$assignments == "small") == 3L))
  expect_true(all(rowSums(plan_large_first$assignments == "small") == 3L))
  expect_true(all(apply(plan_small_first$assignments, 1, function(x) {
    length(unique(small_first[x == "small"])) == 2L
  })))
  expect_true(all(apply(plan_large_first$assignments, 1, function(x) {
    length(unique(large_first[x == "small"])) == 2L
  })))
})

test_that("sorted empirical p-values equal the direct count", {
  set.seed(9)
  null <- c(rnorm(5000), 0, -2, 2, 2, NA[0])
  observed <- c(rnorm(300, sd = 2), 0, 2, -2, 10, -10)
  direct <- vapply(observed, function(statistic) {
    (sum(abs(null) >= abs(statistic)) + 1) / (length(null) + 1)
  }, numeric(1))
  expect_identical(MscoreDMR:::.empirical_p_values(observed, null), direct)
  expect_identical(MscoreDMR:::.empirical_p_values(c(1, -3), numeric()), c(1, 1))
})
