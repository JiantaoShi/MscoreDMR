#' Detect DMRs with CpG-aligned M-score GLS and permutation FDR
#'
#' Runs true-CpG candidate smoothing, two-stage GLS fitting, and pooled
#' permutation-null FDR estimation. mHap files are streamed once into exact
#' locus counts. Region queries reuse original tabix indexes when present;
#' unindexed inputs use bounded private read shards. Observed and permuted
#' analyses share locus caches, exact mappings and full-read global states.
#'
#' @param mhap_files Named character vector or list of mHap file paths.
#' @param cpg_files One shared CpG annotation path, or a named vector/list with
#'   one path per mHap sample. The first two headerless columns are chromosome
#'   and 1-based CpG position.
#' @param design_data Data frame with one row per sample and variables used by
#'   `formula`.
#' @param formula One-sided design formula passed to [stats::model.matrix()].
#' @param coef Character coefficient name, one-based coefficient index, or a
#'   numeric contrast vector with one value per design column.
#' @param B Positive number of group-label permutations. Defaults to `10`.
#' @param threshold Per-locus candidate cutoff. Defaults to `0.10`.
#' @param max_gap_smooth Maximum true-CpG gap allowed within a smoothing block.
#'   Defaults to `2500` bp.
#' @param max_gap_region Maximum true-CpG gap allowed within one candidate
#'   region. Defaults to `1000` bp.
#' @param min_cpgs Minimum number of valid true-CpG loci retained per candidate
#'   after two edge refinements. Defaults to `5`, at least `3`.
#' @param bp_span Local regression distance bandwidth in bp. Defaults to `1000`.
#' @param min_in_span Adaptive nearest-neighbor span parameter. Defaults to `30`.
#' @param chunk_nrows Maximum mHap records per input batch. Defaults to `100000`.
#' @param cache_dir Parent directory for a private temporary disk cache. Defaults
#'   to [tempdir()]. The private cache is removed when the call returns.
#' @param workers Number of forked worker processes used to run the observed
#'   pass and the permutation passes concurrently. Defaults to `1`
#'   (sequential). Results are identical for any value: permutations are
#'   drawn before the passes start and the passes use no random numbers.
#'   Each worker is single-threaded and holds the data of one chromosome at a
#'   time; forking is unavailable on Windows, where `1` is used.
#' @param smoother Local regression used for candidate smoothing. `"direct"`
#'   (default) fits the weighted local quadratic regression exactly at every
#'   CpG in compiled code. `"locfit"` uses [locfit::locfit.raw()], which fits
#'   at the vertices of an adaptive tree and interpolates between them (the
#'   dmrseq behaviour); both use the same tricube kernel and bandwidth
#'   `max(bp_span, distance to the k-th nearest CpG)`, so results differ only
#'   by locfit's interpolation error.
#' @details Loci are screened once for pooled coverage in both observed groups.
#'   Permutations retain this screened coordinate universe and recompute all group-dependent
#'   signals. Smoothing failures stop; they are not replaced by raw effects.
#'
#' @return A [GenomicRanges::GRanges] of observed candidates with Wald
#'   statistics, empirical p-values, and BH FDR in `mcols()`. Global alignment
#'   and coverage QC is stored in `metadata(result)$qc`.
#' @seealso [find_mscore_candidates()] for candidate discovery only and
#'   [fit_mscore_gls()] for the region-level model.
#' @export
#'
#' @examples
#' extdata <- system.file("extdata", package = "MscoreDMR")
#' samples <- c(paste0("control", 1:3), paste0("case", 1:3))
#' mhap_files <- stats::setNames(
#'   file.path(extdata, paste0(samples, ".mhap.gz")), samples
#' )
#' design_data <- data.frame(
#'   sample = samples,
#'   group = factor(rep(c("control", "case"), each = 3),
#'                  levels = c("control", "case"))
#' )
#' \donttest{
#' set.seed(1)
#' dmrs <- dmr_mscore(
#'   mhap_files, file.path(extdata, "toy_cpg.tsv"), design_data,
#'   formula = ~ group, coef = "groupcase", B = 2,
#'   threshold = 0.05, bp_span = 1200, min_in_span = 20
#' )
#' dmrs
#' }
dmr_mscore <- function(mhap_files, cpg_files, design_data, formula, coef,
                       B = 10, threshold = 0.10,
                       max_gap_smooth = 2500, max_gap_region = 1000,
                       min_cpgs = 5, bp_span = 1000, min_in_span = 30,
                       chunk_nrows = 100000L, cache_dir = NULL,
                       workers = 1L, smoother = c("direct", "locfit")) {
  smoother <- match.arg(smoother)
  mhap_files <- .normalise_mhap_files(mhap_files)
  sample_names <- names(mhap_files)
  if (!is.data.frame(design_data) || nrow(design_data) != length(mhap_files)) {
    stop("`design_data` must have one row per mHap file.", call. = FALSE)
  }
  terms_object <- if (inherits(formula, "formula")) {
    stats::terms(formula)
  } else NULL
  if (is.null(terms_object) || attr(terms_object, "response") != 0L) {
    stop("`formula` must be a one-sided model formula.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != as.integer(B)) {
    stop("`B` must be one positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers < 1 || workers != as.integer(workers)) {
    stop("`workers` must be one positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)
  if (workers > 1L && .Platform$OS.type == "windows") {
    message("Forked workers are unavailable on Windows; using workers = 1.")
    workers <- 1L
  }
  .validate_smoothing_parameters(
    threshold, max_gap_smooth, max_gap_region, min_cpgs, bp_span, min_in_span
  )

  design_data <- .align_design_data(design_data, sample_names)
  group_column <- .find_group_column(design_data, formula, coef)
  original_group <- design_data[[group_column]]
  group_levels <- if (is.factor(original_group)) {
    levels(droplevels(original_group))
  } else unique(as.character(original_group))
  if (length(group_levels) != 2L || anyNA(original_group)) {
    stop("The inferred grouping variable must contain exactly two non-missing levels.",
         call. = FALSE)
  }
  design_data[[group_column]] <- factor(
    as.character(original_group), levels = group_levels
  )
  observed_design <- stats::model.matrix(formula, data = design_data)
  if (nrow(observed_design) <= ncol(observed_design) ||
      qr(observed_design)$rank < ncol(observed_design)) {
    stop("The observed design matrix must be full rank with positive residual degrees of freedom.",
         call. = FALSE)
  }
  contrast <- .make_mscore_contrast(coef, colnames(observed_design))

  message("Loading and CpG-aligning mHap data...")
  prepared <- .prepare_mscore_store(mhap_files, cpg_files, chunk_nrows, cache_dir)
  on.exit(unlink(prepared$root, recursive = TRUE), add = TRUE)
  qc <- prepared$qc
  observed_groups <- as.character(design_data[[group_column]])
  named_observed_groups <- stats::setNames(observed_groups, sample_names)
  prepared <- .screen_mscore_store(prepared, named_observed_groups,
                                   group_levels, min_in_span)
  qc <- prepared$qc

  smooth <- function(groups, check_group_coverage = TRUE) {
    .smooth_mscore_store(
      prepared, groups, group_levels, threshold = threshold,
      max_gap_smooth = max_gap_smooth, max_gap_region = max_gap_region,
      min_cpgs = min_cpgs, bp_span = bp_span, min_in_span = min_in_span,
      check_group_coverage = check_group_coverage, smoother = smoother
    )
  }
  # One pass is smoothing plus region matrices (and, for permutations, GLS).
  observed_pass <- function() {
    message("Finding observed candidate regions...")
    regions <- smooth(named_observed_groups)
    input <- if (length(regions)) {
      .region_matrices_store(prepared, regions, sample_names)
    }
    list(regions = regions, input = input)
  }
  permutation_pass <- function(permutation) {
    message("Permutation ", permutation, " of ", effective_B, "...")
    permuted_data <- design_data
    permuted_groups <- permutation_plan$assignments[permutation, ]
    permuted_data[[group_column]] <- factor(
      permuted_groups, levels = group_levels
    )
    permuted_design <- stats::model.matrix(formula, data = permuted_data)
    if (!identical(colnames(permuted_design), colnames(observed_design)) ||
        qr(permuted_design)$rank < ncol(permuted_design)) {
      message("  Skipping permutation with a rank-deficient design matrix.")
      return(numeric())
    }
    null_regions <- smooth(stats::setNames(permuted_groups, sample_names),
                           check_group_coverage = FALSE)
    if (!length(null_regions)) return(numeric())
    null_input <- .region_matrices_store(
      prepared, null_regions, sample_names
    )
    if (!length(null_input$regions)) return(numeric())
    null_fit <- tryCatch(
      fit_mscore_gls(
        null_input$S_matrix, null_input$T_matrix,
        null_input$N_sq_matrix, permuted_design, contrast
      ),
      error = function(error) {
        message("  Skipping failed null fit: ", conditionMessage(error))
        NULL
      }
    )
    if (is.null(null_fit)) return(numeric())
    values <- unname(null_fit$statistic)
    values[is.finite(values)]
  }
  make_plan <- function() {
    plan <- .generate_unique_permutations(observed_groups, B)
    available <- nrow(plan$assignments)
    if (available < B) {
      message(
        "Requested ", B, " permutations, but only ", available,
        " unique non-redundant permutation(s) are available; reducing B to ",
        available, "."
      )
    }
    plan
  }

  if (workers == 1L) {
    observed <- observed_pass()
  } else {
    # All passes run concurrently, so the plan is drawn first. The passes
    # themselves use no random numbers, hence the results do not change.
    permutation_plan <- make_plan()
    effective_B <- nrow(permutation_plan$assignments)
    message("Running the observed pass and ", effective_B,
            " permutation(s) on ", workers, " worker processes...")
    passes <- .run_parallel_passes(
      c(list(observed_pass),
        lapply(seq_len(effective_B), function(permutation) {
          force(permutation)
          function() permutation_pass(permutation)
        })),
      workers
    )
    observed <- passes[[1L]]
  }
  observed_regions <- observed$regions
  if (!length(observed_regions)) {
    message("No observed candidate regions were found.")
    return(.attach_mscore_metadata(
      observed_regions, qc, prepared$cpg_files
    ))
  }
  observed_input <- observed$input
  qc$reads_without_local_cpg <- observed_input$qc$reads_without_local_cpg
  qc$regions_dropped_incomplete_coverage <-
    observed_input$qc$regions_dropped_incomplete_coverage
  observed_regions <- observed_input$regions
  if (!length(observed_regions)) {
    message("No observed candidates had coverage in every sample.")
    return(.attach_mscore_metadata(
      observed_regions, qc, prepared$cpg_files
    ))
  }
  observed_fit <- fit_mscore_gls(
    observed_input$S_matrix, observed_input$T_matrix,
    observed_input$N_sq_matrix, observed_design, contrast
  )
  observed_statistics <- unname(observed_fit$statistic)

  if (workers == 1L) {
    permutation_plan <- make_plan()
    effective_B <- nrow(permutation_plan$assignments)
    null_statistics <- lapply(seq_len(effective_B), permutation_pass)
  } else {
    null_statistics <- passes[-1L]
  }
  qc$permutations_requested <- B
  qc$nonredundant_permutations_available <- permutation_plan$available
  qc$permutations_generated <- effective_B
  qc$permutations_successful <- sum(lengths(null_statistics) > 0L)

  pooled_null <- unlist(null_statistics, use.names = FALSE)
  pooled_null <- pooled_null[is.finite(pooled_null)]
  if (!length(pooled_null)) {
    message("No valid null statistics were produced; empirical p-values are one.")
  }
  empirical_p <- .empirical_p_values(observed_statistics, pooled_null)
  fdr <- stats::p.adjust(empirical_p, method = "BH")
  S4Vectors::mcols(observed_regions)$wald_statistic <- observed_statistics
  S4Vectors::mcols(observed_regions)$empirical_p_value <- empirical_p
  S4Vectors::mcols(observed_regions)$fdr <- fdr
  S4Vectors::mcols(observed_regions)$null_pool_size <- length(pooled_null)
  .attach_mscore_metadata(observed_regions, qc, prepared$cpg_files)
}

# Evaluate argument-free functions in forked workers. Messages and warnings
# are collected in each worker and replayed in task order, so callers that
# log conditions see the same messages as in a sequential run.
.run_parallel_passes <- function(tasks, workers) {
  results <- parallel::mclapply(seq_along(tasks), function(index) {
    conditions <- list()
    value <- withCallingHandlers(
      tasks[[index]](),
      message = function(condition) {
        conditions[[length(conditions) + 1L]] <<- condition
        invokeRestart("muffleMessage")
      },
      warning = function(condition) {
        conditions[[length(conditions) + 1L]] <<- condition
        invokeRestart("muffleWarning")
      }
    )
    list(value = value, conditions = conditions)
  }, mc.cores = workers, mc.preschedule = FALSE)
  lapply(results, function(result) {
    if (inherits(result, "try-error")) stop(attr(result, "condition"))
    if (is.null(result)) {
      stop("A worker process ended without a result (possibly out of memory); ",
           "retry with fewer `workers`.", call. = FALSE)
    }
    for (condition in result$conditions) {
      if (inherits(condition, "warning")) warning(condition) else message(condition)
    }
    result$value
  })
}

.generate_unique_permutations <- function(observed_groups, B,
                                          enumeration_limit = 5e5) {
  observed_names <- names(observed_groups)
  observed_groups <- as.character(observed_groups)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != as.integer(B)) {
    stop("`B` must be one positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  group_sizes <- table(observed_groups)
  group_levels <- names(group_sizes)
  if (length(group_levels) != 2L || anyNA(observed_groups)) {
    stop("Permutation generation requires exactly two non-missing groups.",
         call. = FALSE)
  }

  n_samples <- length(observed_groups)
  balanced <- length(unique(as.integer(group_sizes))) == 1L
  if (balanced) {
    # A deterministic level is sufficient because complementary balanced
    # assignments are canonicalised to the same unlabelled partition below.
    small_level <- sort(group_levels)[1L]
    large_level <- setdiff(group_levels, small_level)
  } else {
    small_level <- names(group_sizes)[which.min(group_sizes)]
    large_level <- names(group_sizes)[which.max(group_sizes)]
  }
  n_small <- unname(group_sizes[small_level])
  n_large <- n_samples - n_small
  labelled_total <- choose(n_samples, n_small)

  assignment_from_indices <- function(indices) {
    assignment <- rep(large_level, n_samples)
    assignment[indices] <- small_level
    names(assignment) <- observed_names
    assignment
  }
  assignment_key <- function(assignment) paste(assignment, collapse = "\r")
  canonicalise <- function(assignment) {
    if (!balanced) return(assignment)
    reversed <- ifelse(
      assignment == small_level, large_level, small_level
    )
    names(reversed) <- names(assignment)
    if (assignment_key(assignment) <= assignment_key(reversed)) {
      assignment
    } else {
      reversed
    }
  }
  observed_key <- assignment_key(canonicalise(observed_groups))
  empty_assignments <- function() {
    result <- matrix(character(), 0L, n_samples)
    colnames(result) <- observed_names
    result
  }
  first_level_similarity <- function(candidate) {
    max(table(observed_groups[candidate == small_level]))
  }
  select_by_similarity <- function(candidates, similarity, number) {
    if (number >= length(candidates)) return(candidates)
    selected_indices <- integer()
    for (level in sort(unique(similarity))) {
      remaining <- number - length(selected_indices)
      if (remaining <= 0L) break
      tier <- which(similarity == level)
      chosen <- if (length(tier) <= remaining) {
        tier
      } else {
        sample(tier, remaining, replace = FALSE)
      }
      selected_indices <- c(selected_indices, chosen)
    }
    candidates[selected_indices]
  }

  if (is.finite(labelled_total) && labelled_total <= enumeration_limit) {
    combinations <- utils::combn(seq_len(n_samples), n_small)
    candidates <- list()
    similarity <- numeric()
    seen <- new.env(parent = emptyenv(), hash = TRUE)
    assign(observed_key, TRUE, envir = seen)
    for (index in seq_len(ncol(combinations))) {
      small_indices <- combinations[, index]
      # Pure subsets preserve an original biological condition and are not
      # suitable null comparisons, following dmrseq's enumerated two-group
      # permutation filter.
      if (length(unique(observed_groups[small_indices])) == 1L) next
      candidate <- canonicalise(assignment_from_indices(small_indices))
      key <- assignment_key(candidate)
      if (exists(key, envir = seen, inherits = FALSE)) next
      assign(key, TRUE, envir = seen)
      candidates[[length(candidates) + 1L]] <- candidate
      # Lower values mean that the new small group more evenly mixes the two
      # original conditions. Equal-score candidates remain randomly sampled.
      similarity <- c(similarity, first_level_similarity(candidate))
    }
    available <- length(candidates)
    effective_B <- min(B, available)
    if (!effective_B) {
      return(list(assignments = empty_assignments(), available = available))
    }
    assignments <- select_by_similarity(
      candidates, similarity, effective_B
    )
  } else {
    # For a large combination space, the exact number of valid label
    # assignments is available analytically even though they are not
    # materialised. Balanced complementary assignments count only once.
    available <- if (balanced) {
      labelled_total / 2 - 1
    } else {
      labelled_total - 1 - choose(n_large, n_small)
    }
    available <- max(0, available)
    effective_B <- as.integer(min(as.double(B), available))
    if (!effective_B) {
      return(list(assignments = empty_assignments(), available = available))
    }

    # Exhaustive tier ranking is impossible here. Build a bounded oversampled
    # pool of unique, mixed assignments and apply the same first-level score as
    # an explicit approximation to the enumerated dmrseq selection rule.
    pool_target <- as.integer(min(
      available, max(20 * as.double(effective_B),
                     as.double(effective_B) + 100), 10000
    ))
    max_attempts <- max(1000, 100 * pool_target)
    candidates <- list()
    similarity <- numeric()
    seen <- new.env(parent = emptyenv(), hash = TRUE)
    assign(observed_key, TRUE, envir = seen)
    attempts <- 0L
    while (length(candidates) < pool_target && attempts < max_attempts) {
      attempts <- attempts + 1L
      small_indices <- sample.int(n_samples, n_small, replace = FALSE)
      if (length(unique(observed_groups[small_indices])) == 1L) next
      candidate <- canonicalise(assignment_from_indices(small_indices))
      key <- assignment_key(candidate)
      if (exists(key, envir = seen, inherits = FALSE)) next
      assign(key, TRUE, envir = seen)
      candidates[[length(candidates) + 1L]] <- candidate
      similarity <- c(similarity, first_level_similarity(candidate))
    }
    if (length(candidates) < effective_B) {
      warning(
        "Only ", length(candidates), " valid unique permutation(s) were ",
        "generated before the large-space attempt limit was reached.",
        call. = FALSE
      )
      effective_B <- length(candidates)
    }
    if (!effective_B) {
      return(list(assignments = empty_assignments(), available = available))
    }
    assignments <- select_by_similarity(
      candidates, similarity, effective_B
    )
  }

  assignment_matrix <- do.call(rbind, assignments)
  colnames(assignment_matrix) <- observed_names
  list(assignments = assignment_matrix, available = available)
}

# (#{|null| >= |observed|} + 1) / (#null + 1) for every observed statistic.
# Counting on the sorted null is exact and O((m + n) log n) instead of O(m n).
.empirical_p_values <- function(observed, null) {
  sorted_null <- sort(abs(null))
  at_least <- length(sorted_null) -
    findInterval(abs(observed), sorted_null, left.open = TRUE)
  (at_least + 1) / (length(sorted_null) + 1)
}

.attach_mscore_metadata <- function(regions, qc, cpg_files) {
  S4Vectors::metadata(regions)$qc <- qc
  S4Vectors::metadata(regions)$cpg_annotation_files <-
    unique(unname(cpg_files))
  S4Vectors::metadata(regions)$coordinate_system <- "1-based inclusive"
  regions
}

.align_design_data <- function(design_data, sample_names) {
  if ("sample" %in% names(design_data) &&
      setequal(as.character(design_data$sample), sample_names)) {
    design_data <- design_data[
      match(sample_names, as.character(design_data$sample)), , drop = FALSE
    ]
  } else if (!is.null(rownames(design_data)) &&
             setequal(rownames(design_data), sample_names)) {
    design_data <- design_data[match(sample_names, rownames(design_data)),
                               , drop = FALSE]
  }
  rownames(design_data) <- sample_names
  design_data
}

.find_group_column <- function(design_data, formula, coef) {
  available <- intersect(all.vars(formula), names(design_data))
  is_binary <- vapply(available, function(variable) {
    length(unique(design_data[[variable]][!is.na(design_data[[variable]])])) == 2L
  }, logical(1))
  binary_variables <- available[is_binary]
  preferred <- intersect(c("groups", "group"), binary_variables)
  if (length(preferred) == 1L) return(preferred)
  if (is.character(coef) && length(coef) == 1L) {
    matching <- binary_variables[vapply(binary_variables, function(variable) {
      grepl(variable, coef, fixed = TRUE)
    }, logical(1))]
    if (length(matching) == 1L) return(matching)
  }
  if (length(binary_variables) == 1L) return(binary_variables)
  stop(
    "Could not uniquely infer the two-level grouping variable from `formula`; ",
    "name it `group` or `groups` in `design_data`.", call. = FALSE
  )
}

.make_mscore_contrast <- function(coef, coefficient_names) {
  n_coef <- length(coefficient_names)
  if (is.character(coef) && length(coef) == 1L && !is.na(coef)) {
    coefficient_index <- match(coef, coefficient_names)
    if (is.na(coefficient_index)) {
      stop("`coef` was not found. Available coefficients: ",
           paste(coefficient_names, collapse = ", "), call. = FALSE)
    }
    contrast <- numeric(n_coef)
    contrast[coefficient_index] <- 1
    return(contrast)
  }
  if (is.numeric(coef) && length(coef) == 1L && is.finite(coef) &&
      coef == as.integer(coef) && coef >= 1L && coef <= n_coef) {
    contrast <- numeric(n_coef)
    contrast[as.integer(coef)] <- 1
    return(contrast)
  }
  if (is.numeric(coef) && length(coef) == n_coef &&
      all(is.finite(coef)) && any(coef != 0)) return(as.numeric(coef))
  stop("`coef` must identify one coefficient or supply a valid contrast vector.",
       call. = FALSE)
}

.build_region_matrices <- function(master_read_dt, cpg_cache, regions,
                                   sample_names) {
  empty_result <- function(reads_without_local_cpg = 0L,
                           dropped_regions = length(regions)) {
    list(
      regions = regions[0],
      S_matrix = matrix(numeric(), 0, length(sample_names)),
      T_matrix = matrix(numeric(), 0, length(sample_names)),
      N_sq_matrix = matrix(numeric(), 0, length(sample_names)),
      qc = list(
        reads_without_local_cpg = as.integer(reads_without_local_cpg),
        regions_dropped_incomplete_coverage = as.integer(dropped_regions)
      )
    )
  }
  if (!length(regions) || !nrow(master_read_dt)) return(empty_result())
  region_table <- data.table::data.table(
    chr = as.character(GenomicRanges::seqnames(regions)),
    region_start = as.numeric(IRanges::start(regions)),
    region_end = as.numeric(IRanges::end(regions)),
    region_id = seq_along(regions)
  )
  read_ranges <- master_read_dt[, .(
    chr, overlap_start = read_start, overlap_end = read_end,
    sample, annotation_key, count, is_methylated_global,
    cpg_index_start, cpg_index_end
  )]
  data.table::setkey(region_table, chr, region_start, region_end)
  hits <- data.table::foverlaps(
    read_ranges, region_table,
    by.x = c("chr", "overlap_start", "overlap_end"),
    by.y = c("chr", "region_start", "region_end"),
    type = "any", nomatch = NULL
  )
  if (!nrow(hits)) return(empty_result())
  hits[, ncpg_local := 0L]
  hit_groups <- unique(hits[, .(annotation_key, chr)])
  for (group_index in seq_len(nrow(hit_groups))) {
    key_value <- hit_groups$annotation_key[group_index]
    chromosome <- hit_groups$chr[group_index]
    row_index <- which(
      hits$annotation_key == key_value & hits$chr == chromosome
    )
    positions <- cpg_cache[[key_value]][[chromosome]]
    if (is.null(positions) || !length(positions)) next
    region_first <- findInterval(hits$region_start[row_index] - 1, positions) + 1L
    region_last <- findInterval(hits$region_end[row_index], positions)
    local_first <- pmax(hits$cpg_index_start[row_index], region_first)
    local_last <- pmin(hits$cpg_index_end[row_index], region_last)
    hits$ncpg_local[row_index] <- pmax(0L, local_last - local_first + 1L)
  }
  reads_without_local_cpg <- sum(hits$ncpg_local == 0L)
  hits <- hits[ncpg_local > 0L]
  if (!nrow(hits)) return(empty_result(reads_without_local_cpg))
  summaries <- hits[, .(
    S = sum(count * ncpg_local * as.integer(is_methylated_global)),
    T = sum(count * ncpg_local),
    N_sq = sum(count * ncpg_local^2)
  ), by = .(region_id, sample)]

  matrix_dim <- c(length(regions), length(sample_names))
  matrix_names <- list(paste0("region", seq_along(regions)), sample_names)
  S_matrix <- matrix(NA_real_, matrix_dim[1L], matrix_dim[2L],
                     dimnames = matrix_names)
  T_matrix <- S_matrix
  N_sq_matrix <- S_matrix
  matrix_index <- cbind(
    summaries$region_id, match(summaries$sample, sample_names)
  )
  S_matrix[matrix_index] <- summaries$S
  T_matrix[matrix_index] <- summaries$T
  N_sq_matrix[matrix_index] <- summaries$N_sq
  complete <- rowSums(
    is.finite(S_matrix) & is.finite(T_matrix) & is.finite(N_sq_matrix) &
      T_matrix > 0 & N_sq_matrix > 0
  ) == length(sample_names)
  dropped <- sum(!complete)
  if (dropped) {
    message("  Dropping ", dropped,
            " region(s) without CpG coverage in every sample.")
  }
  list(
    regions = regions[complete],
    S_matrix = S_matrix[complete, , drop = FALSE],
    T_matrix = T_matrix[complete, , drop = FALSE],
    N_sq_matrix = N_sq_matrix[complete, , drop = FALSE],
    qc = list(
      reads_without_local_cpg = as.integer(reads_without_local_cpg),
      regions_dropped_incomplete_coverage = as.integer(dropped)
    )
  )
}
