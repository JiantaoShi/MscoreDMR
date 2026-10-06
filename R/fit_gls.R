#' Fit the two-stage M-score GLS model
#'
#' Fits the Beta-Bernoulli generalized least-squares model described for
#' candidate regions. Rows of the summary matrices are regions and columns are
#' samples. The effective coverage is calculated as
#' `T_matrix^2 / N_sq_matrix`, where `N_sq_matrix` contains
#' `sum_t(N_trd^2)` for each region and sample.
#'
#' @param S_matrix Numeric matrix of weighted methylated-read counts.
#' @param T_matrix Numeric matrix of total read weights. Must have the same
#'   dimensions as `S_matrix`.
#' @param N_sq_matrix Numeric matrix containing the sum of squared per-read CpG
#'   coverages. Must have the same dimensions as `S_matrix`.
#' @param design Numeric design matrix with samples in rows. It must be full
#'   column rank and have fewer columns than rows.
#' @param contrast Numeric vector defining the Wald contrast, with one value per
#'   column of `design`.
#' @param pseudo_count Non-negative pseudo-count used in the arcsine transform.
#'   The default is `0.5`.
#' @param phi_bounds Length-two numeric vector giving the lower and upper bounds
#'   for the dispersion estimate. Defaults to `c(0.001, 0.999)`.
#'
#' @return A list containing the transformed observations and effective
#'   coverages, Stage 1 and Stage 2 coefficient matrices, region-specific
#'   dispersions, Stage 2 covariance matrices, and Wald statistics and
#'   two-sided normal-reference p-values.
#' @export
#'
#' @examples
#' S <- matrix(c(2, 3, 7, 8, 1, 2, 6, 7), nrow = 2, byrow = TRUE)
#' T <- matrix(10, nrow = 2, ncol = 4)
#' N_sq <- matrix(20, nrow = 2, ncol = 4)
#' X <- cbind(intercept = 1, group = c(0, 0, 1, 1))
#' fit_mscore_gls(S, T, N_sq, X, c(0, 1))
fit_mscore_gls <- function(S_matrix, T_matrix, N_sq_matrix, design, contrast,
                           pseudo_count = 0.5,
                           phi_bounds = c(0.001, 0.999)) {
  summary_matrices <- list(
    S_matrix = S_matrix,
    T_matrix = T_matrix,
    N_sq_matrix = N_sq_matrix
  )

  if (!all(vapply(summary_matrices, is.matrix, logical(1)))) {
    stop("`S_matrix`, `T_matrix`, and `N_sq_matrix` must be matrices.",
         call. = FALSE)
  }
  if (!all(vapply(summary_matrices, is.numeric, logical(1)))) {
    stop("Summary matrices must be numeric.", call. = FALSE)
  }
  matrix_dims <- lapply(summary_matrices, dim)
  if (!all(vapply(matrix_dims[-1], identical, logical(1), matrix_dims[[1]]))) {
    stop("Summary matrices must have identical dimensions.", call. = FALSE)
  }
  if (nrow(S_matrix) < 1L || ncol(S_matrix) < 1L) {
    stop("Summary matrices must contain at least one region and sample.",
         call. = FALSE)
  }
  if (any(!is.finite(S_matrix)) || any(!is.finite(T_matrix)) ||
      any(!is.finite(N_sq_matrix))) {
    stop("Summary matrices cannot contain missing or infinite values.",
         call. = FALSE)
  }
  if (any(T_matrix <= 0) || any(N_sq_matrix <= 0)) {
    stop("`T_matrix` and `N_sq_matrix` must be strictly positive.",
         call. = FALSE)
  }
  if (any(S_matrix < 0) || any(S_matrix > T_matrix)) {
    stop("`S_matrix` values must lie between zero and `T_matrix`.",
         call. = FALSE)
  }

  if (!is.matrix(design) || !is.numeric(design) || any(!is.finite(design))) {
    stop("`design` must be a finite numeric matrix.", call. = FALSE)
  }
  if (nrow(design) != ncol(S_matrix)) {
    stop("The rows of `design` must match the samples in the summary matrices.",
         call. = FALSE)
  }
  n_samples <- nrow(design)
  n_coef <- ncol(design)
  if (n_coef < 1L || n_samples <= n_coef) {
    stop("`design` must have at least one column and positive residual degrees of freedom.",
         call. = FALSE)
  }
  if (qr(design)$rank < n_coef) {
    stop("`design` must have full column rank.", call. = FALSE)
  }

  if (!is.numeric(contrast) || length(contrast) != n_coef ||
      any(!is.finite(contrast))) {
    stop("`contrast` must be a finite numeric vector with one value per design column.",
         call. = FALSE)
  }
  contrast <- as.numeric(contrast)
  if (all(contrast == 0)) {
    stop("`contrast` must contain at least one non-zero value.", call. = FALSE)
  }
  if (!is.numeric(pseudo_count) || length(pseudo_count) != 1L ||
      !is.finite(pseudo_count) || pseudo_count < 0) {
    stop("`pseudo_count` must be one finite, non-negative number.",
         call. = FALSE)
  }
  if (!is.numeric(phi_bounds) || length(phi_bounds) != 2L ||
      any(!is.finite(phi_bounds)) || phi_bounds[1L] < 0 ||
      phi_bounds[2L] > 1 || phi_bounds[1L] > phi_bounds[2L]) {
    stop("`phi_bounds` must be two ordered values within [0, 1].",
         call. = FALSE)
  }

  kappa <- T_matrix^2 / N_sq_matrix
  if (any(!is.finite(kappa)) || any(kappa < 1 - sqrt(.Machine$double.eps))) {
    stop("Effective coverages computed from `T_matrix` and `N_sq_matrix` must be at least one.",
         call. = FALSE)
  }
  kappa <- pmax(kappa, 1)

  probability <- (S_matrix + pseudo_count) /
    (T_matrix + 2 * pseudo_count)
  transform_argument <- 2 * probability - 1
  transform_argument <- pmin(pmax(transform_argument, -1), 1)
  Z <- asin(transform_argument)

  region_names <- rownames(S_matrix)
  if (is.null(region_names)) {
    region_names <- paste0("region", seq_len(nrow(S_matrix)))
  }
  coefficient_names <- colnames(design)
  if (is.null(coefficient_names)) {
    coefficient_names <- paste0("beta", seq_len(n_coef))
  }

  fit_weighted <- function(response, weights) {
    information <- crossprod(design, design * weights)
    covariance <- tryCatch(
      solve(information),
      error = function(e) {
        stop("A region-specific GLS system is singular.", call. = FALSE)
      }
    )
    beta <- drop(covariance %*% crossprod(design, response * weights))
    list(beta = beta, covariance = covariance)
  }

  stage1_beta <- matrix(
    NA_real_, nrow = nrow(S_matrix), ncol = n_coef,
    dimnames = list(region_names, coefficient_names)
  )
  stage2_beta <- stage1_beta
  covariance <- array(
    NA_real_, dim = c(n_coef, n_coef, nrow(S_matrix)),
    dimnames = list(coefficient_names, coefficient_names, region_names)
  )
  phi <- statistic <- numeric(nrow(S_matrix))
  residual_df <- n_samples - n_coef

  for (region in seq_len(nrow(S_matrix))) {
    response <- Z[region, ]
    region_kappa <- kappa[region, ]

    stage1 <- fit_weighted(response, region_kappa)
    stage1_beta[region, ] <- stage1$beta
    residual <- response - drop(design %*% stage1$beta)
    sigma_sq <- sum(region_kappa * residual^2) / residual_df
    phi_denominator <- sum(region_kappa - 1)
    phi_raw <- if (phi_denominator <= sqrt(.Machine$double.eps)) {
      phi_bounds[1L]
    } else {
      n_samples * (sigma_sq - 1) / phi_denominator
    }
    phi[region] <- min(phi_bounds[2L], max(phi_bounds[1L], phi_raw))

    stage2_weights <- region_kappa /
      (1 + (region_kappa - 1) * phi[region])
    stage2 <- fit_weighted(response, stage2_weights)
    stage2_beta[region, ] <- stage2$beta
    covariance[, , region] <- stage2$covariance

    contrast_variance <- drop(
      crossprod(contrast, stage2$covariance %*% contrast)
    )
    statistic[region] <- drop(crossprod(contrast, stage2$beta)) /
      sqrt(contrast_variance)
  }

  names(phi) <- names(statistic) <- region_names
  p_value <- 2 * stats::pnorm(-abs(statistic))
  names(p_value) <- region_names
  dimnames(Z) <- dimnames(kappa) <- list(region_names, colnames(S_matrix))

  structure(
    list(
      coefficients = stage2_beta,
      covariance = covariance,
      statistic = statistic,
      p_value = p_value,
      phi = phi,
      stage1_coefficients = stage1_beta,
      transformed = Z,
      kappa = kappa,
      residual_df = residual_df,
      contrast = contrast
    ),
    class = "mscore_gls_fit"
  )
}
