test_that("two-stage GLS returns stable estimates on mock data", {
  S_matrix <- rbind(
    region_a = c(2, 3, 7, 8),
    region_b = c(1, 2, 6, 7)
  )
  T_matrix <- matrix(
    c(10, 12, 10, 12, 8, 10, 8, 10),
    nrow = 2, byrow = TRUE,
    dimnames = dimnames(S_matrix)
  )
  N_sq_matrix <- T_matrix^2 / matrix(
    c(5, 6, 5, 6, 4, 5, 4, 5),
    nrow = 2, byrow = TRUE
  )
  design <- cbind(intercept = 1, group = c(0, 0, 1, 1))

  fit <- fit_mscore_gls(
    S_matrix, T_matrix, N_sq_matrix,
    design = design, contrast = c(0, 1)
  )

  expect_s3_class(fit, "mscore_gls_fit")
  expect_equal(fit$kappa, T_matrix^2 / N_sq_matrix)
  expect_equal(dim(fit$coefficients), c(2, 2))
  expect_equal(dim(fit$covariance), c(2, 2, 2))
  expect_true(all(fit$phi >= 0.001 & fit$phi <= 0.999))
  expect_true(all(fit$coefficients[, "group"] > 0))
  expect_true(all(is.finite(fit$statistic)))

  expected_z <- asin(2 * (S_matrix + 0.5) / (T_matrix + 1) - 1)
  expect_equal(fit$transformed, expected_z)

  stage1_covariance <- solve(crossprod(design, design * fit$kappa[1, ]))
  expected_stage1_beta <- drop(stage1_covariance %*%
    crossprod(design, expected_z[1, ] * fit$kappa[1, ]))
  stage1_residual <- expected_z[1, ] - drop(design %*% expected_stage1_beta)
  expected_sigma_sq <- sum(fit$kappa[1, ] * stage1_residual^2) /
    fit$residual_df
  expected_phi <- 4 * (expected_sigma_sq - 1) /
    sum(fit$kappa[1, ] - 1)
  expected_phi <- min(0.999, max(0.001, expected_phi))

  expect_equal(unname(fit$stage1_coefficients[1, ]), unname(expected_stage1_beta))
  expect_equal(unname(fit$phi[1]), unname(expected_phi))

  stage2_weights <- fit$kappa[1, ] /
    (1 + (fit$kappa[1, ] - 1) * fit$phi[1])
  expected_covariance <- solve(crossprod(design, design * stage2_weights))
  expected_beta <- drop(expected_covariance %*%
    crossprod(design, expected_z[1, ] * stage2_weights))
  expected_statistic <- expected_beta[2] / sqrt(expected_covariance[2, 2])

  expect_equal(unname(fit$coefficients[1, ]), unname(expected_beta))
  expect_equal(unname(fit$covariance[, , 1]), unname(expected_covariance))
  expect_equal(unname(fit$statistic[1]), unname(expected_statistic))
})

test_that("fit_mscore_gls validates incompatible inputs", {
  design <- cbind(1, c(0, 0, 1, 1))

  expect_error(
    fit_mscore_gls(matrix(1, 1, 4), matrix(2, 1, 3),
                   matrix(1, 1, 4), design, c(0, 1)),
    "identical dimensions"
  )
  expect_error(
    fit_mscore_gls(matrix(3, 1, 4), matrix(2, 1, 4),
                   matrix(2, 1, 4), design, c(0, 1)),
    "between zero"
  )
})
