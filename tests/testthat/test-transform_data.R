# Test for the functions in transform_data.R

# Test for t_sfs ----------------------------------------------------------

test_that("t_sfs throws error for invalid inputs", {
  # Non-numeric input
  expect_error(t_sfs("abc"), "sfs must be a numeric vector")
  expect_error(t_sfs(list(1, 2, 3)), "sfs must be a numeric vector")
  expect_error(t_sfs(NULL), "sfs must be a numeric vector")
})

test_that("t_sfs returns correct transformed SFS", {
  # Test with example from documentation
  sfs <- c(100, 50, 30, 20, 3)
  expected_t_sfs <- log(sfs) -0.05/sfs + 4.6/(sfs^2) -5.4/(sfs^3)
  expect_equal(t_sfs(sfs), expected_t_sfs, tolerance = 1e-6)

  # Test with single-value SFS
  sfs_single <- c(50)
  expected_single <- log(50) -0.05/50 + 4.6/(50^2) -5.4/(50^3)
  expect_equal(t_sfs(sfs_single), expected_single, tolerance = 1e-6)

  # Test with large values
  sfs_large <- c(1000, 2000, 3000)
  expected_large <- log(sfs_large) -0.05/sfs_large + 4.6/(sfs_large^2) -5.4/(sfs_large^3)
  expect_equal(t_sfs(sfs_large), expected_large, tolerance = 1e-6)

  # Test with small values (but > 0)
  sfs_small <- c(0.1, 0.5, 1)
  expected_small <- log(sfs_small) -0.05/sfs_small + 4.6/(sfs_small^2) -5.4/(sfs_small^3)
  expect_equal(t_sfs(sfs_small), expected_small, tolerance = 1e-6)
})

test_that("t_sfs handles edge cases", {
  # Test with very small positive values (close to 0)
  sfs_near_zero <- c(0.001, 0.01, 0.1)
  expect_silent(t_sfs(sfs_near_zero))  # Should not throw an error
  expect_true(all(is.finite(t_sfs(sfs_near_zero))))  # Should return finite values

  # Test with very large values
  sfs_very_large <- c(1e6, 1e7, 1e8)
  expect_silent(t_sfs(sfs_very_large))
  expect_true(all(is.finite(t_sfs(sfs_very_large))))
})

# Test for t_variance ------------------------------------------------------

test_that("t_variance throws error for invalid inputs", {
  # Non-numeric input
  expect_error(t_variance("abc"), "sfs must be a numeric vector")
  expect_error(t_variance(list(1, 2, 3)), "sfs must be a numeric vector")
  expect_error(t_variance(NULL), "sfs must be a numeric vector")
})

test_that("t_variance returns correct expected variance", {
  # Test with example from documentation
  sfs <- c(100, 50, 30, 20, 3)
  expected_t_var <- 1/sfs + 0.1/(sfs)^2 + c(0,0,0,0,0.43) # adding the correcting factor
  expect_equal(t_variance(sfs), expected_t_var, tolerance = 1e-6)

  # Test with single-value SFS
  sfs_single <- c(50)
  expected_single_var <- 1/50 + 0.1/(50)^2
  expect_equal(t_variance(sfs_single), expected_single_var, tolerance = 1e-6)

  # Test with large values
  sfs_large <- c(1000, 2000, 3000)
  expected_large_var <- 1/sfs_large + 0.1/(sfs_large^2)
  expect_equal(t_variance(sfs_large), expected_large_var, tolerance = 1e-6)
})

test_that("t_variance handles edge cases", {
  # Test with very small positive values (close to 0)
  sfs_near_zero <- c(0.001, 0.01, 0.1)
  expect_silent(t_variance(sfs_near_zero))
  expect_true(all(is.finite(t_variance(sfs_near_zero))))

  # Test with very large values
  sfs_very_large <- c(1e6, 1e7, 1e8)
  expect_silent(t_variance(sfs_very_large))
  expect_true(all(is.finite(t_variance(sfs_very_large))))
})

# Integration tests --------------------------------------------------------

test_that("t_sfs and t_variance work together", {
  # Test the example from the documentation
  sfs_ws <- c(100, 50, 30, 20, 3)
  sfs_sw <- c(300, 140, 80, 50, 7)

  # Compute transformed SFS
  t_ratio <- t_sfs(sfs_ws) - t_sfs(sfs_sw)

  # Compute weights (inverse of variance)
  weight <- 1 / (t_variance(sfs_ws) + t_variance(sfs_sw))

  # Check that weight is positive and finite
  expect_true(all(weight > 0))
  expect_true(all(is.finite(weight)))

  # Check that t_ratio is finite
  expect_true(all(is.finite(t_ratio)))
})
