# Test of the funcvtions in simulation.R

# --- Tests for expected_sfs_constant ---
test_that("expected_sfs_constant returns correct structure for valid inputs", {
  n <- 20
  theta <- 100
  B <- 1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  result <- expected_sfs_constant(n, theta, B, mut_bias, GC, eWS, eSW, vect_r)
  expect_type(result, "list")
  expect_equal(names(result), c("WS", "SW"))
  expect_type(result$WS, "double")
  expect_type(result$SW, "double")
  expect_length(result$WS, n-1)
  expect_length(result$SW, n-1)
})

test_that("expected_sfs_constant throws error for invalid inputs", {
  n <- 20
  theta <- 100
  B <- 1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  expect_error(expected_sfs_constant("abc", theta, B, mut_bias, GC, eWS, eSW, vect_r), "n must be a positive integer higher than 1")
  expect_error(expected_sfs_constant(1, theta, B, mut_bias, GC, eWS, eSW,vect_r), "n must be a positive integer higher than 1")
  expect_error(expected_sfs_constant(n, -1, B, mut_bias, GC, eWS, eSW, vect_r), "theta must be a positive numerical value")
  expect_error(expected_sfs_constant(n, theta, "abc", mut_bias, GC, eWS, eSW, vect_r), "B must be a numerical value")
  expect_error(expected_sfs_constant(n, theta, B, -1, GC, eWS, eSW, vect_r), "mut_bias must be a positive numerical value")
  expect_error(expected_sfs_constant(n, theta, B, mut_bias, 0, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(expected_sfs_constant(n, theta, B, mut_bias, 1, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(expected_sfs_constant(n, theta, B, mut_bias, GC, -0.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_constant(n, theta, B, mut_bias, GC, 1.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_constant(n, theta, B, mut_bias, GC, eWS, -0.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_constant(n, theta, B, mut_bias, GC, eWS, 1.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_constant(n, theta, B, mut_bias, GC, eWS, eSW, c(1,1,1)), "vect_r must be a vector of positive numerical values of length n-1",fixed=T)
})

# --- Tests for expected_sfs_hotspot ---
test_that("expected_sfs_hotspot returns correct structure for valid inputs", {
  n <- 20
  theta <- 100
  B0 <- 3
  B1 <- 0.1
  f <- 0.1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  result <- expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, eSW, vect_r)
  expect_type(result, "list")
  expect_equal(names(result), c("WS", "SW"))
  expect_type(result$WS, "double")
  expect_type(result$SW, "double")
  expect_length(result$WS, n-1)
  expect_length(result$SW, n-1)
})

test_that("expected_sfs_hotspot throws error for invalid inputs", {
  n <- 20
  theta <- 100
  B0 <- 3
  B1 <- 0.1
  f <- 0.1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  expect_error(expected_sfs_hotspot("abc", theta, B0, B1, f, mut_bias, GC, eWS, eSW, vect_r), "n must be a positive integer higher than 1")
  expect_error(expected_sfs_hotspot(n, -1, B0, B1, f, mut_bias, GC, eWS, eSW, vect_r), "theta must be a positive numerical value")
  expect_error(expected_sfs_hotspot(n, theta, "abc", B1, f, mut_bias, GC, eWS, eSW, vect_r), "B0 must be a numerical value")
  expect_error(expected_sfs_hotspot(n, theta, B0, "abc", f, mut_bias, GC, eWS, eSW, vect_r), "B1 must be a numerical value")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, -0.1, mut_bias, GC, eWS, eSW, vect_r), "f must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, 1.1, mut_bias, GC, eWS, eSW, vect_r), "f must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, -1, GC, eWS, eSW, vect_r), "mut_bias must be a positive numerical value")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, 0, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, 1, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, -0.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, 1.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, -0.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, 1.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(expected_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, eSW, c(1, 1)), "vect_r must be a vector of positive numerical values of length n-1")
})

# --- Tests for simulated_sfs_constant ---
test_that("simulated_sfs_constant returns correct structure for valid inputs", {
  set.seed(123)  # For reproducibility
  n <- 20
  theta <- 100
  B <- 1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  # Test "free" option
  result <- simulated_sfs_constant(n, theta, B, mut_bias, GC, eWS, eSW, vect_r, fix_snp = "free")
  expect_type(result, "list")
  expect_equal(names(result), c("WS", "SW"))
  expect_type(result$WS, "integer")
  expect_type(result$SW, "integer")

  # Test "fixed_sfs" option
  result <- simulated_sfs_constant(n, theta, B, mut_bias, GC, eWS, eSW, vect_r, fix_snp = "fixed_sfs")
  expect_type(result, "list")
  expect_equal(names(result), c("WS", "SW"))

  # Test "fixed_tot" option
  result <- simulated_sfs_constant(n, theta, B, mut_bias, GC, eWS, eSW, vect_r, fix_snp = "fixed_tot")
  expect_type(result, "list")
  expect_equal(names(result), c("WS", "SW"))
})

test_that("simulated_sfs_constant throws error for invalid inputs", {
  n <- 20
  theta <- 100
  B <- 1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  expect_error(simulated_sfs_constant("abc", theta, B, mut_bias, GC, eWS, eSW, vect_r), "n must be a positive integer higher than 1")
  expect_error(simulated_sfs_constant(n, -1, B, mut_bias, GC, eWS, eSW, vect_r), "theta must be a positive numerical value")
  expect_error(simulated_sfs_constant(n, theta, "abc", mut_bias, GC, eWS, eSW, vect_r), "B must be a numerical value")
  expect_error(simulated_sfs_constant(n, theta, B, -1, GC, eWS, eSW, vect_r), "mut_bias must be a positive numerical value")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, 0, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, 1, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, GC, -0.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, GC, 1.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, GC, eWS, -0.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, GC, eWS, 1.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, GC, eWS, eSW, c(1, 1)), "vect_r must be a vector of positive numerical values of length n-1")
  expect_error(simulated_sfs_constant(n, theta, B, mut_bias, GC, eWS, eSW, vect_r, fix_snp = "invalid"), "The chosen option is not available: must be free, fixed_sfs or fixed_tot")
})

# --- Tests for simulated_sfs_hotspot ---
test_that("simulated_sfs_hotspot returns correct structure for valid inputs", {
  set.seed(123)  # For reproducibility
  n <- 20
  theta <- 100
  B0 <- 3
  B1 <- 0.1
  f <- 0.1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  # Test "free" option
  result <- simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, eSW, vect_r, fix_snp = "free")
  expect_type(result, "list")
  expect_equal(names(result), c("WS", "SW"))
  expect_type(result$WS, "integer")
  expect_type(result$SW, "integer")

  # Test "fixed_sfs" option
  result <- simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, eSW, vect_r, fix_snp = "fixed_sfs")
  expect_type(result, "list")
  expect_equal(names(result), c("WS", "SW"))
})

test_that("simulated_sfs_hotspot throws error for invalid inputs", {
  n <- 20
  theta <- 100
  B0 <- 3
  B1 <- 0.1
  f <- 0.1
  mut_bias <- 2
  GC <- 0.5
  eWS <- 0.01
  eSW <- 0.02
  vect_r <- rep(1, n-1)

  expect_error(simulated_sfs_hotspot("abc", theta, B0, B1, f, mut_bias, GC, eWS, eSW, vect_r), "n must be a positive integer higher than 1")
  expect_error(simulated_sfs_hotspot(n, -1, B0, B1, f, mut_bias, GC, eWS, eSW, vect_r), "theta must be a positive numerical value")
  expect_error(simulated_sfs_hotspot(n, theta, "abc", B1, f, mut_bias, GC, eWS, eSW, vect_r), "B0 must be a numerical value")
  expect_error(simulated_sfs_hotspot(n, theta, B0, "abc", f, mut_bias, GC, eWS, eSW, vect_r), "B1 must be a numerical value")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, -0.1, mut_bias, GC, eWS, eSW, vect_r), "f must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, 1.1, mut_bias, GC, eWS, eSW, vect_r), "f must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, -1, GC, eWS, eSW, vect_r), "mut_bias must be a positive numerical value")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, 0, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, 1, eWS, eSW, vect_r), "GC must be a numerical value strictly higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, -0.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, 1.1, eSW, vect_r), "eWS must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, -0.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, 1.1, vect_r), "eSW must be a numerical value higher than 0 and lower than 1")
  expect_error(simulated_sfs_hotspot(n, theta, B0, B1, f, mut_bias, GC, eWS, eSW, c(1, 1)), "vect_r must be a vector of positive numerical values of length n-1")
})

