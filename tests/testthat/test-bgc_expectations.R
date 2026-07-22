# Tests of the functions in bgc_expectations.R

# --- Tests for ens_neutral ---
test_that("ens_neutral returns correct expected SFS for neutral model", {
  theta <- 100
  j <- 1:5
  expected <- theta / j
  expect_equal(ens_neutral(theta, j), expected, tolerance = 1e-6)
})

test_that("ens_neutral throws error for invalid inputs", {
  expect_error(ens_neutral(-10, 1), "theta must be a positive numeric value")
  expect_error(ens_neutral(10, 0), "j must be a positive integer")
  expect_error(ens_neutral(10, "abc"), "j must be a positive integer")
})

# --- Tests for ens_constant ---
test_that("ens_constant returns correct expected SFS for constant gBGC model", {
  theta <- 100
  B <- 1.2
  n <- 10
  j <- 2
  expected <- ((n * theta) / (j * (n - j))) * ((1 - exp(-B * (1 - j / n))) / (1 - exp(-B)))
  expect_equal(ens_constant(theta, B, n, j), expected, tolerance = 1e-6)
})

test_that("ens_constant returns ens_neutral when B=0", {
  theta <- 100
  n <- 10
  j <- 1:5
  expect_equal(ens_constant(theta, 0, n, j), ens_neutral(theta, j), tolerance = 1e-6)
})

test_that("ens_constant throws error for invalid inputs", {
  expect_error(ens_constant(-10, 1, 10, 1), "theta must be a positive numeric value")
  expect_error(ens_constant(10, "abc", 10, 1), "B must be a positive numeric value")
  expect_error(ens_constant(10, 1, 1, 1), "n must be an integer higher than 1")
  expect_error(ens_constant(10, 1, 10, 0), "j must be a positive integer lower than n")
  expect_error(ens_constant(10, 1, 10, 10), "j must be a positive integer lower than n")
})

# --- Tests for ens_hotspot1 ---
test_that("ens_hotspot1 returns correct expected SFS for first hotspot model", {
  theta <- 100
  B <- 1.2
  f <- 0.2
  n <- 10
  j <- 1:5
  expected <- f * ens_constant(theta, B, n, j) + (1 - f) * ens_neutral(theta, j)
  expect_equal(ens_hotspot1(theta, B, f, n, j), expected, tolerance = 1e-6)
})

test_that("ens_hotspot1 returns ens_neutral when B=0", {
  theta <- 100
  f <- 0.2
  n <- 10
  j <- 1:5
  expect_equal(ens_hotspot1(theta, 0, f, n, j), ens_neutral(theta, j), tolerance = 1e-6)
})

test_that("ens_hotspot1 throws error for invalid inputs", {
  expect_error(ens_hotspot1(-10, 1, 0.5, 10, 1), "theta must be a positive numeric value")
  expect_error(ens_hotspot1(10, "abc", 0.5, 10, 1), "B must be a positive numeric value")
  expect_error(ens_hotspot1(10, 1, -0.1, 10, 1), "f must be a numeric value between 0 and 1")
  expect_error(ens_hotspot1(10, 1, 1.1, 10, 1), "f must be a numeric value between 0 and 1")
  expect_error(ens_hotspot1(10, 1, 0.5, 1, 1), "n must be an integer higher than 1")
  expect_error(ens_hotspot1(10, 1, 0.5, 10, 0), "j must be a positive integer lower than n")
})

# --- Tests for ens_hotspot2 ---
test_that("ens_hotspot2 returns correct expected SFS for second hotspot model", {
  theta <- 100
  B0 <- 0.2
  B1 <- 3.5
  f <- 0.2
  n <- 10
  j <- 1:5
  expected <- f * ens_constant(theta, B1, n, j) + (1 - f) * ens_constant(theta, B0, n, j)
  expect_equal(ens_hotspot2(theta, B0, B1, f, n, j), expected, tolerance = 1e-6)
})

test_that("ens_hotspot2 returns ens_neutral when B0=0 and B1=0", {
  theta <- 100
  f <- 0.2
  n <- 10
  j <- 1:5
  expect_equal(ens_hotspot2(theta, 0, 0, f, n, j), ens_neutral(theta, j), tolerance = 1e-6)
})

test_that("ens_hotspot2 throws error for invalid inputs", {
  expect_error(ens_hotspot2(-10, 1, 1, 0.5, 10, 1), "theta must be a positive numeric value")
  expect_error(ens_hotspot2(10, "abc", 1, 0.5, 10, 1), "B0 must be a positive numeric value")
  expect_error(ens_hotspot2(10, 1, "abc", 0.5, 10, 1), "B1 must be a positive numeric value")
  expect_error(ens_hotspot2(10, 1, 1, -0.1, 10, 1), "f must be a numeric value between 0 and 1")
  expect_error(ens_hotspot2(10, 1, 1, 0.5, 1, 1), "n must be an integer higher than 1")
})

# --- Tests for ens_neutral_err ---
test_that("ens_neutral_err returns correct expected SFS for neutral model with polarization error", {
  theta1 <- 100
  theta2 <- 200
  e1 <- 0.01
  e2 <- 0.02
  n <- 10
  j <- 1:5
  expected <- (1 - e1) * ens_neutral(theta1, j) + e2 * ens_neutral(theta2, n - j)
  expect_equal(ens_neutral_err(theta1, theta2, e1, e2, n, j), expected, tolerance = 1e-6)
})

test_that("ens_neutral_err throws error for invalid inputs", {
  expect_error(ens_neutral_err(-10, 10, 0.1, 0.1, 10, 1), "theta1 must be a positive numeric value")
  expect_error(ens_neutral_err(10, -10, 0.1, 0.1, 10, 1), "theta2 must be a positive numeric value")
  expect_error(ens_neutral_err(10, 10, -0.1, 0.1, 10, 1), "e1 must be a numeric value between 0 and 1")
  expect_error(ens_neutral_err(10, 10, 1.1, 0.1, 10, 1), "e1 must be a numeric value between 0 and 1")
  expect_error(ens_neutral_err(10, 10, 0.1, -0.1, 10, 1), "e2 must be a numeric value between 0 and 1")
  expect_error(ens_neutral_err(10, 10, 0.1, 1.1, 10, 1), "e2 must be a numeric value between 0 and 1")
  expect_error(ens_neutral_err(10, 10, 0.1, 0.1, 1, 1), "n must be an integer higher than 1")
})

# --- Tests for ens_constant_err ---
test_that("ens_constant_err returns correct expected SFS for constant gBGC with polarization error", {
  theta1 <- 100
  theta2 <- 200
  B <- 1.2
  e1 <- 0.01
  e2 <- 0.02
  n <- 10
  j <- 1:5
  expected <- (1 - e1) * ens_constant(theta1, B, n, j) + e2 * ens_constant(theta2, -B, n, n - j)
  expect_equal(ens_constant_err(theta1, theta2, B, e1, e2, n, j), expected, tolerance = 1e-6)
})

test_that("ens_constant_err throws error for invalid inputs", {
  expect_error(ens_constant_err(-10, 10, 1, 0.1, 0.1, 10, 1), "theta1 must be a positive numeric value")
  expect_error(ens_constant_err(10, -10, 1, 0.1, 0.1, 10, 1), "theta2 must be a positive numeric value")
  expect_error(ens_constant_err(10, 10, "abc", 0.1, 0.1, 10, 1), "B0 must be a positive numeric value")
  expect_error(ens_constant_err(10, 10, 1, -0.1, 0.1, 10, 1), "e1 must be a numeric value between 0 and 1")
  expect_error(ens_constant_err(10, 10, 1, 0.1, -0.1, 10, 1), "e2 must be a numeric value between 0 and 1")
  expect_error(ens_constant_err(10, 10, 1, 0.1, 0.1, 1, 1), "n must be an integer higher than 1")
})

# --- Tests for ens_hotspot1_err ---
test_that("ens_hotspot1_err returns correct expected SFS for first hotspot model with polarization error", {
  theta1 <- 100
  theta2 <- 200
  B <- 1.2
  f <- 0.2
  e1 <- 0.01
  e2 <- 0.02
  n <- 10
  j <- 1:5
  expected <- (1 - e1) * ens_hotspot1(theta1, B, f, n, j) + e2 * ens_hotspot1(theta2, -B, f, n, n - j)
  expect_equal(ens_hotspot1_err(theta1, theta2, B, f, e1, e2, n, j), expected, tolerance = 1e-6)
})

test_that("ens_hotspot1_err throws error for invalid inputs", {
  expect_error(ens_hotspot1_err(-10, 10, 1, 0.5, 0.1, 0.1, 10, 1), "theta1 must be a positive numeric value")
  expect_error(ens_hotspot1_err(10, -10, 1, 0.5, 0.1, 0.1, 10, 1), "theta2 must be a positive numeric value")
  expect_error(ens_hotspot1_err(10, 10, "abc", 0.5, 0.1, 0.1, 10, 1), "B must be a positive numeric value")
  expect_error(ens_hotspot1_err(10, 10, 1, -0.1, 0.1, 0.1, 10, 1), "f must be a numeric value between 0 and 1")
  expect_error(ens_hotspot1_err(10, 10, 1, 0.5, -0.1, 0.1, 10, 1), "e1 must be a numeric value between 0 and 1")
  expect_error(ens_hotspot1_err(10, 10, 1, 0.5, 0.1, -0.1, 10, 1), "e2 must be a numeric value between 0 and 1")
  expect_error(ens_hotspot1_err(10, 10, 1, 0.5, 0.1, 0.1, 1, 1), "n must be an integer higher than 1")
})

# --- Tests for ens_hotspot2_err ---
test_that("ens_hotspot2_err returns correct expected SFS for second hotspot model with polarization error", {
  theta1 <- 100
  theta2 <- 200
  B0 <- 0.2
  B1 <- 3.5
  f <- 0.2
  e1 <- 0.01
  e2 <- 0.02
  n <- 10
  j <- 1:5
  expected <- (1 - e1) * ens_hotspot2(theta1, B0, B1, f, n, j) +
    e2 * ens_hotspot2(theta2, -B0, -B1, f, n, n - j)
  expect_equal(ens_hotspot2_err(theta1, theta2, B0, B1, f, e1, e2, n, j), expected, tolerance = 1e-6)
})

test_that("ens_hotspot2_err throws error for invalid inputs", {
  expect_error(ens_hotspot2_err(-10, 10, 1, 1, 0.5, 0.1, 0.1, 10, 1), "theta1 must be a positive numeric value")
  expect_error(ens_hotspot2_err(10, -10, 1, 1, 0.5, 0.1, 0.1, 10, 1), "theta2 must be a positive numeric value")
  expect_error(ens_hotspot2_err(10, 10, "abc", 1, 0.5, 0.1, 0.1, 10, 1), "B0 must be a positive numeric value")
  expect_error(ens_hotspot2_err(10, 10, 1, "abc", 0.5, 0.1, 0.1, 10, 1), "B1 must be a positive numeric value")
  expect_error(ens_hotspot2_err(10, 10, 1, 1, -0.1, 0.1, 0.1, 10, 1), "f must be a numeric value between 0 and 1")
  expect_error(ens_hotspot2_err(10, 10, 1, 1, 0.5, -0.1, 0.1, 10, 1), "e1 must be a numeric value between 0 and 1")
  expect_error(ens_hotspot2_err(10, 10, 1, 1, 0.5, 0.1, -0.1, 10, 1), "e2 must be a numeric value between 0 and 1")
  expect_error(ens_hotspot2_err(10, 10, 1, 1, 0.5, 0.1, 0.1, 1, 1), "n must be an integer higher than 1")
})
