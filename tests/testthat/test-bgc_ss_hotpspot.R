# Test of the functions in bgc_ss_hotspot.R

# Mock helper functions for testing
ratio_B <- function(B, x) {
  # Mock implementation: returns a simple linear function of x
  return(1 + B * x)
}

d_hotspot1 <- function(B, f, x) {
  # Mock implementation: returns zero gradients for testing
  n <- length(x)
  return(list(dB = rep(0, n), df = rep(0, n)))
}

d_hotspot2 <- function(B0, B1, f, x) {
  # Mock implementation: returns zero gradients for testing
  n <- length(x)
  return(list(dB0 = rep(0, n), dB1 = rep(0, n), df = rep(0, n)))
}

d_expected_log_ratio <- function(WS, SW, e1, e2) {
  # Mock implementation: returns zero gradients for testing
  n <- length(WS)
  return(list(d1 = rep(0, n), d2 = rep(0, n)))
}

# --- Tests for sum_of_squares_hotspot1 ---
test_that("sum_of_squares_hotspot1 returns non-negative value for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(1, 0.1, 2, 0.02, 0.01)
  GC <- 0.5
  result <- sum_of_squares_hotspot1(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("sum_of_squares_hotspot1 throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02), sfsWS, sfsSW, 0.5), "A five values vector must be given as par")
  expect_error(sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02, 0.01), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
  expect_error(sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW, 1), "GC content must be strictly between 0 and 1")
})

# --- Tests for gr_sum_of_squares_hotspot1 ---
test_that("gr_sum_of_squares_hotspot1 returns gradient vector of length 5 for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(1, 0.1, 2, 0.02, 0.01)
  GC <- 0.5
  result <- gr_sum_of_squares_hotspot1(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_length(result, 5)
})

test_that("gr_sum_of_squares_hotspot1 throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(gr_sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02), sfsWS, sfsSW, 0.5), "A five values vector must be given as par")
  expect_error(gr_sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02, 0.01), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(gr_sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(gr_sum_of_squares_hotspot1(c(1, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for sum_of_squares_hotspot2 ---
test_that("sum_of_squares_hotspot2 returns non-negative value for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(0.1, 5, 0.1, 2, 0.02, 0.01)
  GC <- 0.5
  result <- sum_of_squares_hotspot2(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("sum_of_squares_hotspot2 throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02), sfsWS, sfsSW, 0.5), "A six values vector must be given as par")
  expect_error(sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02, 0.01), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for gr_sum_of_squares_hotspot2 ---
test_that("gr_sum_of_squares_hotspot2 returns gradient vector of length 6 for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(0.1, 5, 0.1, 2, 0.02, 0.01)
  GC <- 0.5
  result <- gr_sum_of_squares_hotspot2(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_length(result, 6)
})

test_that("gr_sum_of_squares_hotspot2 throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(gr_sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02), sfsWS, sfsSW, 0.5), "A six values vector must be given as par")
  expect_error(gr_sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02, 0.01), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(gr_sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(gr_sum_of_squares_hotspot2(c(0.1, 5, 0.1, 2, 0.02, 0.01), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for sum_of_squares_hotspot2bis ---
test_that("sum_of_squares_hotspot2bis returns non-negative value for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(0.1, 5, 2, 0.02, 0.01)
  GC <- 0.5
  f <- 0.1
  result <- sum_of_squares_hotspot2bis(param, sfsWS, sfsSW, GC, f)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("sum_of_squares_hotspot2bis throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02), sfsWS, sfsSW, 0.5, 0.1), "A five values vector must be given as par")
  expect_error(sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), "abc", sfsSW, 0.5, 0.1), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), sfsWS, sfsSW[-1], 0.5, 0.1), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), sfsWS, sfsSW, 0, 0.1), "GC content must be strictly between 0 and 1")
  expect_error(sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), sfsWS, sfsSW, 0.5, 0), "f must be strictly between 0 and 1/2")
  expect_error(sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), sfsWS, sfsSW, 0.5, 0.6), "f must be strictly between 0 and 1/2")
})

# --- Tests for gr_sum_of_squares_hotspot2bis ---
test_that("gr_sum_of_squares_hotspot2bis returns gradient vector of length 5 for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(0.1, 5, 2, 0.02, 0.01)
  GC <- 0.5
  f <- 0.1
  result <- gr_sum_of_squares_hotspot2bis(param, sfsWS, sfsSW, GC, f)
  expect_type(result, "double")
  expect_length(result, 5)
})

test_that("gr_sum_of_squares_hotspot2bis throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(gr_sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02), sfsWS, sfsSW, 0.5, 0.1), "A five values vector must be given as par")
  expect_error(gr_sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), "abc", sfsSW, 0.5, 0.1), "The two SFSs must have positive numeric values")
  expect_error(gr_sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), sfsWS, sfsSW[-1], 0.5, 0.1), "The two SFSs, WS and SW, must have the same length")
  expect_error(gr_sum_of_squares_hotspot2bis(c(0.1, 5, 2, 0.02, 0.01), sfsWS, sfsSW, 0, 0.1), "GC content must be strictly between 0 and 1")
})
