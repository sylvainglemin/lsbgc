# Tests for the functions in bgc_optimization.R

# File: test_least_square.R

library(testthat)
library(lsbgc)  # Replace with your package name if different

# Mock constants for testing
MMIN <- -5
MMAX <- 5
BMIN <- -100
BMAX <- 100
MAXIT <- 100
FACTR <- 1e7
LMM <- 5
VERBOSE <- 0
USEGR <- FALSE
ONE <- 0.999
ZERO <- 0.001

# --- Tests for AICls ---
test_that("AICls returns correct AIC for valid inputs", {
  n <- 100
  np <- 3
  SSres <- 123
  expected_AIC <- n * log(SSres / n) + 2 * (np + 1) * (n + 2) / (n - np)
  expect_equal(AICls(n, np, SSres), expected_AIC, tolerance = 1e-6)
})

test_that("AICls handles edge cases for n, np, and SSres", {
  expect_equal(AICls(1, 1, 1), 1 * log(1 / 1) + 2 * (1 + 1) * (1 + 2) / (1 - 1), tolerance = 1e-6)
  expect_equal(AICls(10, 5, 10), 10 * log(10 / 10) + 2 * (5 + 1) * (10 + 2) / (10 - 5), tolerance = 1e-6)
})

test_that("AICls prints warnings for invalid inputs", {
  expect_error(AICls(0, 3, 123), "n must be strictly positive")
  expect_error(AICls(100, -1, 123), "np must be positive and lower than n")
  expect_error(AICls(100, 3, -1), "SSres must be positive")
})

# --- Tests for least_square_M ---
test_that("least_square_M returns correct structure for valid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  GC <- 0.5
  result <- least_square_M(sfsWS, sfsSW, GC)
  expect_type(result, "list")
  expect_equal(names(result), c("param", "criteria"))
  expect_equal(names(result$param), c("mutbias"))
  expect_equal(names(result$criteria), c("SStot", "SSres", "R2", "AIC"))
})

test_that("least_square_M throws error for invalid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  expect_error(least_square_M("abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(least_square_M(sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(least_square_M(sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
  expect_error(least_square_M(sfsWS, sfsSW, 1), "GC content must be strictly between 0 and 1")
})

# --- Tests for least_square_B ---
test_that("least_square_B returns correct structure for valid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  GC <- 0.5
  result <- least_square_B(sfsWS, sfsSW, GC)
  expect_type(result, "list")
  expect_equal(names(result), c("param", "criteria"))
  expect_equal(names(result$param), c("B"))
  expect_equal(names(result$criteria), c("SStot", "SSres", "R2", "AIC"))
})

test_that("least_square_B throws error for invalid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  expect_error(least_square_B("abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(least_square_B(sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(least_square_B(sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for least_square_BM ---
test_that("least_square_BM returns correct structure for valid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  GC <- 0.5
  result <- least_square_BM(sfsWS, sfsSW, GC)
  expect_type(result, "list")
  expect_equal(names(result), c("param", "criteria"))
  expect_equal(names(result$param), c("B", "mutbias"))
  expect_equal(names(result$criteria), c("SStot", "SSres", "R2", "AIC"))
})

test_that("least_square_BM throws error for invalid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  expect_error(least_square_BM("abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(least_square_BM(sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(least_square_BM(sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for least_square_hotspot1 ---
test_that("least_square_hotspot1 returns correct structure for valid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  GC <- 0.5
  result <- least_square_hotspot1(sfsWS, sfsSW, GC)
  expect_type(result, "list")
  expect_equal(names(result), c("param", "criteria"))
  expect_equal(names(result$param), c("B", "f", "mutbias"))
  expect_equal(names(result$criteria), c("SStot", "SSres", "R2", "AIC"))
})

test_that("least_square_hotspot1 throws error for invalid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  expect_error(least_square_hotspot1("abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(least_square_hotspot1(sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(least_square_hotspot1(sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for least_square_hotspot2 ---
test_that("least_square_hotspot2 returns correct structure for valid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  GC <- 0.5
  result <- least_square_hotspot2(sfsWS, sfsSW, GC)
  expect_type(result, "list")
  expect_equal(names(result), c("param", "criteria"))
  expect_equal(names(result$param), c("B0", "B1", "f", "mutbias"))
  expect_equal(names(result$criteria), c("SStot", "SSres", "R2", "AIC"))
})

test_that("least_square_hotspot2 throws error for invalid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  expect_error(least_square_hotspot2("abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(least_square_hotspot2(sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(least_square_hotspot2(sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for least_square_hotspot2bis ---
test_that("least_square_hotspot2bis returns correct structure for valid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  GC <- 0.5
  f <- 0.2
  result <- least_square_hotspot2bis(sfsWS, sfsSW, GC, f)
  expect_type(result, "list")
  expect_equal(names(result), c("param", "criteria"))
  expect_equal(names(result$param), c("B0", "B1", "mutbias"))
  expect_equal(names(result$criteria), c("SStot", "SSres", "R2", "AIC"))
})

test_that("least_square_hotspot2bis throws error for invalid inputs", {
  sfsWS <- c(1000, 500, 300, 200, 80, 50)
  sfsSW <- c(2000, 800, 400, 150, 50, 10)
  expect_error(least_square_hotspot2bis("abc", sfsSW, 0.5, 0.2), "The two SFSs must have positive numeric values")
  expect_error(least_square_hotspot2bis(sfsWS, sfsSW[-1], 0.5, 0.2), "The two SFSs, WS and SW, must have the same length")
  expect_error(least_square_hotspot2bis(sfsWS, sfsSW, 0, 0.2), "GC content must be strictly between 0 and 1")
  expect_error(least_square_hotspot2bis(sfsWS, sfsSW, 0.5, -0.1), "f must be between 0 and 1/2")
  expect_error(least_square_hotspot2bis(sfsWS, sfsSW, 0.5, 0.6), "f must be between 0 and 1/2")
})

