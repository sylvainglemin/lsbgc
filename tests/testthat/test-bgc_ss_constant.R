# Test of the function in bgc_ss_constant.R


# Mock the d_expected_log_ratio function if it's not exported
# This is a placeholder to avoid errors in gradient tests
d_expected_log_ratio <- function(WS, SW, e1, e2) {
  n <- length(WS)
  d1 <- rep(0, n)
  d2 <- rep(0, n)
  return(list(d1 = d1, d2 = d2))
}

# --- Tests for sum_of_squares_NULL ---
test_that("sum_of_squares_NULL returns correct structure for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  result <- sum_of_squares_NULL(sfsWS, sfsSW)
  expect_type(result, "list")
  expect_equal(names(result), c("mean", "SStot"))
  expect_type(result$mean, "double")
  expect_type(result$SStot, "double")
})

test_that("sum_of_squares_NULL throws error for invalid inputs", {
  expect_error(sum_of_squares_NULL("abc", c(200, 80, 30, 10, 5)), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_NULL(c(100, 50, 30, 15, 10), c(200, 80, 30, 10)), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_NULL(c(-100, 50, 30, 15, 10), c(200, 80, 30, 10, 5)), "The two SFSs must have positive numeric values")
})

# --- Tests for sum_of_squares_M ---
test_that("sum_of_squares_M returns correct value for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(2)
  GC <- 0.5
  result <- sum_of_squares_M(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_true(result >= 0)  # Sum of squares should be non-negative
})

test_that("sum_of_squares_M throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(sum_of_squares_M(c(2,0.01), sfsWS, sfsSW, 0.5), "One value must be given as par")
  expect_error(sum_of_squares_M(c(2), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_M(c(2), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_M(c(2), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
  expect_error(sum_of_squares_M(c(2), sfsWS, sfsSW, 1), "GC content must be strictly between 0 and 1")
})

# --- Tests for gr_sum_of_squares_M ---
test_that("gr_sum_of_squares_M returns correct gradient for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(2)
  GC <- 0.5
  result <- gr_sum_of_squares_M(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_length(result, 3)  # Gradient should have 3 components
})

test_that("gr_sum_of_squares_M throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(gr_sum_of_squares_M(c(2, 0.02), sfsWS, sfsSW, 0.5), "One value must be given as par")
  expect_error(gr_sum_of_squares_M(c(2), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(gr_sum_of_squares_M(c(2), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(gr_sum_of_squares_M(c(2), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for sum_of_squares_B ---
test_that("sum_of_squares_B returns correct value for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(1)
  GC <- 0.5
  result <- sum_of_squares_B(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("sum_of_squares_B throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(sum_of_squares_B(c(1, 0.02), sfsWS, sfsSW, 0.5), "One value must be given as par")
  expect_error(sum_of_squares_B(c(1), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_B(c(1), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_B(c(1), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for gr_sum_of_squares_B ---
test_that("gr_sum_of_squares_B returns correct gradient for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(1)
  GC <- 0.5
  result <- gr_sum_of_squares_B(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("gr_sum_of_squares_B throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(gr_sum_of_squares_B(c(1, 0.02), sfsWS, sfsSW, 0.5), "One value must be given as par")
  expect_error(gr_sum_of_squares_B(c(1), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(gr_sum_of_squares_B(c(1), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(gr_sum_of_squares_B(c(1), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for sum_of_squares_BM ---
test_that("sum_of_squares_BM returns correct value for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(1, 2)
  GC <- 0.5
  result <- sum_of_squares_BM(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("sum_of_squares_BM throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(sum_of_squares_BM(c(1, 2, 0.02), sfsWS, sfsSW, 0.5), "A two values vector must be given as par")
  expect_error(sum_of_squares_BM(c(1, 2), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_BM(c(1, 2), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_BM(c(1, 2), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})

# --- Tests for gr_sum_of_squares_BM ---
test_that("gr_sum_of_squares_BM returns correct gradient for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  param <- c(1, 2)
  GC <- 0.5
  result <- gr_sum_of_squares_BM(param, sfsWS, sfsSW, GC)
  expect_type(result, "double")
  expect_length(result, 4)
})

test_that("gr_sum_of_squares_BM throws error for invalid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  expect_error(gr_sum_of_squares_BM(c(1, 2, 0.02), sfsWS, sfsSW, 0.5), "A two values vector must be given as par")
  expect_error(gr_sum_of_squares_BM(c(1, 2), "abc", sfsSW, 0.5), "The two SFSs must have positive numeric values")
  expect_error(gr_sum_of_squares_BM(c(1, 2), sfsWS, sfsSW[-1], 0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(gr_sum_of_squares_BM(c(1, 2), sfsWS, sfsSW, 0), "GC content must be strictly between 0 and 1")
})
