# Test of the functions in derivatives.R



# --- Tests for d_expected_log_ratio ---
test_that("d_expected_log_ratio returns correct structure for valid inputs", {
  WS <- c(100, 50, 30)
  SW <- c(200, 80, 30)
  e1 <- 0.01
  e2 <- 0.02
  result <- d_expected_log_ratio(WS, SW, e1, e2)
  expect_type(result, "list")
  expect_equal(names(result), c("d1", "d2"))
  expect_type(result$d1, "double")
  expect_type(result$d2, "double")
  expect_length(result$d1, length(WS))
  expect_length(result$d2, length(WS))
})

test_that("d_expected_log_ratio throws error for invalid inputs", {
  WS <- c(100, 50, 30)
  SW <- c(200, 80, 30)
  expect_error(d_expected_log_ratio("abc", SW, 0.01, 0.02), "The two SFSs must have positive numeric values")
  expect_error(d_expected_log_ratio(WS, SW[-1], 0.01, 0.02), "The two SFSs, WS and SW, must have the same length")
  expect_error(d_expected_log_ratio(c(-100, 50, 30), SW, 0.01, 0.02), "The two SFSs must have positive numeric values")
  expect_error(d_expected_log_ratio(WS, SW, -0.1, 0.02), "e1 content must be between 0 and 1")
  expect_error(d_expected_log_ratio(WS, SW, 1.1, 0.02), "e1 content must be between 0 and 1")
  expect_error(d_expected_log_ratio(WS, SW, 0.01, -0.1), "e2 content must be between 0 and 1")
  expect_error(d_expected_log_ratio(WS, SW, 0.01, 1.1), "e2 content must be between 0 and 1")
})

# --- Tests for ratio_B ---
test_that("ratio_B returns correct value for valid inputs", {
  B <- 1.0
  x <- c(0.1, 0.5, 0.9)
  result <- ratio_B(B, x)
  expect_type(result, "double")
  expect_length(result, length(x))
  # Test for B=0
  expect_equal(ratio_B(0, x), 1 - x, tolerance = 1e-6)
})

test_that("ratio_B throws error for invalid inputs", {
  expect_error(ratio_B("abc", 0.5), "B must be a numerical value")
  expect_error(ratio_B(1.0, -0.1), "x must be between 0 and 1")
  expect_error(ratio_B(1.0, 1.1), "x must be between 0 and 1")
})

# --- Tests for d_ratio_B ---
test_that("d_ratio_B returns correct value for valid inputs", {
  B <- 1.0
  x <- c(0.1, 0.5, 0.9)
  result <- d_ratio_B(B, x)
  expect_type(result, "double")
  expect_length(result, length(x))
  # Test for B=0
  expected <- x * (1 - x) / 2
  expect_equal(d_ratio_B(0, x), expected, tolerance = 1e-6)
})

test_that("d_ratio_B throws error for invalid inputs", {
  expect_error(d_ratio_B("abc", 0.5), "B must be a numerical value")
  expect_error(d_ratio_B(1.0, -0.1), "x must be between 0 and 1")
  expect_error(d_ratio_B(1.0, 1.1), "x must be between 0 and 1")
})

# --- Tests for d_hotspot1 ---
test_that("d_hotspot1 returns correct structure for valid inputs", {
  B <- 1.0
  f <- 0.2
  x <- c(0.1, 0.5, 0.9)
  result <- d_hotspot1(B, f, x)
  expect_type(result, "list")
  expect_equal(names(result), c("dB", "df"))
  expect_type(result$dB, "double")
  expect_type(result$df, "double")
  expect_length(result$dB, length(x))
  expect_length(result$df, length(x))
})

test_that("d_hotspot1 throws error for invalid inputs", {
  expect_error(d_hotspot1("abc", 0.2, 0.5), "B must be a numerical value")
  expect_error(d_hotspot1(1.0, -0.1, 0.5), "f must be between 0 and 1")
  expect_error(d_hotspot1(1.0, 1.1, 0.5), "f must be between 0 and 1")
  expect_error(d_hotspot1(1.0, 0.2, -0.1), "x must be between 0 and 1")
  expect_error(d_hotspot1(1.0, 0.2, 1.1), "x must be between 0 and 1")
})

# --- Tests for d_hotspot2 ---
test_that("d_hotspot2 returns correct structure for valid inputs", {
  B0 <- 0.1
  B1 <- 5.0
  f <- 0.2
  x <- c(0.1, 0.5, 0.9)
  result <- d_hotspot2(B0, B1, f, x)
  expect_type(result, "list")
  expect_equal(names(result), c("dB0", "dB1", "df"))
  expect_type(result$dB0, "double")
  expect_type(result$dB1, "double")
  expect_type(result$df, "double")
  expect_length(result$dB0, length(x))
  expect_length(result$dB1, length(x))
  expect_length(result$df, length(x))
})

test_that("d_hotspot2 throws error for invalid inputs", {
  expect_error(d_hotspot2("abc", 5.0, 0.2, 0.5), "B0 must be a numerical value")
  expect_error(d_hotspot2(0.1, "abc", 0.2, 0.5), "B1 must be a numerical value")
  expect_error(d_hotspot2(0.1, 5.0, -0.1, 0.5), "f must be between 0 and 1")
  expect_error(d_hotspot2(0.1, 5.0, 1.1, 0.5), "f must be between 0 and 1")
  expect_error(d_hotspot2(0.1, 5.0, 0.2, -0.1), "x must be between 0 and 1")
  expect_error(d_hotspot2(0.1, 5.0, 0.2, 1.1), "x must be between 0 and 1")
})

