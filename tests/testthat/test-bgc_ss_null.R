# Test of the function in bgc_ss_null.R


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


# --- Tests for sum_of_squares_NULL_GC ---
test_that("sum_of_squares_NULL_GC returns correct structure for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  GC <- 0.5
  result <- sum_of_squares_NULL_GC(sfsWS, sfsSW,GC)
  expect_type(result, "list")
  expect_equal(names(result), c("SStot"))
  expect_type(result$SStot, "double")
})

test_that("sum_of_squares_NULL_GC throws error for invalid inputs", {
  expect_error(sum_of_squares_NULL_GC("abc", c(200, 80, 30, 10, 5),0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_NULL_GC(c(100, 50, 30, 15, 10), c(200, 80, 30, 10),0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_NULL_GC(c(-100, 50, 30, 15, 10), c(200, 80, 30, 10, 5),0.5), "The two SFSs must have positive numeric values")
})


# --- Tests for sum_of_squares_NULL_GC ---
test_that("sum_of_squares_NULL_GC_err returns correct structure for valid inputs", {
  sfsWS <- c(100, 50, 30, 15, 10)
  sfsSW <- c(200, 80, 30, 10, 5)
  GC <- 0.5
  param <- c(0.01,0.02)
  result <- sum_of_squares_NULL_GC_err(param,sfsWS, sfsSW,GC)
  expect_type(result, "double")
})

test_that("sum_of_squares_NULL_GC_err throws error for invalid inputs", {
  expect_error(sum_of_squares_NULL_GC_err(c(0.01,0.01),"abc", c(200, 80, 30, 10, 5),0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_NULL_GC_err(c(0.01,0.01),c(100, 50, 30, 15, 10), c(200, 80, 30, 10),0.5), "The two SFSs, WS and SW, must have the same length")
  expect_error(sum_of_squares_NULL_GC_err(c(0.01,0.01),c(-100, 50, 30, 15, 10), c(200, 80, 30, 10, 5),0.5), "The two SFSs must have positive numeric values")
  expect_error(sum_of_squares_NULL_GC_err(c(0.01),c(100, 50, 30, 15, 10), c(200, 80, 30, 10, 5),0.5), "A two values vector must be given as par")
})


