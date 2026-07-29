# Tests of the functions in utils.R

# --- Tests for hn ---
test_that("hn returns correct harmonic numbers", {
  expect_equal(hn(1), 1, tolerance = 1e-6)
  expect_equal(hn(2), 1.5, tolerance = 1e-6)
  expect_equal(hn(10), sum(1/(1:10)), tolerance = 1e-6)
})

test_that("hn throws error for invalid inputs", {
  expect_error(hn(0), "n must be a positive integer")
  expect_error(hn(-5), "n must be a positive integer")
  expect_error(hn(3.5), "n must be a positive integer")
  expect_error(hn("abc"), "n must be a positive integer")
})

# --- Tests for skewness_sfs ---
test_that("skewness_sfs returns correct skewness and p-value", {
  sfs <- c(100, 50, 30, 20, 30, 40, 80)
  result <- skewness_sfs(sfs)
  expect_type(result, "list")
  expect_equal(names(result), c("skewness", "p.value"))
  expect_type(result$skewness, "double")
  expect_type(result$p.value, "double")
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("skewness_sfs throws error for invalid inputs", {
  expect_error(skewness_sfs(c(100)), "sfs must be a numeric vector of size greater or equal to 2")
  expect_error(skewness_sfs("abc"), "sfs must be a numeric vector of size greater or equal to 2")
})

# --- Tests for project_sfs ---
test_that("project_sfs returns correct projected SFS", {
  sfs <- c(100, 50, 30, 20, 30, 40, 80)
  sfs_proj <- c(98.75000,47.85714,44.28571,80.89286)
  projected <- project_sfs(sfs, 5)
  expect_type(projected, "double")
  expect_length(projected, 4)
  expect_equal(projected, sfs_proj, tolerance = 1e-6)
})

test_that("project_sfs throws error for invalid inputs", {
  sfs <- c(100, 50, 30, 20, 30, 40, 80)
  expect_error(project_sfs(sfs, 9), "m must be a positive integer lower or equal to the length of SFS + 1",fixed=T)
  expect_error(project_sfs(sfs, 3.5), "m must be a positive integer lower or equal to the length of SFS + 1",fixed=T)
  expect_error(project_sfs("abc", 5), "sfs must be a numeric vector of size greater or equal to 2")
})

test_that("project_sfs returns same SFS when m = n", {
  sfs <- c(100, 50, 30, 20, 30, 40, 80)
  n <- length(sfs) + 1
  projected <- project_sfs(sfs, n)
  expect_equal(projected, sfs, tolerance = 1e-6)
})

