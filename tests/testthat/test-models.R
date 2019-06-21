context("models")

skip_on_cran()

res1 <- 3

test_that("Expected classes", {
  expect_equal(class(res1), "numeric")
})
