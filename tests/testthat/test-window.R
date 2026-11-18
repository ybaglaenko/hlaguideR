test_that("the editing window is inclusive at both bounds", {
  d <- data.frame(position = c(2L, 3L, 10L, 11L))
  expect_equal(nrow(filter_editing_window(d, 3, 10)), 2L)
})

test_that("missing positions are kept unless asked otherwise", {
  d <- data.frame(position = c(5L, NA_integer_))
  expect_equal(nrow(filter_editing_window(d, 3, 10)), 2L)
  expect_equal(nrow(filter_editing_window(d, 3, 10, keep_na = FALSE)), 1L)
})

test_that("a missing column is an error, not a silent pass", {
  expect_error(filter_editing_window(data.frame(x = 1), 3, 10), "not found")
})
