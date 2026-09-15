test_that("spatial grid has correct dimensions", {

  spatial <- create_spatial_params(
    x_min = 0,
    x_max = 1,
    y_min = 0,
    y_max = 1,
    nx = 50,
    ny = 40
  )

  expect_equal(spatial$nx, 50)
  expect_equal(spatial$ny, 40)

  expect_equal(nrow(spatial$grid), 50 * 40)

  expect_equal(length(spatial$x), 50)
  expect_equal(length(spatial$y), 40)

  expect_true(spatial$dx > 0)
  expect_true(spatial$dy > 0)

  expect_true(is.list(spatial$grid_2D))
})
