test_that("movement setup produces finite values", {

  set.seed(1)

  roosts <- create_roosts(
    num_roosts = 10,
    num_clusters = 2
  )

  roosts <- initialize_roost_populations(roosts)

  spatial <- create_spatial_params(
    nx = 25,
    ny = 25
  )

  movement <- setup_movement(
    roost_info = roosts,
    spatial_params = spatial
  )

  expect_true(all(is.finite(movement$dist_matrix)))
  expect_true(all(is.finite(movement$cell_roost_distances)))

  expect_true(all(is.finite(movement$D_roost)))
  expect_true(all(is.finite(movement$D_field)))
  expect_true(all(is.finite(movement$D_x)))
  expect_true(all(is.finite(movement$D_y)))
})


test_that("nighttime diffusion conserves population", {

  set.seed(1)

  roosts <- create_roosts(
    num_roosts = 10,
    num_clusters = 2
  )

  roosts <- initialize_roost_populations(roosts)

  spatial <- create_spatial_params(
    nx = 25,
    ny = 25
  )

  movement <- setup_movement(
    roost_info = roosts,
    spatial_params = spatial,
    times_night = seq(0, 2, 0.2)
  )

  night <- run_night(
    roost_info = roosts,
    spatial_params = spatial,
    movement_params = movement
  )

  N_start <- sum(roosts$N)

  N_S <- sum(night$u_S) * spatial$dx * spatial$dy
  N_I <- sum(night$u_I) * spatial$dx * spatial$dy
  N_R <- sum(night$u_R) * spatial$dx * spatial$dy

  N_final <- N_S + N_I + N_R

  expect_equal(
    N_final,
    N_start,
    tolerance = 1e-6
  )

  expect_true(all(night$u_S >= 0))
  expect_true(all(night$u_I >= 0))
  expect_true(all(night$u_R >= 0))
})
