test_that("dawn return conserves population and respects capacity", {

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

  result <- run_dawn_return(
    night_output = night,
    roost_info = roosts,
    cell_roost_distances = movement$cell_roost_distances,
    alpha = movement$alpha,
    epsilon = movement$epsilon,
    lambda = movement$lambda
  )

  new_roosts <- result$roosts

  expect_equal(
    sum(new_roosts$N),
    sum(roosts$N)
  )

  expect_true(
    all(
      new_roosts$N ==
        new_roosts$S +
        new_roosts$I +
        new_roosts$R
    )
  )

  expect_true(
    all(new_roosts$N <= new_roosts$N_max)
  )

  expect_true(all(new_roosts$S >= 0))
  expect_true(all(new_roosts$I >= 0))
  expect_true(all(new_roosts$R >= 0))
})
