test_that("complete simulation conserves population", {

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

  epidemic <- create_epidemic_params(
    beta = 0.5,
    gamma = 0.1,
    model = "SIR"
  )

  sim <- run_simulation(
    roost_info = roosts,
    spatial_params = spatial,
    movement_params = movement,
    epidemic_params = epidemic,
    n_days = 5
  )

  total_N <- sapply(
    sim,
    function(x) sum(x$roosts$N)
  )

  expect_true(
    all(total_N == total_N[1])
  )
})


test_that("simulation maintains population accounting", {

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

  epidemic <- create_epidemic_params(
    beta = 0.5,
    gamma = 0.1
  )

  sim <- run_simulation(
    roost_info = roosts,
    spatial_params = spatial,
    movement_params = movement,
    epidemic_params = epidemic,
    n_days = 5
  )

  for (day in sim) {

    roosts_day <- day$roosts

    expect_true(
      all(
        roosts_day$N ==
          roosts_day$S +
          roosts_day$I +
          roosts_day$R
      )
    )

    expect_true(all(roosts_day$N <= roosts_day$N_max))

    expect_true(all(roosts_day$S >= 0))
    expect_true(all(roosts_day$I >= 0))
    expect_true(all(roosts_day$R >= 0))
  }
})
