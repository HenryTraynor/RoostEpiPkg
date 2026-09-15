test_that("create_roosts creates valid roosts", {

  set.seed(1)

  roosts <- create_roosts(
    num_roosts = 25,
    num_clusters = 5,
    sd = 0.05,
    avg_N_max = 70,
    min_N_max = 10
  )

  expect_equal(nrow(roosts), 25)

  expect_true(all(roosts$roost_id == 1:25))

  expect_true(all(roosts$x >= 0 & roosts$x <= 1))
  expect_true(all(roosts$y >= 0 & roosts$y <= 1))

  expect_true(all(roosts$N_max >= 10))

  expect_true(all(roosts$S == 0))
  expect_true(all(roosts$I == 0))
  expect_true(all(roosts$R == 0))
  expect_true(all(roosts$N == 0))
})


test_that("initialize_roost_populations creates valid populations", {

  set.seed(1)

  roosts <- create_roosts(
    num_roosts = 25,
    num_clusters = 5
  )

  roosts <- initialize_roost_populations(
    roosts,
    initial_occupancy = 0.25,
    initial_infected = 0.25
  )

  expect_true(all(roosts$N <= roosts$N_max))

  expect_true(all(roosts$S >= 0))
  expect_true(all(roosts$I >= 0))
  expect_true(all(roosts$R >= 0))

  expect_true(
    all(roosts$N == roosts$S + roosts$I + roosts$R)
  )
})
