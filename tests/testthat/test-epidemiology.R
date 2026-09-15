test_that("SIR step conserves population", {

  set.seed(1)

  y <- c(
    S = 50,
    I = 10,
    R = 0
  )

  result <- SIR_step(
    y = y,
    beta = 0.5,
    gamma = 0.1
  )

  expect_true(all(is.finite(result)))

  expect_equal(
    sum(result),
    sum(y)
  )

  expect_true(all(result >= 0))
})


test_that("SIRS step conserves population", {

  set.seed(1)

  y <- c(
    S = 40,
    I = 10,
    R = 20
  )

  result <- SIRS_step(
    y = y,
    beta = 0.5,
    gamma = 0.1,
    omega = 0.05
  )

  expect_true(all(is.finite(result)))

  expect_equal(
    sum(result),
    sum(y)
  )

  expect_true(all(result >= 0))
})


test_that("daytime epidemiology conserves population", {

  set.seed(1)

  roosts <- create_roosts(
    num_roosts = 10,
    num_clusters = 2
  )

  roosts <- initialize_roost_populations(roosts)

  initial_N <- sum(roosts$N)

  result <- run_daytime(
    roost_info = roosts,
    beta = 0.5,
    gamma = 0.1,
    model = "SIR"
  )

  expect_equal(
    sum(result$N),
    initial_N
  )

  expect_true(
    all(result$N == result$S + result$I + result$R)
  )

  expect_true(all(result$N <= result$N_max))
})
