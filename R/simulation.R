# simulation.R


#' Run one complete model day
#'
#' @param roost_info Roost population data frame.
#' @param spatial_params Spatial model parameters.
#' @param movement_params Movement model parameters.
#' @param epidemic_params Epidemic model parameters.
#'
#' @return A list containing roost populations, nighttime densities,
#'   and dawn return probabilities.
#' @export
run_day <- function(
    roost_info,
    spatial_params,
    movement_params,
    epidemic_params
) {

  # ---- Daytime epidemiology ----

  roosts_after_day <- run_daytime(
    roost_info = roost_info,
    beta = epidemic_params$beta,
    gamma = epidemic_params$gamma,
    model = epidemic_params$model,
    omega = epidemic_params$omega
  )


  # ---- Nighttime movement ----

  night_output <- run_night(
    roost_info = roosts_after_day,
    spatial_params = spatial_params,
    movement_params = movement_params
  )


  # ---- Dawn return ----

  dawn_output <- run_dawn_return(
    night_output = night_output,
    roost_info = roosts_after_day,
    cell_roost_distances = movement_params$cell_roost_distances,
    alpha = movement_params$alpha,
    epsilon = movement_params$epsilon,
    lambda = movement_params$lambda
  )


  list(
    roosts = dawn_output$roosts,
    night = night_output,
    roost_probability = dawn_output$roost_probability
  )
}


#' Run a multi-day simulation
#'
#' @param roost_info Initial roost populations.
#' @param spatial_params Spatial model parameters.
#' @param movement_params Movement model parameters.
#' @param epidemic_params Epidemic model parameters.
#' @param n_days Number of simulation days.
#'
#' @return A list containing model output for each day.
#' @export
run_simulation <- function(
    roost_info,
    spatial_params,
    movement_params,
    epidemic_params,
    n_days = 30
) {

  simulation <- vector("list", n_days + 1)

  simulation[[1]] <- list(
    roosts = roost_info
  )

  for (day in seq_len(n_days)) {

    simulation[[day + 1]] <- run_day(
      roost_info = simulation[[day]]$roosts,
      spatial_params = spatial_params,
      movement_params = movement_params,
      epidemic_params = epidemic_params
    )
  }

  names(simulation) <- paste0("day_", 0:n_days)

  simulation
}
