#' Create spatial model parameters
#'
#' @param x_min Minimum x coordinate.
#' @param x_max Maximum x coordinate.
#' @param y_min Minimum y coordinate.
#' @param y_max Maximum y coordinate.
#' @param nx Number of grid cells in the x direction.
#' @param ny Number of grid cells in the y direction.
#'
#' @return A list containing the spatial grid and associated parameters.
#' @export
create_spatial_params <- function(
    x_min = 0,
    x_max = 1,
    y_min = 0,
    y_max = 1,
    nx = 100,
    ny = 100
) {

  spatial <- create_spatial_grid(
    x_min = x_min,
    x_max = x_max,
    y_min = y_min,
    y_max = y_max,
    nx = nx,
    ny = ny
  )

  spatial
}


#' Create movement parameters
#'
#' @param lambda Distance-decay parameter.
#' @param sigma Initial spatial density spread around roosts.
#' @param D_base Baseline diffusion coefficient.
#' @param gamma_D Controls the effect of roost isolation on diffusion.
#' @param alpha Controls population-based roost attraction.
#' @param epsilon Small value preventing zero attraction.
#' @param times_night Time sequence for the nighttime diffusion simulation.
#'
#' @return A list of movement parameters.
#' @export
create_movement_params <- function(
    lambda = 5,
    sigma = 0.01,
    D_base = 0.0005,
    gamma_D = 1,
    alpha = 1,
    epsilon = 0.01,
    times_night = seq(0, 10, 0.1)
) {

  list(
    lambda = lambda,
    sigma = sigma,
    D_base = D_base,
    gamma_D = gamma_D,
    alpha = alpha,
    epsilon = epsilon,
    times_night = times_night
  )
}


#' Create epidemic parameters
#'
#' @param beta Transmission rate.
#' @param gamma Recovery rate.
#' @param model Epidemic model, either `"SIR"` or `"SIRS"`.
#' @param omega Rate of loss of immunity for the SIRS model.
#'
#' @return A list of epidemic parameters.
#' @export
create_epidemic_params <- function(
    beta = 0.5,
    gamma = 0.1,
    model = "SIR",
    omega = NULL
) {

  model <- match.arg(model, c("SIR", "SIRS"))

  if (model == "SIRS" && is.null(omega)) {
    stop("omega must be supplied when model = 'SIRS'.")
  }

  list(
    beta = beta,
    gamma = gamma,
    model = model,
    omega = omega
  )
}

#' Set up spatial movement parameters
#'
#' Calculates roost connectivity, diffusion coefficients, the diffusion
#' field, and cell-to-roost distances.
#'
#' @param roost_info Roost information.
#' @param spatial_params Spatial model parameters.
#' @param lambda Distance-decay parameter.
#' @param sigma Initial spatial density spread.
#' @param D_base Baseline diffusion coefficient.
#' @param gamma_D Effect of roost isolation on diffusion.
#' @param alpha Population-attraction exponent.
#' @param epsilon Small value preventing zero attraction.
#' @param times_night Nighttime integration times.
#'
#' @return A list containing movement parameters and derived spatial objects.
#' @export
setup_movement <- function(
    roost_info,
    spatial_params,
    lambda = 5,
    sigma = 0.01,
    D_base = 0.0005,
    gamma_D = 1,
    alpha = 1,
    epsilon = 0.01,
    times_night = seq(0, 10, 0.1)
) {

  # Roost distances
  dist_matrix <- roost_distance_matrix(roost_info)

  cell_roost_dist <- cell_roost_distances(
    grid = spatial_params$grid,
    roost_info = roost_info
  )

  # Connectivity and roost-specific diffusion
  diffusion_coefficients <- calculate_diffusion_coefficients(
    distance_matrix = dist_matrix,
    lambda = lambda,
    D_base = D_base,
    gamma_D = gamma_D
  )

  # Continuous diffusion field
  D_field_vector <- make_diffusion_field(
    cell_roost_distances = cell_roost_dist,
    D_roost = diffusion_coefficients$D_roost,
    lambda = lambda
  )

  D_field <- matrix(
    D_field_vector,
    nrow = spatial_params$nx,
    ncol = spatial_params$ny
  )

  # Diffusion coefficients at cell faces
  faces <- make_diffusion_faces(D_field)

  list(
    lambda = lambda,
    sigma = sigma,
    D_base = D_base,
    gamma_D = gamma_D,
    alpha = alpha,
    epsilon = epsilon,
    times_night = times_night,

    dist_matrix = dist_matrix,
    cell_roost_distances = cell_roost_dist,

    connectivity_matrix =
      diffusion_coefficients$connectivity_matrix,

    connectivity =
      diffusion_coefficients$connectivity,

    isolation =
      diffusion_coefficients$isolation,

    D_roost =
      diffusion_coefficients$D_roost,

    D_field = D_field,
    D_x = faces$D_x,
    D_y = faces$D_y
  )
}
