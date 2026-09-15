# movement.R

#' Create initial spatial density
#'
#' Distributes individuals from each roost across the spatial grid
#' using a Gaussian kernel centered on the roost location.
#'
#' @param roost_info Data frame containing roost coordinates and the
#'   requested compartment.
#' @param compartment Name of the compartment column, e.g. "S", "I", or "R".
#' @param grid Data frame containing grid-cell coordinates with columns
#'   `x` and `y`.
#' @param nx Number of grid cells in the x direction.
#' @param ny Number of grid cells in the y direction.
#' @param dx Width of a grid cell in the x direction.
#' @param dy Width of a grid cell in the y direction.
#' @param sigma Standard deviation of the Gaussian spatial kernel.
#'
#' @return A matrix containing the spatial density for the requested
#'   compartment.
#'
#' @export
make_initial_density <- function(
    roost_info,
    compartment,
    grid,
    nx,
    ny,
    dx,
    dy,
    sigma = 0.01
) {

  density <- matrix(
    0,
    nrow = nx,
    ncol = ny
  )

  for (i in seq_len(nrow(roost_info))) {

    distance_squared <-
      (grid$x - roost_info$x[i])^2 +
      (grid$y - roost_info$y[i])^2

    gaussian <- exp(
      -distance_squared / (2 * sigma^2)
    )

    # Normalize so that the integral over the grid equals 1
    gaussian <- gaussian /
      (sum(gaussian) * dx * dy)

    # Weight by the number of individuals in the compartment
    density <- density +
      matrix(
        gaussian * roost_info[[compartment]][i],
        nrow = nx,
        ncol = ny
      )
  }

  density
}


#' Convert a spatial density to cell probabilities
#'
#' Converts a spatial density surface into probabilities that sum to one.
#'
#' @param cell_density Matrix of spatial density values.
#'
#' @return Matrix of probabilities with the same dimensions as
#'   `cell_density`.
#'
#' @export
make_cell_probability <- function(cell_density) {

  cell_density <- pmax(cell_density, 0)

  total <- sum(cell_density)

  if (total == 0) {
    return(
      matrix(
        0,
        nrow = nrow(cell_density),
        ncol = ncol(cell_density)
      )
    )
  }

  cell_density / total
}


#' Run nighttime movement
#'
#' Simulates nighttime movement of susceptible, infected, and recovered
#' individuals using spatial diffusion.
#'
#' @param roost_info Roost population data frame.
#' @param spatial_params Spatial model parameters.
#' @param movement_params Movement model parameters.
#'
#' @return A list containing the final spatial densities for S, I, and R.
#'
#' @importFrom deSolve ode.2D
#' @export
run_night <- function(
    roost_info,
    spatial_params,
    movement_params
) {

  u_S <- make_initial_density(
    roost_info,
    "S",
    spatial_params$grid,
    spatial_params$nx,
    spatial_params$ny,
    spatial_params$dx,
    spatial_params$dy,
    movement_params$sigma
  )

  u_I <- make_initial_density(
    roost_info,
    "I",
    spatial_params$grid,
    spatial_params$nx,
    spatial_params$ny,
    spatial_params$dx,
    spatial_params$dy,
    movement_params$sigma
  )

  u_R <- make_initial_density(
    roost_info,
    "R",
    spatial_params$grid,
    spatial_params$nx,
    spatial_params$ny,
    spatial_params$dx,
    spatial_params$dy,
    movement_params$sigma
  )

  state <- c(
    as.vector(u_S),
    as.vector(u_I),
    as.vector(u_R)
  )

  diffusion_parms <- list(
    nx = spatial_params$nx,
    ny = spatial_params$ny,
    D_x = movement_params$D_x,
    D_y = movement_params$D_y,
    grid_2D = spatial_params$grid_2D
  )

  night <- deSolve::ode.2D(
    y = state,
    times = movement_params$times_night,
    func = diffusion_SIR,
    parms = diffusion_parms,
    dimens = c(
      spatial_params$nx,
      spatial_params$ny
    ),
    lrw = 5000000
  )

  n_cells <- spatial_params$nx * spatial_params$ny

  final_state <- night[nrow(night), ]

  u_S_final <- matrix(
    final_state[2:(n_cells + 1)],
    spatial_params$nx,
    spatial_params$ny
  )

  u_I_final <- matrix(
    final_state[
      (n_cells + 2):(2 * n_cells + 1)
    ],
    spatial_params$nx,
    spatial_params$ny
  )

  u_R_final <- matrix(
    final_state[
      (2 * n_cells + 2):(3 * n_cells + 1)
    ],
    spatial_params$nx,
    spatial_params$ny
  )

  list(
    u_S = u_S_final,
    u_I = u_I_final,
    u_R = u_R_final
  )
}
