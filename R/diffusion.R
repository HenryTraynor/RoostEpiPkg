# diffusion.R

#' Calculate roost-specific diffusion coefficients
#'
#' Calculates roost connectivity, isolation, and diffusion coefficients
#' based on the distance-decay kernel between roosts.
#'
#' @param distance_matrix Matrix of pairwise distances between roosts.
#' @param lambda Distance-decay parameter.
#' @param D_base Baseline diffusion coefficient.
#' @param gamma_D Strength of the effect of isolation on diffusion.
#'
#' @return A list containing connectivity, isolation, and diffusion
#'   coefficients for each roost.
#'
#' @export
calculate_diffusion_coefficients <- function(
    distance_matrix,
    lambda = 5,
    D_base = 0.0005,
    gamma_D = 1
) {
  
  # Distance-decay connectivity between roosts
  connectivity_matrix <- exp(-lambda * distance_matrix)
  
  # Roosts do not contribute to their own connectivity
  diag(connectivity_matrix) <- 0
  
  # Total connectivity of each roost
  connectivity <- rowSums(connectivity_matrix)
  
  # Isolation is the inverse of connectivity
  isolation <- 1 / connectivity
  
  # Diffusion decreases with isolation
  D_roost <- D_base * exp(-gamma_D * isolation)
  
  list(
    connectivity_matrix = connectivity_matrix,
    connectivity = connectivity,
    isolation = isolation,
    D_roost = D_roost
  )
}


#' Create a spatially varying diffusion field
#'
#' Calculates the diffusion coefficient at each spatial grid cell by
#' weighting roost-specific diffusion coefficients according to the
#' distance-decay relationship between cells and roosts.
#'
#' @param cell_roost_distances Matrix of distances between grid cells
#'   and roosts.
#' @param D_roost Vector of roost-specific diffusion coefficients.
#' @param lambda Distance-decay parameter.
#'
#' @return A vector containing the diffusion coefficient for each grid cell.
#'
#' @export
make_diffusion_field <- function(
    cell_roost_distances,
    D_roost,
    lambda = 5
) {
  
  distance_decay <- exp(-lambda * cell_roost_distances)
  
  weights <- rowSums(distance_decay)
  
  D_field <- rowSums(
    distance_decay *
      matrix(
        D_roost,
        nrow = nrow(distance_decay),
        ncol = length(D_roost),
        byrow = TRUE
      )
  ) / weights
  
  D_field
}


#' Create diffusion coefficients at grid-cell faces
#'
#' Converts a cell-centered diffusion field into x- and y-direction
#' diffusion coefficients required by ReacTran.
#'
#' @param D_field Matrix of diffusion coefficients at grid-cell centers.
#'
#' @return A list containing x- and y-direction face coefficients.
#'
#' @export
make_diffusion_faces <- function(D_field) {
  
  nx <- nrow(D_field)
  ny <- ncol(D_field)
  
  D_x <- matrix(
    NA_real_,
    nrow = nx + 1,
    ncol = ny
  )
  
  D_x[2:nx, ] <-
    (D_field[1:(nx - 1), ] +
       D_field[2:nx, ]) / 2
  
  # Zero-gradient boundaries
  D_x[1, ] <- D_field[1, ]
  D_x[nx + 1, ] <- D_field[nx, ]
  
  
  D_y <- matrix(
    NA_real_,
    nrow = nx,
    ncol = ny + 1
  )
  
  D_y[, 2:ny] <-
    (D_field[, 1:(ny - 1)] +
       D_field[, 2:ny]) / 2
  
  # Zero-gradient boundaries
  D_y[, 1] <- D_field[, 1]
  D_y[, ny + 1] <- D_field[, ny]
  
  
  list(
    D_x = D_x,
    D_y = D_y
  )
}


#' Calculate diffusion derivatives for SIR compartments
#'
#' Calculates spatial diffusion of susceptible, infected, and recovered
#' individuals using ReacTran.
#'
#' @param t Current simulation time.
#' @param state Vector containing S, I, and R spatial densities.
#' @param parms List containing nx, ny, D_x, D_y, and grid_2D.
#'
#' @return A list containing the derivatives of S, I, and R.
#'
#' @importFrom ReacTran tran.2D
#' @export
diffusion_SIR <- function(t, state, parms) {
  
  nx <- parms$nx
  ny <- parms$ny
  D_x <- parms$D_x
  D_y <- parms$D_y
  grid_2D <- parms$grid_2D
  
  n_cells <- nx * ny
  
  u_S <- matrix(
    state[1:n_cells],
    nrow = nx,
    ncol = ny
  )
  
  u_I <- matrix(
    state[(n_cells + 1):(2 * n_cells)],
    nrow = nx,
    ncol = ny
  )
  
  u_R <- matrix(
    state[(2 * n_cells + 1):(3 * n_cells)],
    nrow = nx,
    ncol = ny
  )
  
  
  dS <- ReacTran::tran.2D(
    C = u_S,
    D.x = D_x,
    D.y = D_y,
    v.x = 0,
    v.y = 0,
    grid = grid_2D
  )
  
  dI <- ReacTran::tran.2D(
    C = u_I,
    D.x = D_x,
    D.y = D_y,
    v.x = 0,
    v.y = 0,
    grid = grid_2D
  )
  
  dR <- ReacTran::tran.2D(
    C = u_R,
    D.x = D_x,
    D.y = D_y,
    v.x = 0,
    v.y = 0,
    grid = grid_2D
  )
  
  
  list(
    c(
      dS$dC,
      dI$dC,
      dR$dC
    )
  )
}