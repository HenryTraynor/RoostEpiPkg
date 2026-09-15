#' Create a two-dimensional spatial grid
#'
#' Creates a regular two-dimensional grid and the corresponding
#' ReacTran grid used for spatial diffusion.
#'
#' @param x_min Minimum x-coordinate.
#' @param x_max Maximum x-coordinate.
#' @param y_min Minimum y-coordinate.
#' @param y_max Maximum y-coordinate.
#' @param nx Number of cells in the x direction.
#' @param ny Number of cells in the y direction.
#'
#' @return A list containing the spatial grid, ReacTran grid,
#'   cell dimensions, and grid dimensions.
#'
#' @importFrom ReacTran setup.grid.1D setup.grid.2D
#' @export
create_spatial_grid <- function(
    x_min = 0,
    x_max = 1,
    y_min = 0,
    y_max = 1,
    nx = 100,
    ny = 100
) {
  
  dx <- (x_max - x_min) / nx
  dy <- (y_max - y_min) / ny
  
  # Cell-centered coordinates
  x <- seq(
    x_min + dx / 2,
    x_max - dx / 2,
    length.out = nx
  )
  
  y <- seq(
    y_min + dy / 2,
    y_max - dy / 2,
    length.out = ny
  )
  
  grid <- expand.grid(
    x = x,
    y = y
  )
  
  n_cells <- nrow(grid)
  
  # ReacTran grid
  x_grid <- ReacTran::setup.grid.1D(
    x.up = x_min,
    x.down = x_max,
    N = nx
  )
  
  y_grid <- ReacTran::setup.grid.1D(
    x.up = y_min,
    x.down = y_max,
    N = ny
  )
  
  grid_2D <- ReacTran::setup.grid.2D(
    x.grid = x_grid,
    y.grid = y_grid
  )
  
  list(
    grid = grid,
    grid_2D = grid_2D,
    x = x,
    y = y,
    nx = nx,
    ny = ny,
    dx = dx,
    dy = dy,
    n_cells = n_cells
  )
}