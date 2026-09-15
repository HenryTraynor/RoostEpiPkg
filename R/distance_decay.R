#' Calculate distances from grid cells to roosts
#'
#' @param grid Spatial grid containing x and y coordinates.
#' @param roost_info Data frame containing roost coordinates.
#'
#' @return A matrix of distances between grid cells and roosts.
cell_roost_distances <- function(grid, roost_info) {

  n_cells <- nrow(grid)
  n_roosts <- nrow(roost_info)

  cell_roost_dist <- matrix(
    NA_real_,
    nrow = n_cells,
    ncol = n_roosts
  )

  for (j in seq_len(n_roosts)) {

    cell_roost_dist[, j] <-
      sqrt(
        (grid$x - roost_info$x[j])^2 +
          (grid$y - roost_info$y[j])^2
      )
  }

  return(cell_roost_dist)
}

#' Calculate exponential distance decay
#'
#' @param distance Distance values.
#' @param lambda Distance-decay parameter.
#'
#' @return Distance-decay weights.
distance_decay <- function(distance, lambda) {
  exp(-lambda * distance)
}

#' Calculate roost connectivity
#'
#' @param dist_matrix Matrix of distances between roosts.
#' @param lambda Distance-decay parameter.
#'
#' @return A list containing connectivity information.
roost_connectivity <- function(dist_matrix, lambda) {

  connectivity_matrix <-
    exp(-lambda * dist_matrix)

  # A roost is not connected to itself
  diag(connectivity_matrix) <- 0

  connectivity <- rowSums(connectivity_matrix)

  isolation <- 1 / connectivity

  return(
    list(
      connectivity_matrix = connectivity_matrix,
      connectivity = connectivity,
      isolation = isolation
    )
  )
}
