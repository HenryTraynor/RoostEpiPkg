#' @importFrom stats runif rnorm dist
NULL

#' Create a set of roosts
#'
#' Creates spatially clustered roosts with variable carrying capacities.
#'
#' @param num_roosts Number of roosts.
#' @param num_clusters Number of spatial clusters.
#' @param sd Standard deviation of roost locations around cluster centers.
#' @param avg_N_max Mean roost carrying capacity.
#' @param min_N_max Minimum roost carrying capacity.
#'
#' @return A data frame containing roost locations, carrying capacities,
#'   and population compartments.
#'
#' @export
create_roosts <- function(
    num_roosts = 25,
    num_clusters = 5,
    sd = 0.05,
    avg_N_max = 70,
    min_N_max = 10
) {

  # Cluster centers
  cluster_x <- runif(
    num_clusters,
    min = 0,
    max = 1
  )

  cluster_y <- runif(
    num_clusters,
    min = 0,
    max = 1
  )

  # Assign roosts to clusters
  cluster_id <- sample(
    seq_len(num_clusters),
    size = num_roosts,
    replace = TRUE
  )

  # Generate roost locations around cluster centers
  x <- cluster_x[cluster_id] +
    rnorm(num_roosts, mean = 0, sd = sd)

  y <- cluster_y[cluster_id] +
    rnorm(num_roosts, mean = 0, sd = sd)

  # Keep roosts within the landscape
  x <- pmin(pmax(x, 0), 1)
  y <- pmin(pmax(y, 0), 1)

  # Generate maximum roost capacities
  N_max <- pmax(
    round(
      rnorm(
        num_roosts,
        mean = avg_N_max,
        sd = avg_N_max * 0.25
      )
    ),
    min_N_max
  )

  data.frame(
    roost_id = seq_len(num_roosts),
    x = x,
    y = y,
    N_max = N_max,
    S = 0,
    I = 0,
    R = 0,
    N = 0
  )
}


#' Calculate pairwise roost distances
#'
#' Calculates Euclidean distances between all pairs of roosts.
#'
#' @param roost_info Data frame containing `x` and `y` roost coordinates.
#'
#' @return A matrix containing pairwise Euclidean distances between roosts.
#'
#' @export
roost_distance_matrix <- function(roost_info) {

  coords <- as.matrix(
    roost_info[, c("x", "y")]
  )

  as.matrix(
    dist(coords)
  )
}


#' Initialize roost populations
#'
#' Assigns an initial population to each roost and divides individuals
#' between susceptible and infected compartments.
#'
#' @param roost_info Data frame containing roost capacities.
#' @param initial_occupancy Proportion of capacity initially occupied.
#' @param initial_infected Proportion of the initial population that
#'   is infected.
#'
#' @return Updated roost information data frame.
#'
#' @export
initialize_roost_populations <- function(
    roost_info,
    initial_occupancy = 0.25,
    initial_infected = 0.25
) {

  N <- round(
    roost_info$N_max * initial_occupancy
  )

  I <- round(
    N * initial_infected
  )

  S <- N - I

  roost_info$N <- N
  roost_info$S <- S
  roost_info$I <- I
  roost_info$R <- 0

  roost_info
}
