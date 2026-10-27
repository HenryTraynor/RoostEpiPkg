#' @importFrom stats runif rnorm dist
NULL


#' Create a set of roosts
#'
#' Creates spatially clustered roosts with a prescribed proportion of
#' modern and traditional roosts, each with different carrying capacities.
#'
#' @param num_roosts Number of roosts.
#' @param num_clusters Number of spatial clusters.
#' @param sd Standard deviation of roost locations around cluster centers.
#' @param prop_modern Proportion of roosts that are modern, between 0 and 1.
#' @param min_N_max Minimum roost carrying capacity.
#'
#' @return A data frame containing roost locations, types, carrying
#'   capacities, and population compartments.
#'
#' @export
create_roosts <- function(
    num_roosts = 25,
    num_clusters = 5,
    sd = 0.05,
    prop_modern = 0.5,
    min_N_max = 10
) {

  # Validate inputs
  stopifnot(
    num_roosts > 0,
    num_clusters > 0,
    prop_modern >= 0,
    prop_modern <= 1
  )

  # Assign roost types
  n_modern <- round(num_roosts * prop_modern)

  roost_type <- sample(
    c(
      rep("modern", n_modern),
      rep("traditional", num_roosts - n_modern)
    )
  )

  # Set mean carrying capacity by roost type
  mean_N_max <- ifelse(
    roost_type == "modern",
    100,
    20
  )

  # Cluster centers
  cluster_x <- runif(num_clusters, min = 0, max = 1)
  cluster_y <- runif(num_clusters, min = 0, max = 1)

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

  # Generate carrying capacities
  N_max <- pmax(
    round(rnorm(
      num_roosts,
      mean = mean_N_max,
      sd = mean_N_max * 0.25
    )),
    min_N_max
  )

  data.frame(
    roost_id = seq_len(num_roosts),
    roost_type = roost_type,
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

  if(sum(I)==0 & initial_infected != 0) {
    index <- sample(1:length(I),1)
    I[index] <- 1
  }

  S <- N - I

  roost_info$N <- N
  roost_info$S <- S
  roost_info$I <- I
  roost_info$R <- 0

  roost_info
}


#' Assign modernization order
#'
#' Assigns ranked index for progression of modernization
#'
#' @param roost_info Data frame containing roost capacities.
#'
#' @return Updated roost information data frame.
#'
#' @export
assign_modernization_order <- function(roost_info) {

  n <- nrow(roost_info)

  # Random, fixed order for modernization
  modernization_order <- sample.int(n)

  roost_info$modernization_order <- modernization_order

  roost_info
}

#' Increase landscape prop_modern
#'
#' Keeps roost locations but makes more modern
#'
#' @param roost_info Data frame containing roost capacities.
#' @param prop_modern New prop_modern desired
#' @param min_N_max minimum roost size
#'
#' @return Updated roost information data frame.
#'
#' @export

apply_modernization <- function(
    roost_info,
    prop_modern,
    min_N_max = 10
) {

  stopifnot(
    prop_modern >= 0,
    prop_modern <= 1
  )

  n <- nrow(roost_info)

  # Number of roosts that should be modern
  n_modern <- round(prop_modern * n)

  # Preserve the existing coordinates and modernization order
  roost_info$roost_type <- ifelse(
    roost_info$modernization_order <= n_modern,
    "modern",
    "traditional"
  )

  mean_N_max <- ifelse(
    roost_info$roost_type == "modern",
    100,
    20
  )

  roost_info$N_max <- pmax(
    round(rnorm(
      num_roosts,
      mean = mean_N_max,
      sd = mean_N_max * 0.25
    )),
    min_N_max
  )

  # Reset populations before initialization
  roost_info$S <- 0L
  roost_info$I <- 0L
  roost_info$R <- 0L
  roost_info$N <- 0L

  roost_info
}
