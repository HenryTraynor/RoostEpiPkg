# dawn_return.R


#' Assign individuals to roosts subject to capacity constraints
#'
#' Assigns individuals to roosts based on cell-specific roost probabilities
#' while preventing roost capacities from being exceeded.
#'
#' @param n_bats Number of individuals to assign.
#' @param cell_probability Probability of selecting each spatial cell.
#' @param roost_probability Matrix giving roost probabilities for each
#'   spatial cell.
#' @param N_current Current number of individuals occupying each roost.
#' @param N_max Maximum capacity of each roost.
#'
#' @return Integer vector containing the number of newly assigned
#'   individuals at each roost.
assign_bats <- function(
    n_bats,
    cell_probability,
    roost_probability,
    N_current,
    N_max
) {
  
  n_roosts <- length(N_current)
  n_cells <- length(cell_probability)
  
  assignments <- integer(n_roosts)
  
  if (n_bats <= 0) {
    return(assignments)
  }
  
  # Remaining capacity of each roost
  capacity_remaining <- pmax(
    N_max - N_current,
    0
  )
  
  # Remove any invalid cell probabilities
  cell_probability <- pmax(cell_probability, 0)
  
  if (sum(cell_probability) <= 0) {
    stop("cell_probability must contain positive probability.")
  }
  
  cell_probability <- cell_probability / sum(cell_probability)
  
  # Select a spatial cell for every bat
  selected_cells <- sample(
    seq_len(n_cells),
    size = n_bats,
    replace = TRUE,
    prob = cell_probability
  )
  
  for (cell in selected_cells) {
    
    # Roosts that still have capacity
    available <- which(capacity_remaining > 0)
    
    if (length(available) == 0) {
      break
    }
    
    probabilities <- roost_probability[
      cell,
      available,
      drop = TRUE
    ]
    
    probabilities <- pmax(probabilities, 0)
    
    if (sum(probabilities) <= 0) {
      # If no available roost has positive probability,
      # this individual cannot be assigned.
      next
    }
    
    probabilities <- probabilities / sum(probabilities)
    
    selected_roost <- sample(
      available,
      size = 1,
      prob = probabilities
    )
    
    assignments[selected_roost] <-
      assignments[selected_roost] + 1
    
    capacity_remaining[selected_roost] <-
      capacity_remaining[selected_roost] - 1
  }
  
  assignments
}


#' Calculate roost attraction probabilities
#'
#' Calculates the probability that an individual in each spatial cell
#' returns to each roost based on distance, current roost population,
#' and available roost capacity.
#'
#' @param cell_roost_distances Matrix of distances between cells and roosts.
#' @param N Current roost population sizes.
#' @param N_max Maximum roost capacities.
#' @param lambda Distance-decay parameter.
#' @param alpha Strength of population-based attraction.
#' @param epsilon Small value preventing zero attraction.
#'
#' @return Matrix of normalized roost probabilities for each cell.
#'
#' @export
make_roost_probability <- function(
    cell_roost_distances,
    N,
    N_max,
    lambda = 5,
    alpha = 1,
    epsilon = 0.01
) {
  
  n_roosts <- length(N)
  
  # Population-based attraction
  population_attraction <-
    (
      (N + epsilon) /
        (sum(N) + epsilon * n_roosts)
    )^alpha
  
  # Remaining capacity
  capacity_factor <-
    1 - N / N_max
  
  capacity_factor <- pmax(
    capacity_factor,
    0
  )
  
  attraction <- population_attraction *
    capacity_factor
  
  # Distance-decay kernel
  distance_decay <- exp(
    -lambda * cell_roost_distances
  )
  
  # Combine distance and attraction
  roost_weight <- sweep(
    distance_decay,
    MARGIN = 2,
    STATS = attraction,
    FUN = "*"
  )
  
  # Normalize across roosts for each cell
  row_totals <- rowSums(roost_weight)
  
  roost_probability <- matrix(
    0,
    nrow = nrow(roost_weight),
    ncol = ncol(roost_weight)
  )
  
  valid <- row_totals > 0
  
  roost_probability[valid, ] <-
    roost_weight[valid, , drop = FALSE] /
    row_totals[valid]
  
  roost_probability
}


#' Return bats from the landscape to roosts
#'
#' Converts continuous nighttime spatial densities into discrete roost
#' populations while accounting for distance, population-based attraction,
#' and roost capacity.
#'
#' Individuals are assigned by epidemiological compartment in the order
#' susceptible, infected, and recovered.
#'
#' @param night_output List containing spatial densities `u_S`, `u_I`,
#'   and `u_R`.
#' @param roost_info Data frame containing current roost populations,
#'   locations, and capacities.
#' @param cell_roost_distances Matrix of distances between grid cells
#'   and roosts.
#' @param alpha Strength of population-based attraction.
#' @param epsilon Small value preventing zero population attraction.
#' @param lambda Distance-decay parameter.
#'
#' @return Data frame containing updated roost populations.
#'
#' @export
run_dawn_return <- function(
    night_output,
    roost_info,
    cell_roost_distances,
    alpha = 1,
    epsilon = 0.01,
    lambda = 5
) {
  
  # Convert nighttime densities to cell probabilities
  p_S <- make_cell_probability(
    night_output$u_S
  )
  
  p_I <- make_cell_probability(
    night_output$u_I
  )
  
  p_R <- make_cell_probability(
    night_output$u_R
  )
  
  
  # Roost-return probabilities
  roost_probability <- make_roost_probability(
    cell_roost_distances = cell_roost_distances,
    N = roost_info$N,
    N_max = roost_info$N_max,
    lambda = lambda,
    alpha = alpha,
    epsilon = epsilon
  )
  
  
  # Current populations
  n_S <- sum(roost_info$S)
  n_I <- sum(roost_info$I)
  n_R <- sum(roost_info$R)
  
  
  # Start with empty roosts for the return
  N_current <- rep(
    0,
    nrow(roost_info)
  )
  
  
  # Assign susceptible individuals
  S_new <- assign_bats(
    n_bats = n_S,
    cell_probability = as.vector(p_S),
    roost_probability = roost_probability,
    N_current = N_current,
    N_max = roost_info$N_max
  )
  
  N_current <- S_new
  
  
  # Assign infected individuals
  I_new <- assign_bats(
    n_bats = n_I,
    cell_probability = as.vector(p_I),
    roost_probability = roost_probability,
    N_current = N_current,
    N_max = roost_info$N_max
  )
  
  N_current <- N_current + I_new
  
  
  # Assign recovered individuals
  R_new <- assign_bats(
    n_bats = n_R,
    cell_probability = as.vector(p_R),
    roost_probability = roost_probability,
    N_current = N_current,
    N_max = roost_info$N_max
  )
  
  
  # Construct updated roost populations
  roosts_next <- roost_info
  
  roosts_next$S <- S_new
  roosts_next$I <- I_new
  roosts_next$R <- R_new
  
  roosts_next$N <-
    roosts_next$S +
    roosts_next$I +
    roosts_next$R
  
  
  # Check model invariants
  stopifnot(
    all(roosts_next$N <= roosts_next$N_max),
    sum(roosts_next$S) == n_S,
    sum(roosts_next$I) == n_I,
    sum(roosts_next$R) == n_R,
    all(
      roosts_next$N ==
        roosts_next$S +
        roosts_next$I +
        roosts_next$R
    )
  )
  
  
  list(
    roosts = roosts_next,
    roost_probability = roost_probability
  )
}