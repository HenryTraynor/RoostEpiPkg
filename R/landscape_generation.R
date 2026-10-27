
#' Generate Spatial Landscapes Across Modernization and Spatial Gradients
#'
#' Generates spatial landscapes with balanced coverage across intervals
#' of the Clark-Evans aggregation index and creates a crossed experimental
#' design across proportions of modern roosts.
#'
#' @param n_per_bin Integer. Number of landscapes to retain in each
#'   Clark-Evans R interval.
#' @param R_breaks Numeric vector of strictly increasing interval
#'   boundaries for the Clark-Evans R statistic.
#' @param prop_modern_values Numeric vector of proportions of modern
#'   roosts to simulate. Values must be between 0 and 1.
#' @param num_roosts Integer. Number of roosts in each landscape.
#' @param cluster_values Numeric vector of candidate numbers of clusters
#'   used to generate spatial arrangements.
#' @param sd_values Numeric vector of candidate within-cluster standard
#'   deviations used to generate spatial arrangements.
#' @param modern_N_max Numeric scalar. Mean maximum occupancy capacity
#'   for modern roosts.
#' @param traditional_N_max Numeric scalar. Mean maximum occupancy
#'   capacity for traditional roosts.
#' @param min_N_max Numeric scalar. Minimum allowed roost capacity.
#' @param max_attempts Integer. Maximum number of candidate landscapes
#'   to generate while filling all spatial intervals.
#'
#' @return A list with four components:
#' \describe{
#'   \item{landscapes}{List of selected base landscapes. Each element
#'   contains a landscape identifier, the base roost data frame,
#'   the Clark-Evans R value, the generation parameters,
#'   and the spatial interval index. The base roost data include the
#'   fixed modernization order and potential capacities for both
#'   building types.}
#'   \item{design}{Data frame containing one row for each combination
#'   of landscape and modernization level. The \code{roost_info} column
#'   is a list-column containing the corresponding roost data frame.}
#'   \item{counts}{Integer vector giving the number of accepted
#'   landscapes in each Clark-Evans R interval.}
#'   \item{R_breaks}{The interval boundaries supplied to the function.}
#' }
#'
#' @examples
#' \dontrun{
#' R_breaks <- seq(0.1, 1.3, length.out = 10)
#'
#' landscape_set <- generate_landscapes(
#'   n_per_bin = 10,
#'   R_breaks = R_breaks,
#'   prop_modern_values = seq(0, 1, by = 0.2),
#'   num_roosts = 100,
#'   cluster_values = cluster_seq,
#'   sd_values = sd_seq,
#'   modern_N_max = 100,
#'   traditional_N_max = 20,
#'   min_N_max = 5
#' )
#'
#' simulation_grid <- landscape_set$design
#' }
#'
#' @export

generate_landscapes <- function(
    n_per_bin,
    R_breaks,
    prop_modern_values = seq(0, 1, by = 0.2),
    num_roosts = 100,
    cluster_values,
    sd_values,
    modern_N_max = 100,
    traditional_N_max = 20,
    min_N_max = 5,
    max_attempts = 10000
) {

  n_bins <- length(R_breaks) - 1

  # Validate inputs
  stopifnot(
    n_per_bin > 0,
    length(R_breaks) >= 2,
    all(diff(R_breaks) > 0),
    all(prop_modern_values >= 0 & prop_modern_values <= 1)
  )

  accepted <- vector("list", n_bins)
  counts <- integer(n_bins)

  for (b in seq_len(n_bins)) {
    accepted[[b]] <- list()
  }

  # Generate landscapes until each R bin has enough spatial arrangements
  for (attempt in seq_len(max_attempts)) {

    if (all(counts >= n_per_bin)) {
      break
    }

    nc <- sample(cluster_values, 1)
    cluster_sd <- sample(sd_values, 1)

    # Generate coordinates and base roost information
    roost_info <- create_roosts(
      num_roosts = num_roosts,
      num_clusters = nc,
      sd = cluster_sd,
      prop_modern = 0
    )

    R_value <- calculate_clark_evans_R(
      roost_info$x,
      roost_info$y
    )

    if (!is.finite(R_value)) next

    bin <- findInterval(
      R_value,
      R_breaks,
      rightmost.closed = TRUE
    )

    if (bin < 1 || bin > n_bins) next
    if (counts[bin] >= n_per_bin) next

    # Assign a fixed random order for modernization
    roost_info$modernization_order <- sample.int(
      nrow(roost_info)
    )

    # Generate fixed capacity values for each roost
    roost_info$N_max_traditional <- pmax(
      round(rnorm(
        num_roosts,
        mean = traditional_N_max,
        sd = traditional_N_max * 0.25
      )),
      min_N_max
    )

    roost_info$N_max_modern <- pmax(
      round(rnorm(
        num_roosts,
        mean = modern_N_max,
        sd = modern_N_max * 0.25
      )),
      min_N_max
    )

    counts[bin] <- counts[bin] + 1L

    accepted[[bin]][[counts[bin]]] <- list(
      landscape_id = paste0(
        "B", bin, "_L", counts[bin]
      ),
      base_roosts = roost_info,
      Clark_Evans_R = R_value,
      num_clusters = nc,
      sd = cluster_sd,
      R_bin = bin
    )
  }

  if (any(counts < n_per_bin)) {
    stop(
      "Could not fill every R bin. Counts: ",
      paste(counts, collapse = ", "),
      ". Increase max_attempts or adjust R_breaks."
    )
  }

  # Flatten the selected landscapes
  landscapes <- unlist(
    accepted,
    recursive = FALSE
  )

  # Build cross
  design <- expand.grid(
    landscape_index = seq_along(landscapes),
    prop_modern = prop_modern_values,
    KEEP.OUT.ATTRS = FALSE
  )

  design$landscape_id <- vapply(
    landscapes[design$landscape_index],
    function(x) x$landscape_id,
    character(1)
  )

  design$Clark_Evans_R <- vapply(
    landscapes[design$landscape_index],
    function(x) x$Clark_Evans_R,
    numeric(1)
  )

  design$num_clusters <- vapply(
    landscapes[design$landscape_index],
    function(x) x$num_clusters,
    numeric(1)
  )

  design$sd <- vapply(
    landscapes[design$landscape_index],
    function(x) x$sd,
    numeric(1)
  )

  design$R_bin <- vapply(
    landscapes[design$landscape_index],
    function(x) x$R_bin,
    integer(1)
  )

  design$replicate <- match(
    design$landscape_index,
    unique(design$landscape_index)
  )

  # Apply modernization while preserving coordinates
  design$roost_info <- lapply(
    seq_len(nrow(design)),
    function(i) {

      roost_info <- landscapes[[
        design$landscape_index[i]
      ]]$base_roosts

      n <- nrow(roost_info)

      n_modern <- round(
        design$prop_modern[i] * n
      )

      roost_info$roost_type <- ifelse(
        roost_info$modernization_order <= n_modern,
        "modern",
        "traditional"
      )

      roost_info$N_max <- ifelse(
        roost_info$roost_type == "modern",
        roost_info$N_max_modern,
        roost_info$N_max_traditional
      )

      # Reset populations; initialize separately before
      # running each epidemic simulation.
      roost_info$S <- 0L
      roost_info$I <- 0L
      roost_info$R <- 0L
      roost_info$N <- 0L

      roost_info
    }
  )

  list(
    landscapes = landscapes,
    design = design,
    counts = counts,
    R_breaks = R_breaks
  )
}
