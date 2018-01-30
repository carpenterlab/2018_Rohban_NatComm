## borrowed from Claire's code in issues/115 cytominer. Will use her code once the pull request is confirmed. 
library(Matrix)

generate_component_matrix <- function(n_features, n_components, density) {
  # Generate nonzero elements - follows the binomial distribution:
  #
  #   1 (nonzero) with probability density
  #   0           with probability 1 - density
  #
  # (Requires n_features * n_components space)
  nonzero_elements <- rbinom(n_features * n_components, 1, density)
  indices <- which(nonzero_elements != 0)
  
  # Generate sign of the nonzero elements - the probability of
  # positive or negative sign is equal (0.5).
  signs <- 2 * rbinom(length(indices), 1, 0.5) - 1
  
  # Compute the value of the nonzero elements
  k <- sqrt(1.0 / (density * n_components))
  
  # Generate the component matrix
  component_matrix <- sparseMatrix(
    i = (indices - 1) %% n_features + 1,
    j = (indices - 1) %/% n_features + 1,
    x = k * signs,
    dims = c(n_features, n_components)
  )
  
  component_matrix
}