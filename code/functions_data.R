
# function to generate data 

# NOTE:  THINK ABOUT HOW BETA INTERACTS WITH ERROR FOR EFFECT SIZES FOR X AND COVS
generate_data <- function(beta_x, beta_covs, n_obs) {

  # constants
  mean_covs <- 0
  sd_covs <- 1
  sd_error <- 1
 
  # make x 
  x <- c(rep(0, n_obs*0.5), rep(1, n_obs*0.5))
 
  # make covs 
  n_covs <- length(beta_covs) 
  covs <- MASS::mvrnorm(n_obs, mu = rep(mean_covs, n_covs), Sigma = diag(n_covs))
  
  # error
  error <- rnorm(n = n_obs, mean = 0, sd = sd_error)
  
  # make y
  y <- beta_x * x + covs %*% beta_covs + error
  y <- drop(y)
 
  # combine all into tibble 
  covs <- covs |>  
    tibble::as_tibble(.name_repair = "universal")
  names(covs) <- paste0("c", 1:n_covs)
 tibble::tibble(x = x, y = y) |> 
   dplyr::bind_cols(covs)
}













