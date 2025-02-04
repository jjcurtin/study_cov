
# function to generate data 


generate_data <- function(beta_covs, beta_x, n_obs, n_covs, mean_covs = 0, sd_covs = 1) {
  
  covs <- MASS::mvrnorm(n_obs, mu = rep(mean_covs, n_covs), Sigma = diag(n_covs))
  
  x <- c(rep(0, n_obs*0.5), rep(1, n_obs*0.5)) |> sample()
  
  error <- rnorm(n = n_obs, mean = 0, sd = 1)
  
  y <- beta_x * x + covs %*% beta_covs + error
  
  data <- tibble::tibble(
    !!!purrr::set_names(tibble::as_tibble(covs), paste0("c", 1:n_covs)),
    x = x,
    y = drop(y))
  
  return(data)
  
}














