################################################################################
# General functions (across methods)
################################################################################

# function to generate data 

generate_data <- function(n_obs, b_x, n_covs, r_ycov, p_good_covs, r_cov) {

  # generates y and covs with mean=0, variance=1 
  # generates dichotomous x coded -.5, .5
  # n_obs = sample size
  # b_x = x effect
  # n_covs = number of covariates
  # r_ycov = correlation between y and covariates
  # p_good_covs = proportion of "good" covariates (nonzero effect)
  # r_cov = correlation among "good" covariates
 
  # make x 
  x <- c(rep(0, n_obs/2), rep(1, n_obs/2))

  # make sigma for covs and y
  sigma <- diag(n_covs + 1)
  
  # make n good covs
  n_good_covs <- n_covs * p_good_covs
  
  # make correlation matrix of good predictors
  corr_matrix <- matrix(r_cov, 
                        nrow = n_good_covs, 
                        ncol = n_good_covs)
  diag(corr_matrix) <- 1
 
  # superimpose corr_matrix onto sigma 
  sigma[2:(nrow(corr_matrix) + 1), 2:(ncol(corr_matrix) + 1)] <- corr_matrix
  
  # add correlations between y and good covs
  sigma[1, 2:(n_good_covs + 1)] <- rep(r_ycov, n_good_covs)
  sigma[2:(n_good_covs + 1), 1] <- rep(r_ycov, n_good_covs)
  
  # make covs + y pre-manipulation of x
  # all with mean = 0 and variance = 1
  ycovs <- MASS::mvrnorm(n_obs, mu = rep(0, n_covs + 1), Sigma = sigma)
  y <- ycovs[, 1]
  covs <- ycovs[, -1]
  
  # Add x effect into y
  y <- y + b_x * x 
 
  # combine all into tibble
  covs <- covs |>  
    tibble::as_tibble(.name_repair = "minimal")
  
  names(covs) <- stringr::str_c("c", 1:n_covs)
  
  tibble::tibble(y = y, x = x) |> 
   dplyr::bind_cols(covs)
}


# LAUREN: Add in proportion of good covariates found (hits) and proportion 
#   of bad covariates included (false alarms)
# LAUREN: I renamed model to method both as the input to the function and in the 
#   tibble that is returned to be clearer on what that string describes
#   I fixed this in fit_cov.R but it likely will also affect affect your qmd.  
#   You will need to update it there.

# make results tibble : b, SE, p-value
get_results <- function(model, method, sim) {
  # model: an lm object
  # method: The name of the method used to select covariates (e.g., ) 
  # sim: simulation number

  ndf <- model |> broom::glance() |> dplyr::pull(df)
  ddf <- model |> broom::glance() |> dplyr::pull(df.residual)
  output <- model |> broom::tidy() |> dplyr::filter(term == "x")
  
  tibble::tibble(method = method, 
         simulation_id = sim,
         estimate = output$estimate,
         SE = output$std.error,
         p_value = output$p.value,
         ndf = ndf,
         ddf = ddf) 
}


###############################################################################
# Fit functions
###############################################################################

# fit no covariate model
fit_no_covs <- function(d) {
  lm(y ~ x, data = d)
}


# fit all covariate model
fit_all_covs <- function(d) {
  lm(y ~ ., data = d)
}


# fit p-hacked model
fit_p_hacked <- function(d, n_covs) {
  
  # making base model
  lm_base <- lm(y ~ x, data = d) 
  p_base <- lm_base |> 
    broom::tidy() |> 
    dplyr::filter(term == "x") |> 
    dplyr::pull(p.value)
  
  # empty vector of covariates added to model
  covs_added <- character(0) 
  
  for(i in 1:n_covs) {
    
    ci <- grep("^c", names(d), value = TRUE)[i]
    
    formula_1cov <- reformulate(termlabels = c("x", ci), response = "y")
    lm_1cov <- lm(formula = formula_1cov, data = d)
    p_1cov <- lm_1cov |> 
      broom::tidy() |> 
      dplyr::filter(term == "x") |> 
      dplyr::pull(p.value)
    
    if(p_1cov < p_base) {
      covs_added <- c(covs_added, ci)
    }
  }

  # fit model with hacked covariates 
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  lm(formula = formula_final, data = d)
}


# fit partial r
fit_partial_r <- function(d, n_covs, alpha = 0.05) {
  
  # empty vector of covariates that are significant on y, controlling x
  covs_added <- character(0) 
  
  for(i in 1:n_covs) {
    ci <- grep("^c", names(d), value = TRUE)[i]
    
    formula_1cov <- reformulate(termlabels = c("x", ci), response = "y")
    lm_1cov <- lm(formula = formula_1cov, data = d)
    p_1cov <- lm_1cov |> 
      broom::tidy() |> 
      dplyr::filter(term == ci) |> 
      dplyr::pull(p.value)
    
    if(p_1cov < alpha) {
      covs_added <- c(covs_added, ci)
    }
  }
  
  # fit model with partial r covariates 
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  lm(formula = formula_final, data = d)
}

# fit LASSO model
fit_lasso <- function(d, n_covs) {
  
  splits_boot <- d |> rsample::bootstraps(times = 100)
  grid_penalty <- tidyr::expand_grid(penalty = exp(seq(-8, 8, length.out = 1000)))
  
  # tune lasso
  fits_lasso <-
    parsnip::linear_reg(penalty = tune(), mixture = 1) |> 
    parsnip::set_engine("glmnet",
                        penalty.factor = c(0, rep(1, n_covs))) |> 
    tune::tune_grid(preprocessor = recipes::recipe(y ~ ., data = d),
                    resamples = splits_boot, 
                    grid = grid_penalty, 
                    metrics = yardstick::metric_set(yardstick::rmse))
  
  # select best lasso
  fit_best_lasso <-
    parsnip::linear_reg(penalty = tune::select_best(fits_lasso, metric = "rmse")$penalty, 
                        mixture = 1) |> 
    parsnip::set_engine("glmnet",
                        penalty.factor = c(0, rep(1, n_covs))) |> 
    parsnip::fit(y ~ ., data = d)
  
  # select nonzero covariates
  covs_added <- fit_best_lasso |> 
    broom::tidy() |> 
    dplyr::filter(estimate != 0) |> 
    dplyr::filter(stringr::str_starts(term, "c")) |> 
    dplyr::pull(term)
  
  # build formula with nonzero covs and return model
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  
  lm(formula = formula_final, data = d)
}