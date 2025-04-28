
# function to generate data 

# NOTE:  THINK ABOUT HOW BETA INTERACTS WITH ERROR FOR EFFECT SIZES FOR X AND COVS


# generate_data <- function(n_obs, b_x, n_covs, b_cov, p_good_covs, r_cov) {
# 
#   # Currently not using r_cov
#   
#   # constants
#   mean_y <- 0
#   mean_covs <- 0
#   sd_covs <- 1
#   # sd_error <- 1
#  
#   # make x 
#   x <- c(rep(0, n_obs*0.5), rep(1, n_obs*0.5))
#   
#   # make beta covs
#   beta_covs <- c(rep(b_cov, n_covs*p_good_covs), rep(0, n_covs*(1-p_good_covs)))
# 
#   # make sigma
#   sigma = diag(n_covs + 1)
#   
#   # make covs + initial y
#   ycovs <- MASS::mvrnorm(n_obs, mu = c(mean_y, rep(mean_covs, n_covs)), Sigma = sigma)
#   y <- ycovs[,1]
#   # y <- drop(y) # not sure this is needed anymore
#   covs <- ycovs[,-1]
#   
#   
#   # Add x variance into y
#   y <- y + b_x * x 
#  
#   # combine all into tibble and relocate y to first column
#   covs <- covs |>  
#     tibble::as_tibble(.name_repair = "minimal")
#   
#   names(covs) <- paste0("c", 1:n_covs)
#   
#   tibble::tibble(x = x, y = y) |> 
#    dplyr::bind_cols(covs) |> 
#     dplyr::relocate(y)
# }

n_obs <- 100
b_x <- 0.5
n_covs <- 12
r_ycov <- 0.3
p_good_covs <- 0.75
r_cov <- 0

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
  x <- c(rep(0, n_obs*0.5), rep(1, n_obs*0.5))

  # make sigma
  sigma <- diag(n_covs + 1)
  
  # make n good covs
  n_good_covs <- n_covs*p_good_covs
  
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
  
  # make covs + initial y
  ycovs <- MASS::mvrnorm(n_obs, mu = rep(0, n_covs + 1), Sigma = sigma)
  y <- ycovs[,1]
  covs <- ycovs[,-1]
  
  # Add x effect into y
  y <- y + b_x * x 
 
  # combine all into tibble and relocate y to first column
  covs <- covs |>  
    tibble::as_tibble(.name_repair = "minimal")
  
  names(covs) <- paste0("c", 1:n_covs)
  
  tibble::tibble(x = x, y = y) |> 
   dplyr::bind_cols(covs) |> 
    dplyr::relocate(y)
}

# fit no covariate model

get_no_covs_lm <- function(data) {
  lm(y ~ x, data = data)
}


# fit all covariate model

get_all_covs_lm <- function(data) {
  lm(y ~ ., data = data)
}


# fit p-hacked model

get_p_hacked_lm <- function(data) {
  
  lm_base <- lm(y ~ x, data = data) # making base model
  tidy_lm_base <- lm_base |> broom::tidy() |> filter(term == "x")
  p_val_base <- tidy_lm_base$p.value
  
  lm_final <- lm(y ~ x, data = data) # making final model 
  
  n_covs <- grep("^c", names(data), value = TRUE) |> length() # isolate names that start with c
  
  p_val_binary <- rep(0, n_covs) # empty vector of true/false if p-value is smaller
  covs_added <- character(0) # empty vector of covariates added to model
  
  for(i in 1:n_covs) {
    
    ci <- grep("^c", names(data), value = TRUE)[i]
    
    lm_formula <- reformulate(termlabels = c("x", ci), response = "y")
    lm_model <- lm(formula = lm_formula, data = data)
    tidy_lm_model_x <- lm_model |> broom::tidy() |> filter(term == "x")
    
    if(tidy_lm_model_x$p.value < p_val_base) {
      
      p_val_binary[i] <- 1
      
      if(!(ci %in% covs_added)) {
        covs_added <- c(covs_added, ci)
        lm_final <- update(lm_final, paste0(". ~ . + ", ci))
        
      }}
  }
  
  return(lm_final)
}


# fit partial r

get_partial_r_lm <- function(data, alpha = 0.05) {
  
  lm_final <- lm(y ~ x, data = data) # making final model 
  
  n_covs <- grep("^c", names(data), value = TRUE) |> length() # isolate names that start with c
  
  sig_covs <- character(0) # empty vector of covariates that are significant on y
  
  p_vals <- tibble(covariate = character(), 
                   p_value = numeric())
  
  for(i in 1:n_covs) {
    
    ci <- grep("^c", names(data), value = TRUE)[i]
    
    lm_formula <- reformulate(termlabels = c("x", ci), response = "y")
    lm_model <- lm(formula = lm_formula, data = data)
    tidy_lm_model_ci <- lm_model |> broom::tidy() |> filter(term == ci)
    
    p_vals <- add_row(p_vals, covariate = ci, p_value = tidy_lm_model_ci$p.value)
    
    if(tidy_lm_model_ci$p.value < alpha) {
      
      sig_covs <- c(sig_covs, ci)
      lm_final <- update(lm_final, paste0(". ~ . + ", ci))
      
      }
  }
  
  return(lm_final)
}

# fit LASSO model

# packages needed:
# rsample, tidyr, tune, yardstick, parsnip (all in train.sif)

get_lasso_lm <- function(data) {
  
  splits_boot <- data |> rsample::bootstraps(times = 100)
  
  grid_penalty <- tidyr::expand_grid(penalty = exp(seq(-8, 8, length.out = 1000)))
  
  # tune lasso
  
  fits_lasso <-
    parsnip::linear_reg(penalty = tune(), mixture = 1) |> 
    parsnip::set_engine("glmnet",
                        penalty.factor = c(0, rep(1, n_covs))) |> 
    tune::tune_grid(preprocessor = recipes::recipe(y ~ ., data = data),
                    resamples = splits_boot, 
                    grid = grid_penalty, 
                    metrics = yardstick::metric_set(yardstick::rmse))
  
  # select best lasso
  
  fit_best_lasso <-
    parsnip::linear_reg(penalty = tune::select_best(fits_lasso, metric = "rmse")$penalty, 
                        mixture = 1) |> 
    parsnip::set_engine("glmnet",
                        penalty.factor = c(0, rep(1, n_covs))) |> 
    parsnip::fit(y ~ ., data = data)
  
  # select nonzero covariates
  
  nz_covs <- fit_best_lasso |> 
    broom::tidy() |> 
    filter(estimate != 0) |> 
    filter(stringr::str_starts(term, "c")) |> 
    pull(term)
  
  # build formula with nonzero covs and return model
  
  lm_formula <- reformulate(termlabels = c("x", nz_covs), response = "y")
  
  lm(formula = lm_formula, data = data)
  
}


# make results tibble : b, SE, p-value

get_results <- function(lm_model, model_name, sim) {
  
  # make empty results tibble
  
  lm_results <- tibble(model = character(),
                       simulation_id = numeric(),
                       term = character(),
                       estimate = numeric(),
                       SE = numeric(),
                       p_value = numeric())
  
  # tidy output for x
  
  tidy_lm_model_x <- lm_model |> broom::tidy() |> filter(term == "x")
  
  # extract values
  
  lm_results <- bind_rows(lm_results, 
                          tibble(model = model_name, 
                                 simulation_id = sim,
                                 term = "x",
                                 estimate = tidy_lm_model_x$estimate,
                                 SE = tidy_lm_model_x$std.error,
                                 p_value = tidy_lm_model_x$p.value))
  
  return(lm_results)
}






