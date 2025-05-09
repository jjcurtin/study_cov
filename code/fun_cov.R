################################################################################
# Notes
################################################################################
# Script contains functions for generating data and fitting models

# Requires the following packages:
#   dplyr, broom, stringr, rsample, tidyr, tune, yardstick, parsnip, recipes,
#   tibble, MASS, glmnet

# Functions included:
#   generate_data(n_obs, b_x, n_covs, r_ycov, p_good_covs, r_cov)
#   get_results(model, method, n_covs, p_good_covs, sim)
#   fit_no_covs(d)
#   fit_all_covs(d)
#   fit_p_hacked(d)
#   fit_partial_r(d)
#   fit_lasso(d)

################################################################################
# General functions (across covariate selection methods)
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


# make results tibble
get_results <- function(model, method, n_covs, p_good_covs, sim) {
  # model: an lm object
  # method: The name of the method used to select covariates (e.g., ) 
  # sim: simulation number
 
  ndf <- model |> broom::glance() |> dplyr::pull(df)
  ddf <- model |> broom::glance() |> dplyr::pull(df.residual)
  output <- model |> broom::tidy() |> dplyr::filter(term == "x")
 
  # get proportion of good covariates found (true positive rate for covs) and 
  # proportion of bad covariates included (false positive rate for covs)
  good_covs <- stringr::str_c("c", 1:(n_covs * p_good_covs))
  covs_included <- broom::tidy(model) |> 
    dplyr::pull(term) |> 
    stringr::str_subset("^c")
  covs_tpr <- sum(covs_included %in% good_covs) / (n_covs * p_good_covs)
  covs_fpr <- sum(!(covs_included %in% good_covs)) / (n_covs - (n_covs * p_good_covs))
 
  # put it all in a results tibble 
  tibble::tibble(method = method, 
         simulation_id = sim,
         estimate = output$estimate,
         SE = output$std.error,
         p_value = output$p.value,
         ndf = ndf,
         ddf = ddf,
         covs_tpr = covs_tpr,
         covs_fpr = covs_fpr) 
}


###############################################################################
# Function fit lm using various covariate selection methods 
###############################################################################

# fit no covariate model
# Fits linear model with no covariates
fit_no_covs <- function(d) {
  lm(y ~ x, data = d)
}


# fit all covariate model
# Fits linear model with all available covariates
fit_all_covs <- function(d) {
  lm(y ~ ., data = d)
}


# fit p-hacked model
# Fits linear model with covariates that improve p-value for x.  Considers each 
#   covariate one at a time and includes it if it reduces the p-value for the x
#   effect from the p-value for x from the simple x only model
fit_p_hacked <- function(d) {
  
  # making base model
  lm_base <- lm(y ~ x, data = d) 
  p_base <- lm_base |> 
    broom::tidy() |> 
    dplyr::filter(term == "x") |> 
    dplyr::pull(p.value)
  
  # empty vector of covariates added to model
  covs_added <- character(0) 
 
  # calculate n_covs from data 
  n_covs <- names(d) |> 
    stringr::str_subset("^c") |> 
    length()
  
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

  # fit model with selected covariates 
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  lm(formula = formula_final, data = d)
}

# fit r
# Fits linear model with covariates that are significant on y (overall).
#   Considers each covariate one at at a time and includes it if it is significant
#   in a linear model (i.e., if bivarte r is significant).
fit_r <- function(d, alpha = 0.05) {
  
  # empty vector of covariates that are significant on y, controlling x
  covs_added <- character(0) 
  
  # calculate n_covs from data 
  n_covs <- names(d) |> 
    stringr::str_subset("^c") |> 
    length()
  
  for(i in 1:n_covs) {
    ci <- grep("^c", names(d), value = TRUE)[i]
    formula_1cov <- reformulate(termlabels = c(ci), response = "y")
    lm_1cov <- lm(formula = formula_1cov, data = d)
    p_1cov <- lm_1cov |> 
      broom::tidy() |> 
      dplyr::filter(term == ci) |> 
      dplyr::pull(p.value)
    
    if(p_1cov < alpha) {
      covs_added <- c(covs_added, ci)
    }
  }
  
  # fit model with selected covariates 
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  lm(formula = formula_final, data = d)
}

# fit partial r
# Fits linear model with covariates that are significant on y, controlling for x.
#   Considers each covariate one at at a time and includes it if it is significant
#   in a linear model that also includes x.
fit_partial_r <- function(d, alpha = 0.05) {
  
  # empty vector of covariates that are significant on y, controlling x
  covs_added <- character(0) 
  
  # calculate n_covs from data 
  n_covs <- names(d) |> 
    stringr::str_subset("^c") |> 
    length()
  
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
  
  # fit model with selected covariates 
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  lm(formula = formula_final, data = d)
}


# fit full_lm 
# Fits linear model with covariates that are significant on y, controlling for x 
# and all other covariates.
fit_full_lm <- function(d, alpha = 0.05) {
 
  # get vector of sig covs from full model with x and all covs 
  covs_added <- lm(y ~ ., data = d) |>  
    broom::tidy() |> 
    dplyr::filter(stringr::str_starts(term, "c")) |> 
    dplyr::filter(p.value < .05) |> print() |> 
    dplyr::pull(term)
  
  # fit model with selected covariates 
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  lm(formula = formula_final, data = d)
}

# fit LASSO model
# Fits linear model with covariates selected by LASSO.  Tunes (across 100 bootstraps) 
#   and fits a best LASSO model using x and all covariates with no penalty applied 
#   to x (so it is never dropped).  Selects those covariates that had non-zero 
#   coefficients in the best LASSO model and includes them in final linear model 
#   with x.
fit_lasso <- function(d) {
  
  # calculate n_covs from data 
  n_covs <- names(d) |> 
    stringr::str_subset("^c") |> 
    length()

  splits_boot <- d |> rsample::bootstraps(times = 100)
  
  # use a very wide set of penalties/lambda
  grid_penalty <- tidyr::expand_grid(penalty = exp(seq(-8, 8, length.out = 1000)))
  
  # tune lasso
  # no need to standardize covariates because all have same variance and penalty
  # will not be applied to x
  fits_lasso <-
    parsnip::linear_reg(penalty = tune(), mixture = 1) |> # lasso model
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
  
  # fit model with selected covariates 
  formula_final <- reformulate(termlabels = c("x", covs_added), response = "y")
  
  lm(formula = formula_final, data = d)
}