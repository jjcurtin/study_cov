# This is the script that will fit our several approaches to the data to evaluate the effect of x
# It will 
# 1. generate the data
# fit the models
# extract parameter estimate, SE, p-value for X
# extract other info??? (e.g., about covariate effects maybe)
# save the results in a tibble that has separate rows for each approach and columns for the above info

# This script will be run across jobs where each job will have a different dataset that varies on
# the population effect size for x, the sample size, the number of covariates, the size of the covariate effects etc.

args <- commandArgs(trailingOnly = TRUE) 

beta_x <- as.numeric(args[1])
n_obs <- as.numeric(args[2])
n_covs <- as.numeric(args[3])
beta_cov <- as.numeric(args[4])
p_sigcovs <- as.numeric(args[5])

# for testing
# comment out for use on CHTC
beta_x <- .5
n_obs <- 100
n_covs <- 2
p_sigcovs <- .5
setwd("./code") 


library(tidyverse)
source("functions_data.R")

# make beta_covs


# generate data

d <- generate_data(beta_covs = c(0.5, 0.5), beta_x = 0, n_obs = 100, n_covs = 2)


# fit no covariate model

get_no_covs_lm <- function(data) {
  
  lm(y ~ x, data = data)
  
}


# fit all covariate model

get_all_covs_lm <- function(data) {
  
  lm_final <- lm(y ~ x, data = data)
  
  n_covs <- grep("^c", names(data), value = TRUE) |> length() # isolate names that start with c
  
  for(i in 1:n_covs) {
    
    ci <- grep("^c", names(data), value = TRUE)[i]
    
    lm_final <- update(lm_final, paste0(". ~ . + ", ci))
  
  }
  return(lm_final)
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

# -- hold --

# make results tibble
# b, SE, p-value

get_results <- function(lm_model, model_name) {
  
  # make empty results tibble
  
  lm_results <- tibble(model = character(),
                       term = character(),
                       estimate = numeric(),
                       SE = numeric(),
                       p_value = numeric())
  
  # tidy output for x
  
  tidy_lm_model_x <- lm_model |> broom::tidy() |> filter(term == "x")
  
  # extract values
  
  lm_results <- bind_rows(lm_results, 
                          tibble(model = model_name, 
                                 term = "x",
                                 estimate = tidy_lm_model_x$estimate,
                                 SE = tidy_lm_model_x$std.error,
                                 p_value = tidy_lm_model_x$p.value))
  
  return(lm_results)
}

(lm_no_covs <- get_no_covs_lm(data = d))
(lm_all_covs <- get_all_covs_lm(data = d))
(lm_p_hack <- get_p_hacked_lm(data = d))
(lm_partial_r <- get_partial_r_lm(data = d, alpha = 0.05))

get_results(lm_model = lm_no_covs, model_name = "lm no covs")
get_results(lm_model = lm_all_covs, model_name = "lm all covs")
get_results(lm_model = lm_p_hack, model_name = "lm p-hacked")
get_results(lm_model = lm_partial_r, model_name = "lm partial r")

(results_tibble <- bind_rows(get_results(lm_model = lm_no_covs, model_name = "lm no covs"),
                             get_results(lm_model = lm_all_covs, model_name = "lm all covs"),
                             get_results(lm_model = lm_p_hack, model_name = "lm p-hacked"),
                             get_results(lm_model = lm_partial_r, model_name = "lm partial r")))
