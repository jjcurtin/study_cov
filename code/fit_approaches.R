# This is the script that will fit our several approaches to the data to evalute the effect of x
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


d <- generate_data(beta_covs = c(0.5, 0.5), beta_x = 0, n_obs = 100)



# Fit no covariate model




# fit all covariate model


# fit p-hacked model

# fit partial r

# fit LASSO model






# make results tibble