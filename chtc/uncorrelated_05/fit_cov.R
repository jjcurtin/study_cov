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

job_num <- as.numeric(args[1])
n_sims <- as.numeric(args[2])
n_obs <- as.numeric(args[3])
b_x <- as.numeric(args[4])
n_covs <- as.numeric(args[5])
r_ycov <- as.numeric(args[6])
p_good_covs <- as.numeric(args[7])
r_cov <- as.numeric(args[8])

# for testing
# comment out for use on CHTC
# job_num <- 1
# n_sims <- 10
# n_obs <- 100
# b_x <- 0
# n_covs <- 4
# r_ycov <- 2 (previous b_cov)
# p_good_covs <- .5
# r_cov <- 0

suppressMessages(library(dplyr)) 
suppressMessages(library(readr))
suppressMessages(library(stringr))
source("fun_cov.R")


# START LOOP

full_results <- tibble()

set.seed(job_num)

for(i in 1:n_sims) {
  
  di <- generate_data(n_obs, b_x, n_covs, r_ycov, p_good_covs, r_cov)
  
  results <- bind_rows(get_results(lm_model = get_no_covs_lm(data = di), model_name = "lm no covs", sim = i),
                       get_results(lm_model = get_all_covs_lm(data = di), model_name = "lm all covs", sim = i),
                       get_results(lm_model = get_p_hacked_lm(data = di), model_name = "lm p-hacked", sim = i),
                       get_results(lm_model = get_partial_r_lm(data = di), model_name = "lm partial r", sim = i),
                       get_results(lm_model = get_lasso_lm(data = di), model_name = "lm lasso", sim = i))
  
  full_results <- bind_rows(full_results, results)
  
}

# add job_num as first column
full_results <- full_results |> 
  mutate(job_num = job_num) |> 
  relocate(job_num)

# tibble of research setting info
research_setting <- tibble(job_num = job_num,
                           n_obs = n_obs,
                           b_x = b_x,
                           n_covs = n_covs,
                           r_ycov = r_ycov, 
                           p_good_covs = p_good_covs,
                           r_cov = r_cov)

# join full results and write to csv file
full_results |>
  left_join(research_setting, by = "job_num") |> 
  write_csv(str_c("results_", job_num, ".csv"))
