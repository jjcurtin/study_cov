# This script fits linear models using our five covariate 
#   selection methods in a simulated data set that varies characteristics of the research 
#   setting (e.g., sample size, number of covariates, effect size for x; see below)

# It will 
# - generate the data
# - fit the models
# - extract parameter estimate, SE, p-value for X 
# - also extract true positive rate and false positive rate for selecting good covariates 
# - save the results in a tibble that has separate rows for each approach and columns 
#   for the above info

# This script will be run across jobs where each job will have a different dataset that varies 
# the research setting/characteristics of the dataset.  These include:
# - n_obs.  number of observations in the dataset
# - b_x.  The effect size for the x effect
# - n_covs.  The number of covariates in the dataset
# - r_ycov.  The correlation between good covariates and y 
# - p_good_covs.  The proportion of covariates that are "good" (nonzero effect)
# - r_cov.  The corrleations among the good covariates ("bad" covariates are 
#      uncorrelated with y, x, good covariates and other bad covariates.)

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
# job_num <- 1
# n_sims <- 10
# n_obs <- 100
# b_x <- 0
# n_covs <- 4
# r_ycov <- .5 
# p_good_covs <- .5
# r_cov <- .3 

source("fun_cov.R")

# Loop over sims
set.seed(job_num)
full_results <- tibble::tibble()

for(i in 1:n_sims) {
  
  di <- generate_data(n_obs, b_x, n_covs, r_ycov, p_good_covs, r_cov)
  
  results <- dplyr::bind_rows(
                       get_results(model = fit_no_covs(di), 
                                   method = "no_covs", 
                                   n_covs, p_good_covs,
                                   sim = i),
                       get_results(model = fit_all_covs(di), 
                                   method = "all_covs", 
                                   n_covs, p_good_covs,
                                   sim = i),
                       get_results(model = fit_p_hacked(di), 
                                   method = "p_hacked", 
                                   n_covs, p_good_covs,
                                   sim = i),
                       get_results(model = fit_r(di), 
                                   method = "r", 
                                   n_covs, p_good_covs,
                                   sim = i),
                       get_results(model = fit_partial_r(di), 
                                   method = "partial_r", 
                                   n_covs, p_good_covs,
                                   sim = i),
                       get_results(model = fit_full_lm(di), 
                                   method = "full_lm", 
                                   n_covs, p_good_covs,
                                   sim = i),
                       get_results(model = fit_lasso(di), 
                                   method = "lasso", 
                                   n_covs, p_good_covs,
                                   sim = i)
                       )
  
  full_results <- dplyr::bind_rows(full_results, results)
}

# add job_num as first column
full_results <- full_results |> 
  dplyr::mutate(job_num = job_num) |> 
  dplyr::relocate(job_num)

# tibble of research setting info
research_setting <- tibble::tibble(job_num = job_num,
                           n_obs = n_obs,
                           b_x = b_x,
                           n_covs = n_covs,
                           r_ycov = r_ycov, 
                           p_good_covs = p_good_covs,
                           r_cov = r_cov)

# join full results and write to csv file
full_results |>
  dplyr::left_join(research_setting, by = "job_num") |> 
  readr::write_csv(stringr::str_c("results_", job_num, ".csv"))