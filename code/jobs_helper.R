
# R script to write jobs csv for CHTC

suppressMessages(library(tidyverse))

n_obs_lvls <- c(100, 200, 400, 1000)
n_covs_lvls <- c(4, 8, 12)
p_good_covs_lvls <- c(0.25, 0.5, 0.75)

jobs <- expand_grid(n_sims = 500,
                    n_obs = n_obs_lvls,
                    b_x = 0,
                    n_covs = n_covs_lvls,
                    r_ycov = 0.3, # previous b_cov
                    p_good_covs = p_good_covs_lvls,
                    r_cov = 0.3) |> 
  slice(rep(1:n(), each = 40)) |> 
  mutate(job_num = row_number()) |> 
  relocate(job_num)


jobs |> write_csv("code/jobs_corr.csv", col_names = FALSE)


jobs |> select(-job_num) |> distinct() |> nrow()


