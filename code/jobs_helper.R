
# R script to write jobs csv for CHTC

suppressMessages(library(tidyverse))

n_obs_lvls <- c(100, 200, 400, 1000)
n_covs_lvls <- c(4, 8, 16)
p_good_covs_lvls <- c(0.25, 0.5, 0.75)

jobs <- expand_grid(n_sims = 20000,
                    n_obs = n_obs_lvls,
                    b_x = 0,
                    n_covs = n_covs_lvls,
                    b_cov = 2,
                    p_good_covs = p_good_covs_lvls,
                    rcov = 0) |> 
  mutate(job_num = row_number()) |> 
  relocate(job_num)


jobs |> write_csv("code/jobs.csv")




