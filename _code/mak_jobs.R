# R script to write jobs csv for CHTC

library(tidyverse)

n_obs_lvls <- c(50, 100, 150, 200, 300, 400)
n_covs_lvls <- c(4, 8, 12, 16, 20)
p_good_covs_lvls <- c(0.25, 0.5, 0.75)
r_ycov_lvls <- c(0.3, 0.5)    # correlation between y & good covs
b_x_lvls <- c(0, 0.2, 0.5)    # effect of x on y

# we want 40,000 simulations
# batching jobs at 500 sims per job
# 500 * 80 repeats = 40,0000 sims
jobs <- expand_grid(n_sims = 500,
                    n_obs = n_obs_lvls,
                    b_x = b_x_lvls,
                    n_covs = n_covs_lvls,
                    r_ycov = r_ycov_lvls,
                    p_good_covs = p_good_covs_lvls,
                    r_cov = 0.3) |> 
  slice(rep(1:n(), each = 80)) |> 
  mutate(job_num = row_number()) |> 
  relocate(job_num)


jobs |> write_csv("_chtc/jobs.csv", col_names = FALSE)

jobs |> select(-job_num) |> distinct() |> nrow()
