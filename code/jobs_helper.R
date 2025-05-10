
# R script to write jobs csv for CHTC

suppressMessages(library(tidyverse))

n_obs_lvls <- c(100, 150, 200, 300, 400)
n_covs_lvls <- c(4, 8, 12, 16)
p_good_covs_lvls <- c(0.25, 0.5, 0.75)
r_ycov_lvls <- c(0.3, 0.5)  # correlation between y and covs
r_cov_lvls <- 0.3           # correlation among good covs

jobs <- expand_grid(n_sims = 500,
                    n_obs = n_obs_lvls,
                    b_x = 0.5,
                    n_covs = n_covs_lvls,
                    r_ycov = r_ycov_lvls,
                    p_good_covs = p_good_covs_lvls,
                    r_cov = r_cov_lvls) |> 
  slice(rep(1:n(), each = 40)) |> 
  mutate(job_num = row_number()) |> 
  relocate(job_num)


jobs |> write_csv("chtc/new_methods_05/jobs.csv", col_names = FALSE)

jobs |> select(-job_num) |> distinct() |> nrow()
