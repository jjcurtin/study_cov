
# R script to write jobs csv for CHTC

suppressMessages(library(tidyverse))

n_obs_lvls <- c(100, 150, 200, 300, 400)
# n_covs_lvls <- c(4, 8, 12, 16)
p_good_covs_lvls <- c(0.25, 0.5, 0.75)
r_ycov_lvls <- c(0.3, 0.5)    # correlation between y & good covs
# r_cov_lvls <- c(0.3, 0.5)   # correlation among good covs
b_x_lvls <- c(0, 0.3, 0.5)    # effect of x on y

jobs <- expand_grid(n_sims = 500,
                    n_obs = n_obs_lvls,
                    b_x = b_x_lvls,
                    n_covs = 20,
                    r_ycov = r_ycov_lvls,
                    p_good_covs = p_good_covs_lvls,
                    r_cov = 0.3) |> 
  slice(rep(1:n(), each = 40)) |> 
  mutate(job_num = row_number()) |> 
  relocate(job_num)


jobs |> write_csv("chtc/batch_20_covs/jobs.csv", col_names = FALSE)

jobs |> select(-job_num) |> distinct() |> nrow()
