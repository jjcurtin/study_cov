
# R script for analysis
# this script will:
# concatenate all results csv files saved from CHTC into one dataframe
# calculate what percentage of results has p < 0.05
# group data by research setting
# calculate the type I and type II error rates for when b_x is zero vs nonzero
# visualize the sampling distribution of X effect


suppressMessages(devtools::source_url("https://github.com/jjcurtin/lab_support/blob/main/format_path.R?raw=true"))
(rdrive_path <- format_path("studydata/cov/raw_data"))



# load packages + set ggplot theme

suppressMessages(library(dplyr))
library(skimr)
library(purrr)
library(ggplot2)

theme_set(theme_classic())

# read in jobs csv

jobs <- read.csv("jobs.csv", header = FALSE,
                 col.names = c("job_num", "n_sims", "n_obs", "b_x", "n_covs", "b_cov", "p_good_covs", "r_cov"))

# cd to research drive raw_data folder

setwd("/Volumes/jjcurtin/studydata/cov/raw_data")
getwd()
dir()

# read in all (results) csv files to a list and row bind them

data <- map(list.files(path = here::here(), pattern = "*.csv", full.names = TRUE), read.csv) |> list_rbind()

# clean responses

data <- data |> mutate(model = gsub("[ -]", "_", model)) |> glimpse()

# proportion significant

data |> group_by(model) |> summarise(prop_sig = mean(p_value < 0.05))

# group by research setting

data |> group_by(job_num, model) |> summarise(prop_sig = mean(p_value < 0.05))

# sampling distribution beta_x estimate

plot_hist_density <- function(data, model_name) {
  ggplot(data |> filter(model == model_name), aes(x = estimate)) +
    geom_histogram(aes(y = after_stat(density)), bins = 10, fill = "grey", color = "black") +
    geom_density(fill = "blue", alpha = 0.25, lwd = 1) + 
    labs(title = paste("Density Plot:", model_name), x = "Beta_X Estimate", y = "Density")
}

plot_hist_density(data, "lm_no_covs")
plot_hist_density(data, "lm_all_covs")
plot_hist_density(data, "lm_p_hacked")
plot_hist_density(data, "lm_partial_r")




