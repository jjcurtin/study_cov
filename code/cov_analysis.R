
# R script for analysis
# this script will:
# concatenate all results csv files saved from CHTC into one dataframe
# calculate what percentage of results has p < 0.05
# group data by research setting
# calculate the type I and type II error rates for when b_x is zero vs nonzero
# visualize the sampling distribution of X effect

setwd("/Volumes/jjcurtin/studydata/cov/raw_data")
getwd()
dir()

suppressMessages(library(dplyr))
library(purrr)
library(ggplot2)


# read in all csv files to a list and row bind them

data <- map(list.files(path = here::here(), pattern = "*.csv", full.names = TRUE), read.csv) |> list_rbind()

# clean responses

data <- data |> mutate(model = gsub("[ -]", "_", model)) |> glimpse()

# statistics

# mean(data$p_value < 0.05)
# data |> filter(model == "lm_p_hacked") |> summarise(prop_sig = mean(p_value < 0.05))

data |> group_by(model) |> summarise(prop_sig = mean(p_value < 0.05))

# sampling distribution 

data |> pull(b_x) |> hist()

data |> pull(b_x) |> density() |> plot(main = "Density Plot of b_x")

ggplot(data, aes(x = b_x)) +
  geom_density(fill = "blue", alpha = 0.7) + 
  labs(title = "Density Plot of b_x", x = "Values", y = "Density") +
  theme_minimal() 



