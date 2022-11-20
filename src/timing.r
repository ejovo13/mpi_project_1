# Plotting benchmarks

library(tidyverse)
library(pracma)

benchmark_file <- "csv/timing.csv"

df <- tibble(read.csv(benchmark_file))

# Pivot to compare new and old

df <- df |> 
    mutate(least_squares = old, quadrature = new) |>
    select(l, least_squares, quadrature)
    

df.new_old <- df |> pivot_longer(2:3, names_to = "model", values_to = "time")


p <- df.new_old |> 
    mutate(model = factor(model)) |>
    ggplot(aes(l, time, col = model)) + 
    geom_line() + 
    geom_point() +
    labs(y = "time (s)", 
         x = "l_max",
         title = "Average runtime of models for ETOPO1_small.csv",
         subtitle = "n_runs = 5, n_points = 648000")

ggsave("average_runtime.png", p)