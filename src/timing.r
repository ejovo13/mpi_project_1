# Plotting benchmarks

library(tidyverse)
library(pracma)

benchmark_file <- "csv/timing.csv"

df <- tibble(read.csv(benchmark_file))

# Pivot to compare new and old
labls <- c("least squares", "quadrature")

cp <- c("#000000", "#FF0000")
df <- df |> 
    mutate(least_squares = old, quadrature = new) |>
    select(l, least_squares, quadrature)
    

df.new_old <- df |> pivot_longer(2:3, names_to = "Model", values_to = "time") |> mutate(Model = factor(Model))


levels(df.new_old$Model) <- labls

p <- df.new_old |> 
    mutate(Model = factor(Model)) |>
    ggplot(aes(l, time, col = Model)) + 
    geom_line(aes(linetype = Model)) + 
    geom_point() +
    labs(y = "time (s)", 
         x = "degree",
         title = "Average runtime of models for ETOPO1_small.csv") +
        #  subtitle = "n_runs = 5, n_points = 648000") + 
    scale_color_manual(values = cp) +
    theme(text = element_text(size = 20),
        panel.border = element_blank(),
        panel.grid.minor = element_blank())


## Add to the data frame to get 

ggsave("average_runtime.png", p)