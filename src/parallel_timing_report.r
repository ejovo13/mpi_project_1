library(pracma)
library(tidyverse)

csv_path <- "csv/parallel_timing.csv"

df <- as_tibble(read.csv(csv_path))

legend_order <- c("150", "100", "50", "25")

df.speedup <- df |> 
    mutate(t25_speedup = t25[1] / t25) |> 
    mutate(t50_speedup = t50[1] / t50) |>
    mutate(t100_speedup = t100[1] / t100) |>
    mutate(t150_speedup = t150[1] / t150) |>
    select(n_threads, t25_speedup, t50_speedup, t100_speedup, t150_speedup)

df.sp.longer <- df.speedup |>
    pivot_longer(c(2:5), names_to = "model_order", values_to = "speedup")

df.sp.longer |> 
    ggplot(aes(n_threads, speedup)) + 
    geom_line(aes(col = model_order)) + 
    geom_point() +
    scale_x_continuous(breaks = 1:12) +
    scale_colour_discrete(labels = legend_order, name = "Model order") +
    labs(x = "n threads", y = "speedup factor")

# Plot the functions of time

df.longer <- df |>
    pivot_longer(c(2:5), names_to = "model_order", values_to = "time") |>
    mutate(model_order = factor(model_order, levels = c("t150", "t100", "t50", "t25")))



# df.longer |> ggplot(aes(n_threads, time)) + geom_line() + facet_wrap(~ model_order)

df.longer |> 
    ggplot(aes(n_threads, time)) + 
    geom_line(aes(col = model_order)) + 
    geom_point() +
    scale_x_continuous(breaks = 1:12) +
    scale_colour_discrete(labels = legend_order, name = "Model order") +
    labs(x = "n threads", y = "time (s)")

