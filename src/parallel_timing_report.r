library(pracma)
library(tidyverse)

csv_path <- "csv/parallel_timing.csv"

df <- as_tibble(read.csv(csv_path))

df.speedup <- df |> 
    mutate(t25_speedup = t25[1] / t25) |> 
    mutate(t50_speedup = t50[1] / t50) |>
    mutate(t100_speedup = t100[1] / t100) |>
    select(n_threads, t25_speedup, t50_speedup, t100_speedup)