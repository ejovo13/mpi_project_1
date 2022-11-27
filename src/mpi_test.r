library(tidyverse)

file <- "csv/hercule_mpi_timing.csv"

n_thread <- c(1, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200)

df <- as_tibble(read_csv(file))

df.mean <- df |> mutate(across(everything(), mean))

means <- unlist(df.mean[1,-1])

df.mean <- tibble(n = n_thread, mu = means)

df.mean <- df.mean |> mutate(speedup = mu[1] / mu) |> mutate(scalability = speedup / n) |> filter(n != 1)

# Now let's get the scalability
p <- df.mean |> 
    ggplot(aes(n, speedup)) + 
    geom_point() + 
    # geom_line(col = "blue") + 
    # geom_smooth(col = "blue") + 
    geom_smooth() + 
    geom_line(aes(x = n_thread[-1], y = n_thread[-1]), size = 0.9) +
    theme(text = element_text(size = 20)) +
    labs(x = "Number MPI threads", y = "Speedup")

ggsave("speedup_mpi_200.png", p)

p <- df.mean |> 
    ggplot(aes(n, scalability)) + 
    geom_smooth() + 
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    geom_line(aes(x = n_thread[-1], y = 1), size = 0.9) +
    theme(text = element_text(size = 20)) +
    labs(x = "Number MPI threads", y = "Scalability")

ggsave("scalability_mpi_200.png", p)