# Plotting scalability 
library(tidyverse)

n_runs <- 10
degrees <- c(50, 100, 200)

# convert timing values to scala by multiplying by (1 / t0)
f <- function(c) {
    c[1] / c
}

load_df <- function(degree, nruns) {
    filename <- paste("csv/mpi/mpi_", degree, "_", n_runs, ".csv", sep = "")

    df <- as_tibble(read_csv(filename))
    
    df.mean <- df |> mutate(across(everything(), .fn = mean))

    x_mean <- unlist(df.mean[1,2:13])
    np <- 1:12

    tibble(mu = x_mean, p = np, k = x_mean[1] / x_mean)
}

df50 <- load_df(50, n_runs)
df100 <- load_df(100, n_runs)
df200 <- load_df(200, n_runs)

df <- df50 |> 
    mutate(mu50 = mu, mu100 = df100$mu, mu200 = df200$mu) |>
    mutate(k50 = mu50[1] / mu50, k100 = mu100[1] / mu100, k200 = mu200[1] / mu200) |>
    select(c(p, k50, k100, k200, mu50, mu100, mu200))
    
df.mu_long <- df |>
    select(c(p, mu50, mu100, mu200)) |>
    pivot_longer(c(2:4), names_to = "model", values_to = "mu")

df.k_long <- df |>
    select(c(p, k50, k100, k200)) |>
    pivot_longer(c(2:4), names_to = "model", values_to = "k")

replace <- function(mu_model) {

    fn <- function(mu) {

        if (mu == "mu50") { return ("50") }
        if (mu == "mu100") { return ("100") }
        if (mu == "mu200") { return ("200") }
    }

    purrr::map_chr(mu_model, fn)
}

lvl <- c("50", "100", "200")

df.long <- df.mu_long |> 
    mutate(
        k = df.k_long$k, 
        model = replace(model)
    ) |>
    mutate(model = factor(model, levels = rev(lvl)))

# Plot the time of running a model
p <- df.long |> 
    ggplot(aes(p, mu, color = model)) + 
    geom_line() + 
    geom_point() +
    scale_x_continuous(breaks = 1:12) +
    theme(text = element_text(size = 20)) +
    labs(color = "Model degree", x = "p", y = "time (s)", title = "Average runtime of mpi model run with p threads")

ggsave("mpi_runtime.png", p)

# Now plot the speedup as a function of the number of threads
p <- df.long |> 
    ggplot(aes(p, k, color = model)) + 
    geom_line() + 
    scale_x_continuous(breaks = 1:12) +
    theme(text = element_text(size = 20)) +
    labs(x = "p", y = "speedup", title = "Average speedup factor of mpi model run with p threads") +
    geom_line(aes(x = p, y = p), col = "black", linetype = "dashed")

ggsave("mpi_speedup.png", p)

p <- df.long |> 
    mutate(efficiency = k / p) |>
    ggplot(aes(p, efficiency, col = model)) +
    geom_smooth() +
    geom_point() +
    labs(color = "Model degree") +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(limits = c(0, 1.15)) +
    theme(text = element_text(size = 20)) +
    geom_line(aes(x = rep(1:12, 3), y = rep(1, 36)), col = "black") 

ggsave("mpi_efficieny.png", p)
