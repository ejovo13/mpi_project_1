# Statistical analysis of the difference values computed by our model
# library(tidyverse)
library(pracma)
library(tikzDevice)

library(ggplot2)
library(tibble)
library(purrr)
library(dplyr)
library(tidyr)


make_name <- function(degree) {
    paste("csv/diff/diff", degree, "small.csv", sep = "_")
}

# sample_values <- c(5, 10, 20, 30, 50, 100, 150, 175, 200, 300, 400, 500)
sample_values <- c(5, 10, 20, 30, 50, 100, 150, 175, 200, 300)

csv_files <- vector(mode = "character", length = length(sample_values))

for (i in seq_along(sample_values)) {
    csv_files[i] <- make_name(sample_values[i])
}

df <- tibble(read.csv(csv_files[1]))

for (i in seq_len(length(csv_files) - 1)) {
    # print(i)
    tmp <- tibble(read.csv(csv_files[i + 1]))
    df <- df |> add_column(tmp[[1]], .name_repair = "unique")
}


new_col_names <- purrr::map_chr(sample_values, function(x) { paste("L", x, sep = "") } )
colnames(df) <- new_col_names

plot_stuff <- function() {



    # df.long <- df |> mutate(n = 1:64800)
    # df.long <- pivot_longer(df.long, c(1:length(csv_files)), names_to = "model_degree", values_to = "err")


    df.long <- df |> 
        mutate(n = 1:64800) |>
        pivot_longer(c(1:length(csv_files)), names_to = "model_degree", values_to = "err") |>
        mutate(model_degree = factor(model_degree, levels = new_col_names))


    # Density
    # df.long |> ggplot(aes(err)) + geom_density() + facet_wrap(~ model_degree)

    # Histogram
    # df.long |> ggplot(aes(err)) + geom_histogram() + facet_wrap(~ model_degree)

    # tikz(file = "plot_test.text", width = 5, height = 5)

    # Boxplot
    df.long |> 
        ggplot(aes(err, col = model_degree)) + 
        geom_boxplot() +
        coord_flip() 

    ggsave("diff_boxplot.png")
    # dev.off()


    df.abs <- df |> mutate(across(.cols = everything(), .fns = abs))

    # Let's extract the standard deviation and the mean of the absolute error

    df.mean <- df.abs |> mutate(across(.cols = everything(), .fns = mean))
    df.sd   <- df.abs |> mutate(across(.cols = everything(), .fns = sd))


    x_mean <- unlist(df.mean[1,])
    x_std  <- unlist(df.sd[1,])

    df.stat <- tibble(order = factor(new_col_names, levels=new_col_names), mu = x_mean, sigma = x_std)
}

# Bootstrap 
# So currently we have a data set that contains the differences at each point for a variety of
# models. We'd like to study how a different number of samples caputes the true mean for L200

# We could use the frequentist approach and construct confidence intervals by setting a certain value of
# n. Just to start, let's set n = 1000 and see how our confidence intervals perform

n <- 5000 
abs_200 <- abs(df$L200)
mu_200 <- mean(abs_200)

# sample n points from L200 and estimate the mean
x <- sample(abs_200, n)
x_mean <- mean(x)

# Let's compute the sample mean n_sample times and see the distribution
n_sample <- 1000
sample_means <- vector("numeric", n_sample)

for (i in 1:n_sample) {
    sample_means[i] <- mean(sample(abs_200, n))
}

df.sample_means <- tibble(sample_means)

# df.sample_means |> ggplot(aes(sample_means)) + geom_histogram()
# df.sample_means |> ggplot(aes(sample_means)) + geom_density()

# Construct a 95% confidence interval
# The returned interval is a list with fields
# $lower and $upper
# confidence is a number between 0 and 1
ci <- function(population, n, confidence) {
    q <- 1 - (1 - confidence) / 2
    z <- qnorm(q)
    x <- sample(population, n)
    x_bar <- mean(x)
    s <- sd(x)
    plus_or_min <- z * s / sqrt(n)
    list(lower = x_bar - plus_or_min, upper = x_bar + plus_or_min)
}

# Function that returns true if the 95% confidence interval captured the true mean mu_200
confidence <- 0.99

is_mean_in_ci <- function() {
    interv <- ci(abs_200, n, confidence)
    interv$lower <= mu_200 && mu_200 <= interv$upper
}

# Now run 1000 confidence intervals!
n_interval <- 10000
interval_success <- vector("numeric", n_interval)

for (i in 1:n_interval) {
    interval_success[i] <- is_mean_in_ci()
}

success_rate <- sum(interval_success) / n_interval

print(success_rate)