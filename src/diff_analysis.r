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
    # paste("csv/diff/diff", degree, "small.csv", sep = "_")
    paste("csv/diff/diff_small_", degree, ".csv", sep = "")
}

med_name <- function(degree) {
    paste("csv/diff/diff_med_", degree, ".csv", sep = "")
}

hi_name <- function(degree) {
    paste("csv/diff/diff_hi_", degree, ".csv", sep = "")
}

# sample_values <- c(5, 10, 20, 30, 50, 100, 150, 175, 200, 300, 400, 500)
# sample_values <- c(5, 10, 20, 30, 50, 100, 150, 175, 200, 300)
sample_values <- c(0, 5, 10, 20, 50, 100, 150, 200, 300)

csv_files <- vector(mode = "character", length = length(sample_values))
csv_files_med <- map_chr(sample_values, med_name)
csv_files_hi <- map_chr(sample_values, hi_name)

for (i in seq_along(sample_values)) {
    csv_files[i] <- make_name(sample_values[i])
    
}

df <- tibble(read.csv(csv_files[1]))
df.med <- tibble(read.csv(csv_files_med[1]))
df.hi <- tibble(read.csv(csv_files_hi[1]))

for (i in seq_len(length(csv_files) - 1)) {
    # print(i)
    tmp <- tibble(read.csv(csv_files[i + 1]))
    df <- df |> add_column(tmp[[1]], .name_repair = "unique")
}

for (i in seq_len(length(csv_files_med) - 1)) {

    tmp <- tibble(read.csv(csv_files_med[i + 1]))
    df.med <- df.med |> add_column(tmp[[1]], .name_repair = "unique")

}

for (i in seq_len(length(csv_files_hi) - 1)) {

    tmp <- tibble(read.csv(csv_files_hi[i + 1]))
    df.hi <- df.hi |> add_column(tmp[[1]], .name_repair = "unique")

}


new_col_names <- purrr::map_chr(sample_values, function(x) { paste("L", x, sep = "") } )
colnames(df) <- new_col_names
colnames(df.med) <- new_col_names
colnames(df.hi) <- new_col_names



plot_stuff <- function(df, suffix = "", save = TRUE) {



    # df.long <- df |> mutate(n = 1:64800)
    # df.long <- pivot_longer(df.long, c(1:length(csv_files)), names_to = "model_degree", values_to = "err")


    df.long <- df |> 
        mutate(n = 1:nrow(df)) |>
        pivot_longer(c(1:ncol(df)), names_to = "model_degree", values_to = "err") |>
        mutate(model_degree = factor(model_degree, levels = new_col_names))


    # Density
    # df.long |> ggplot(aes(err)) + geom_density() + facet_wrap(~ model_degree)

    # Histogram
    # df.long |> ggplot(aes(err)) + geom_histogram() + facet_wrap(~ model_degree)

    # tikz(file = "plot_test.text", width = 5, height = 5)

    # Boxplot
    p <- df.long |> 
        ggplot(aes(err, col = model_degree)) + 
        geom_boxplot() +
        coord_flip() 

    print(p)
    if (save) {
        ggsave(paste("boxplot_", suffix, ".png", sep = ""), plot = p)
    }

    # Let's draw the densities on the same plot
    p <- df.long |>
        ggplot(aes(err, col = model_degree)) +
        geom_density()

    print(p)
    if (save) {
        ggsave(paste("density_", suffix, ".png", sep = ""), plot = p)
    }

    p <- df.long |>
        ggplot(aes(err, fill = model_degree)) +
        geom_histogram(bins = 1000)

    if (save) {
        ggsave(paste("histogram_", suffix, ".png", sep =""), plot = p)
    }
    # ggsave("diff_boxplot.png")
    # dev.off()

    ## PLOT the density in facets.
    p <- df.long |>
        ggplot(aes(err, fill = model_degree)) +
        geom_histogram(bins = 1000) +
        facet_wrap(~ model_degree)

    if (save) {
        ggsave(paste("faceted_histogram_", suffix, ".png", sep = ""), plot = p)
    }

    df.abs <- df |> mutate(across(.cols = everything(), .fns = abs))

    # Let's extract the standard deviation and the mean of the absolute error

    df.mean <- df.abs |> mutate(across(.cols = everything(), .fns = mean))
    df.sd   <- df.abs |> mutate(across(.cols = everything(), .fns = sd))
    df.min  <- df.abs |> mutate(across(.cols = everything(), .fns = min))
    df.max  <- df.abs |> mutate(across(.cols = everything(), .fns = max))


    x_mean <- unlist(df.mean[1,])
    x_std  <- unlist(df.sd[1,])
    x_min  <- unlist(df.min[1,])
    x_max  <- unlist(df.max[1,])

    df.stat <- tibble(order = factor(new_col_names, levels=new_col_names), mu = x_mean, sigma = x_std, min = x_min, max = x_max)
    df.stat <- df.stat |> mutate(ord = sample_values)

    return (df.stat)
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
n_interval <- 10000
n <- 5000

test_confidence_interval <- function(population, n_interval, confidence) {
    
    is_mean_in_ci <- function() {
        interv <- ci(population, n, confidence)
        interv$lower <= mean(population) & mean(population) <= interv$upper
    }

    interval_success <- vector("numeric", n_interval)

    for (i in 1:n_interval) {
        interval_success[i] <- is_mean_in_ci()
    }

    success_rate <- sum(interval_success) / n_interval

    print(success_rate)
}
# Now run 1000 confidence intervals!
# test_confidence_interval(abs_200, n_interval, confidence)
df.small.stat <- plot_stuff(df, save = FALSE)

df.med.stat <- plot_stuff(df.med, save = FALSE)

df.hi.stat <- plot_stuff(df.hi, save = FALSE)

F_emp_med300 <- ecdf(df.med$L300)
y <- F_emp_med300(df.med$L300)

l300_med_sorted <- sort(df.med$L300)
# Take the mean of only the middle 95%
i_0 <- floor(0.025 * length(l300_med_sorted))
i_f <- floor(0.975 * length(l300_med_sorted))
abs_300 <- abs(l300_med_sorted)
sort_abs_300 <- sort(abs_300)
# Play around with df.med.300

df.med.300 <- df.med |> select(L300) |> mutate(n = 1:nrow(df.med))

df.med.300 |> ggplot(aes(L300)) + geom_density()

mean(abs(df.med.300$L300))


df.stat <- df.hi.stat |> 
    add_column(mu_med = df.med.stat$mu, mu_small = df.small.stat$mu) |>
    select(order = ord, mu_small, mu_med, mu_hi = mu) |>
    pivot_longer(c(2:4), names_to = "dataset_size", values_to = "mu")

p <- df.stat |> 
    ggplot(aes(order, mu, col = dataset_size)) + 
    scale_color_discrete(labels = c("high", "med", "small")) +
    geom_line() + 
    geom_point() + 
    theme(text = element_text(size = 30)) + 
    facet_wrap(~ dataset_size)

