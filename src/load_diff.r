# Write a function that loads in all the differences of a data set so that we can compute the statistics about the differences
library(purrr)
library(tidyverse)

as_str <- function(x) {
    paste("", x, sep = "")
}

# dataset_size: "small" | "med" | "hi"
# degrees: a vector of degree sizes
load_diff <- function(dataset_size, degrees) {

    # Start off by using a custom naming function to find the csv files
    get_csv_path <- function(degree) {
        paste("csv/diff/diff_", dataset_size, "_", degree, ".csv", sep = "")
    }

    csvs <- map_chr(degrees, get_csv_path)
    print(csvs)

    df <- as_tibble(read.table(csvs[1], header = TRUE))

    # for (i in 1:(length(degrees) - 1)) {
    #     tmp <- as_tibble(read.table(csvs[i + 1], header = TRUE))
    #     df <- df |> add_column(tmp[[1]], .name_repair = "unique")
    # }

    # Read all of the csvs
    df_all <- map(csvs, function(x) { as_tibble(read_csv(x)) })
    df_all.first <- df_all[[1]]

    for (i in 1:(length(degrees) - 1)) {
        df_all.first <- df_all.first |> add_column(df_all[[i + 1]]$diff, .name_repair = "unique")
    }

    new_colnames <- map_chr(degrees, function(d) { paste("L", d, sep = "")})
    colnames(df_all.first) <- new_colnames
    
    all_means <- df_all.first |> mutate(across(everything(), function(x) { mean(abs(x)) }))
    df_mean <- tibble(degree = degrees, mu = unlist(all_means[1,]))
    # Compute the means

    L_names <- map_chr(degrees, function(x) { paste("L", x, sep = "")})

    df_all.first.long <- df_all.first |> 
        mutate(n = 1:nrow(df_all.first)) |> 
        pivot_longer(c(1:length(degrees)), names_to = "model", values_to = "error") |> 
        mutate(model = factor(model, levels = L_names))

    list(mean = df_mean, err = df_all.first, long = df_all.first.long)
}

# ========================== Main program =============================#
med_degrees <- c(0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800)
small_degrees <- c(0, 2, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 50, 100, 150, 175, 200, 215, 250, 300)


df.med <- load_diff("med", med_degrees)
df.small <- load_diff("small", small_degrees)

df.med$mean |> 
    ggplot(aes(degree, mu)) + 
    geom_line()  
    # geom_point(aes(x = df.small$mean$degree, y = df.small$mean$mu))

# I want to somehow coerce the small data set to the large one

small_pad <- vector("numeric", length(med_degrees))

for (i in 1:length(med_degrees)) {
    # if med_degrees[i] is in the small degrees, fill the appropriate value of df.small$mean
    mean_vec <- df.small$mean$mu
    ind <- small_degrees == med_degrees[i]
    if (sum(ind) == 0) {
        small_pad[i] <- NaN
    } else {
        small_pad[i] <- mean_vec[ind]
    }
}

p <- df.small$mean |> 
    ggplot(aes(degree, mu)) + 
    geom_line(size = 1) + 
    geom_point(aes(col = factor(degree)), size = 5) +
    theme(text = element_text(size = 20), legend.position = "none") +
    labs(y = "Average Absolute Erorr (AAE)", x = "Model degree")

ggsave("small_err.png", p)

p <- df.med$mean |> 
    ggplot(aes(degree, mu)) + 
    geom_line(aes(y = small_pad), linetype = "dashed", size = 1, color = "red") +
    geom_line(size = 1) + 
    geom_point(aes(col = factor(degree)), size = 5) +
    theme(text = element_text(size = 20), legend.position = "none") +
    labs(y = "Average Absolute Erorr (AAE)", x = "Model degree") +
    scale_y_continuous(breaks = seq(0, 2000, 250)) 

ggsave("med_err.png", p)

p <- df.med$long |> 
    filter(model == c("L0", "L10", "L200") | model == c("L100", "L400", "L600") | model == c("L800", "L50", "L700")) |>
    ggplot(aes(error, fill = model)) +
    # geom_histogram(bins = 1000) +
    geom_density() +
    facet_wrap(~ model)
    # filter(model == c("L0", "L10", "L50", "L100", "L200", "L400", "L600", "L800"))
ggsave("med_err_density.png", p, width = 10, height = 10)



df.med$mean |> 
    ggplot(aes(degree, mu)) + 
    geom_line() +
    geom_line(aes(y = small_pad), col = "red")

p <- df.med$long |> 
    ggplot(aes(error, col = model, fill = model)) +
    coord_flip() + 
    geom_boxplot(alpha = 0.7) +
    theme(text = element_text(size = 30))

ggsave("med_error_box.png", p)

p <- df.small$long |> 
    ggplot(aes(error, col = model, fill = model)) +
    coord_flip() + 
    geom_boxplot(alpha = 0.7)

ggsave("small_error_box.png", p)
# Compute the average absolute error 
