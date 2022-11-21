# Compare the error differences between least squares and numerical integration
library(tidyverse)
library(purrr)
library(pracma)
library(readr)
library(viridis)

old_models <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

old_model_error <- "csv/old_model_error.csv"

new_csv <- function(degree) {
    paste("csv/diff/diff_small_", degree, ".csv", sep = "")
}

csvs <- map_chr(old_models, new_csv) 

df <- as_tibble(read.table(csvs[1], header = TRUE))

for (i in 1:(length(old_models) - 1)) {
    tmp <- as_tibble(read.table(csvs[i + 1], header = TRUE))
    df <- df |> add_column(tmp[[1]], .name_repair = "unique")
}
# Load up all of these values
df <- df |> mutate(n = 1:nrow(df))


col_names <- map_chr(old_models, function(x) { paste("l", x, sep = "")})

colnames(df) <- c(col_names, "n")


df.mean <- df |> mutate(across(everything(), .fn = function(x) { mean(abs(x))})) |> head(1)
df.long <- df |> pivot_longer(c(1:length(old_models)), values_to = "diff", names_to = "degree")

x_mean <- unlist(df.mean[1,1:(ncol(df.mean) - 1)])

# Compute the means of the new model computation
df.mean <- tibble(degree = old_models, mu = x_mean)

df.old <- as_tibble(read_csv(old_model_error))

df.mean <- df.mean |> 
    mutate(mu_new = mu, mu_old = df.old$avg_error) |> 
    select(degree, mu_old, mu_new) |>
    pivot_longer(c(2:3), values_to = "mu", names_to = "model") |>
    mutate(model = factor(model))

cp <- c("#000000", "#FF0000")
labls <- c("least squares", "quadrature")

cp <- rev(cp)
labls <- rev(labls)

# names(cp) <- labls

p <- df.mean |> 
    ggplot(aes(degree, mu, col = model, group = model)) + 
    scale_color_manual(labels = labls, values = cp) +
    # scale_fill_manual(values = cp) +
    # scale_color_discrete(labels = c("least squares", "quadrature")) +
    # scale_color_viridis(option = "magma", discrete = TRUE) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = seq(0, 20, 2)) +
    scale_y_continuous(limits = c(0, 2500), breaks = seq(0, 2500, 250)) +
    # theme_gray() +
    labs(color = "Model", y = "Average abs errror") +
    theme(text = element_text(size = 30),
        panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    


ggsave("diff_error.png", p)