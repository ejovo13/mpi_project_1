# Let's show the 2-d image of the med or small data set.
library(tidyverse)
library(readr)
library(pracma)
# library(imager)
library(viridis)
library(ggdark)
# library(raster)

data.small <- list(size = "small", N = 64800)
data.med   <- list(size = "med", N = 583200)
data.hi    <- list(size = "hi", N = 1080 * 2160)
# data.small <- list(size = "small", N = 64800)


process_dataset <- function(data, reduction = 6) {

    csv <- paste("csv/ETOPO1_", data$size, ".csv", sep = "")
    df <- as_tibble(read.table(csv), header = FALSE, sep = "\t")
    colnames(df) <- c("ph", "th", "r")

    n_th <- sqrt(data$N / 2)
    n_ph <- n_th * 2

    d_th <- pi / (n_th)
    d_ph <- 2*pi / (n_ph)

    # TH <- lin

    # Plot the enitre world
    p <- df |> ggplot(aes(ph, th, fill = r)) +
    scale_colour_brewer(palette = 3) +
    scale_fill_viridis_c(option="magma") + 
    geom_raster() +
    # theme_void() + 
    theme(
        legend.position = "bottom"
    )

    print(p)

    # Plot the top left quadrant
    df.top_left <- df |> filter(ph < 0, th > 0)

    # Plot one third of this data set
    # p <- df.top_left |> ggplot(aes(ph, th, fill = r)) +
    #     scale_colour_brewer(palette = 3) +
    #     scale_fill_viridis_c(option="magma") + 
    #     geom_raster() +
    #     # theme_void() + 
    #     theme(
    #         legend.position = "bottom"
    #     )
    # print(p)

    df |> mutate(i = 1:nrow(df)) |> filter(mod(i, reduction) == 0)

    # Let's imagine we are filling a nth x nph matrix column wise
    get_row <- function(i) {
        i0 <- i - 1 # convert to zero based indexing
        mod(i0, n_th) + 1
    }

    get_col <- function(i) {
        i0 <- i - 1
        floor(i / n_th) + 1
    }

    d_th <- d_th * reduction
    d_ph <- d_th * reduction


    df <- df |> 
        mutate(i = 1:nrow(df)) |> 
        filter(mod(get_row(i), reduction) == 0) |>
        filter(mod(get_col(i), reduction) == 0)

    # Plot a reduced world
    p <- df |> ggplot(aes(ph, th, fill = r)) +
    scale_colour_brewer(palette = 3) +
    scale_fill_viridis_c(option="magma") + 
    geom_tile(width = d_ph, height = d_th) +
    # theme_void() + 
    theme(
        legend.position ="none" 
    )

    ggsave("discretized_reduced.png", p, width = 10, height = 5)

    # print(p)
    df

}

# TH <- linspace(0, pi - d_th, n_th)
# PH <- linspace(0, 2 * pi - d_ph, n_ph)
# x <- matrix(df$r, nrow = n_th, ncol = n_ph, byrow = FALSE)

# x <- matrix(df$r, nrow = n_ph, ncol = n_th, byrow = TRUE)

# m_df <- tibble(as.data.frame(x))
# m_df |> 
#     ggplot(aes(x = x, y = y, fill = layer)) +
#     geom_raster() +
# invert_geom_defaults()

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


df <- process_dataset(data.small)
# process_dataset(data.med)

# use mesh grid of a smaller data set
# discretize -pi/2 and pi/2 into n_theta points

n_th <- 30 
n_ph <- 2 * n_th
TH <- linspace(-pi / 2, pi / 2, n_th)
PH <- linspace(-pi, pi, n_th * 2)

mat <- meshgrid(TH, PH)
th <- c(mat$X)
ph <- c(mat$Y)

p <- tibble(ph, th) |> ggplot(aes(ph, th)) + geom_point()
# ggsave("discretized_grid.png", p, width = 10, height = 5, units = "in")
