library(pracma)
library(tidyverse)
library(readr)

med_values <- c(0, 5, 10, 20, 50, 100, 150, 200, 300)


med_name <- function(degree) {
    paste("csv/models/sph_med_", degree, ".txt", sep = "")
}

hi_name <- function(degree) {
    paste("csv/models/sph_hi_", degree, ".txt", sep = "")
}

plot_C <- function(degree) {

    df <- tibble(read.table(med_name(degree), sep = "\t"))
    colnames(df) <- c("l", "m", "Clm", "Slm")
    df <- df |> mutate(i = 1:nrow(df))
    df |> ggplot(aes(i, Clm)) + geom_point()

}

plot_S <- function(degree) {

    df <- tibble(read.table(med_name(degree), sep = "\t"))
    colnames(df) <- c("l", "m", "Clm", "Slm")
    df <- df |> mutate(i = 1:nrow(df))
    df |> ggplot(aes(i, Slm)) + geom_point()

}


# load in small med and hi
sets <- c("small", "med", "hi")
degrees <- c(0, 5, 10, 20, 50, 100, 150, 200, 300)

# dataset_size is a value in {"small", "med", "hi"}
create_coeff_set <- function(dataset_size) {

    make_name <- function(degree) {
        paste("csv/models/sph_", dataset_size, "_", degree, ".txt", sep = "")
    }

    list_df <- map(degrees, function(degree) { 
        # Load in the data set
        df <- tibble(read.csv(make_name(degree), sep = "\t", header = FALSE))
        colnames(df) <- c("degree", "order", "Clm", "Slm")
        df
    })

    list_df
}

# Create three coeff sets
coeff <- map(sets, create_coeff_set)

coeff.small <- coeff[[1]]
coeff.med   <- coeff[[2]]
coeff.hi    <- coeff[[3]]

# let's join these into one single list of length(orders)
coeff.all <- coeff.small

for (i in 1:(length(coeff) - 1)) {

    for (j in 1:length(coeff[[1]])) {
        coeff.all[[j]] <- coeff.all[[j]] |> add_column(coeff[[i + 1]][[j]][3], .name_repair = "unique") # add Clm to the base
        coeff.all[[j]] <- coeff.all[[j]] |> add_column(coeff[[i + 1]][[j]][4], .name_repair = "unique") # add Slm to the base
    }

    # Let's rename the columns too
    colnames(coeff.all[[j]]) <- c("degree", "order", "Clm_s", "Slm_s", "Clm_m", "Slm_m", "Clm_h", "Slm_h")

}

# difference between previous data set
coeff.300 <- coeff.all[[9]]
df.diff.300 <- coeff.300 |> 
    mutate(Clm_ms = Clm_m - Clm_s) |>
    mutate(Clm_hm = Clm_h - Clm_m) |>
    select(degree, order, Clm_ms, Clm_hm)

# The difference between Medium and small is substantial
df.diff.300 |> ggplot(aes(x = 1:nrow(df.diff.300), y = Clm_ms)) + geom_point(alpha = 0.5)
# But the difference between hi and medium is far more uniform
df.diff.300 |> ggplot(aes(x = 1:nrow(df.diff.300), y = Clm_hm)) + geom_point(alpha = 0.5)

# coeff <- list(list())

# for (s in sets) {
#     coeff[s] <- create_coeff_set(s)
# }




csv_files_med <- map_chr(med_values, med_name)

df.med <- tibble(read.table(csv_files_med[1], sep = "\t"))

df.300 <- tibble(read.table(med_name(300), sep = "\t"))
colnames(df.300) <- c("l", "m", "Clm", "Slm")
df.300 <- df.300 |> mutate(i = 1:nrow(df.300))

df.300 |> ggplot(aes(i, Clm)) + geom_point()

#

# Let's see how for order = 200, how do the coefficients change
# eventually for all orders?

# First step, need to load in all of the coefficients.
# Let's load in 3 data sets, small, med, and hi

# for (i in seq_len(length(csv_files_med) - 1)) {

#     tmp <- tibble(read.csv(csv_files_med[i + 1]))
#     df.med <- df.med |> add_column(tmp[[1]], .name_repair = "unique")

# }