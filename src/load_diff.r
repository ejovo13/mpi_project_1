# Write a function that loads in all the differences of a data set so that we can compute the statistics about the differences
library(purrr)
library(tidyverse)

med_degrees <- c(0, 1, 2, 3, 4, 5, 10, 20, 30)

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


    df_all

}

df.med <- load_diff("med", med_degrees)