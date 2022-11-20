# Let's show the 2-d image of the med or small data set.
library(tidyverse)
library(readr)

csv_small <- "csv/ETOPO1_small.csv"
df <- as_tibble(read.table(csv_small), header = FALSE, sep = "\t")

colnames(df) <- c("", "phi")