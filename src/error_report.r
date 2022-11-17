# Visualize the data held in the files error_report.csv

library(pracma)
library(tidyverse)

filename <- "csv/small_error_report.csv"

df <- tibble(read.csv(filename))

# TODO clean this up for a latex graph

df |> ggplot(aes(l, avg_error)) + geom_line()