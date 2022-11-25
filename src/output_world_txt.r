# Output an ascii text file of the world
library(tidyverse)

brightness <- paste("$@B%8&WM#", "*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\\|", "()1{}[]?-_+~<>i!lI;:,", '"', "^`'.", sep = "")

# scale is from 1 to 69
n <- nchar(brightness)

# load the data in df.small and then convert into 
# a text file

load("discrete_grid.RData")
r <- df.small$r

LOWER_BOUND <- -4000
UPPER_BOUND <- 5000

r[r < LOWER_BOUND] <- LOWER_BOUND
r[r > UPPER_BOUND] <- UPPER_BOUND

m <- min(r)
M <- max(r)


# minimizing via the min and max isnt good enough

r_ <- (r - m) / (M - m)

df.small <- df.small |> mutate(r_ = r_)

min(r - m)

rth_char <- function(r) {
    # index = (floor((1 - r) * n) + 1)
    index = floor((1 - r) * n - 1) 
    c <- substr(brightness, index, index)
    if (c == "") {
        c <- "$"
    }
    return (c)
}

r_char <- map_chr(r_, rth_char)

# print the world

nth <- 60; nph <- 120 

for (t in seq(0, nth, 2)) {
    for (p in 1:nph) {
        cat(r_char[nth * p + t])
    }
    cat("\n")
}

hist(r_)