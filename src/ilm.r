# Play around with the functions needed to convert an iteration value to a (l, m) pair
library(purrr)

i_to_lmax <- function(i) {

    two_LL <- (i + 1) * 2

    lmax <- sqrt(two_LL + 0.25) - 1.5
    ceiling(lmax)
}

i_to_m <- function(i) {

    # first compute lmax
    lmax <- i_to_lmax(i)
    # Then subtract the previous LL - 1
    LL_prev <- (lmax) * (lmax + 1) / 2
    start <- i - LL_prev
    start
}

i <- 0:10

print(map_dbl(i, i_to_lmax))
print(map_dbl(i, i_to_m))