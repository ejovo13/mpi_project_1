library(tidyverse)
library(pracma)


filename <- "testcoeff_series.csv"

df <- tibble(read.csv(filename))

# plot the trajectory of c00

# df |> ggplot(aes(t, c0_0)) + geom_point() + geom_point(aes(df$t, df$c1_0)) 

# Time and a final column get superfluously picked up
two_ll <- ncol(df) - 2
ll <- two_ll / 2
lmax <- sqrt(two_ll + 0.25) - 3 / 2 

# Go ahead and pivot longer

end_cos <- ll + 1
end <- ll + 1
start_sin <- ll + 2
end_sin <- two_ll + 1

df.clm <- df |> select(c(1:end_cos, ncol(df))) # still have mse
df.slm <- df |> select(c(1, c(start_sin:end_sin), ncol(df)))

row_c_f <- df.clm[nrow(df.clm), 2:end]
row_s_f <- df.slm[nrow(df.slm), 2:end]
# Get the order of row_c_f and row_s_f

row_c_order <- order(row_c_f, decreasing = TRUE)
row_s_order <- order(row_s_f, decreasing = TRUE)

c_lev <- colnames(df.clm)
c_lev <- c_lev[1 + row_c_order]

s_lev <- colnames(df.slm)
s_lev <- s_lev[1 + row_s_order]
# c_lev <- c_lev[row_c_order]

# Rearrange the order of the cols
df.clm <- df.clm |> select(c(1, 1 + row_c_order), end + 1)
df.slm <- df.slm |> select(c(1, 1 + row_s_order), end + 1)

# Pivot and set correct factor level
df.clm <- df.clm |> 
    pivot_longer(c(2:end), names_to = "c", values_to = "val") |>
    mutate(c = factor(c, levels = c_lev))

df.slm <- df.slm |>
    pivot_longer(c(2:end), names_to = "s", values_to = "val") |>
    mutate(s = factor(s, levels = s_lev))

p_mse <- df |> ggplot(aes(t, mse)) + geom_line()

ggsave("mse_series.png", plot = p_mse)

# Plotting
p_c <- df.clm |> 
    ggplot(aes(t, val, col = c)) + 
    geom_line() +
    scale_fill_discrete(breaks = c_lev) + 
    theme(legend.position="none")

p_s <- df.slm |> 
    ggplot(aes(t, val, col = s)) + 
    geom_line() +
    scale_fill_discrete(breaks = s_lev) + 
    theme(legend.position="none")

ggsave("Clm_series5_s.png", plot = p_c)
ggsave("Slm_series5_s.png", plot = p_s)

p_c

# Plot the mse



# Let's try and reorder the legend so that the highest values are at the top of the legend
# access the last row of data
