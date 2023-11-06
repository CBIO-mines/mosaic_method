library(ggplot2)
library(dplyr)
library(readr)
if(!require(reticulate)) {
  install.packages("reticulate", repos = "https://cloud.r-project.org/")
  library(reticulate)
}
parsefit <- import("parsefit")


mus <- 5e-9
muc <- 6e-11
delta <- 0.55
dr <- 0.1

mh_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[1]])
}
mc_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[2]])
}
sum_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[1]] + parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[2]])
}

binned_match_df <- read_csv("./binned_mld.csv")
fitted_params <- read_csv("./fitted_params.csv")
par1_fit <- unlist(c(fitted_params[1, "log10tau"], fitted_params[1, "log10rho"]))
list_args_mcmh <- list(par1 = par1_fit, dr = dr, mus = mus, muc = muc, d = delta, L0 = as.numeric(fitted_params[1, "L0"]))

fitted_curve_plot <- ggplot(binned_match_df, aes(x = match_length, y = freq)) +
  geom_point() +
  geom_function(fun = mc_fun, color = "blue", args = list_args_mcmh) +
  geom_function(fun = mh_fun, color = "red", args = list_args_mcmh) +
  geom_function(fun = sum_fun, color = "black", args = list_args_mcmh) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "log10(Match length)", y = "log10(Frequency)")

ggsave("plot_fig2.png", fitted_curve_plot)
