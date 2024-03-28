library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(reticulate)
use_condaenv("test_florian")


fit <- import("fit")


mus <- snakemake@params[["mus"]]
muc <- snakemake@params[["muc"]]
delta <- snakemake@params[["delta"]]
results_dir <- snakemake@params[["results_dir"]]
# the small difference for numerical derivation
dr <- 0.1

mh_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(fit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[1]])
}
mc_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(fit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[2]])
}
sum_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(fit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[1]] + fit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[2]])
}

binned_match_df <- read_csv(snakemake@input[["binned_mld"]])
fitted_params <- read_csv(snakemake@input[["fitted_params"]])
species <- str_split_1(snakemake@params[["species"]], ",")
par1_fit <- unlist(c(fitted_params[1, "log10tau"], fitted_params[1, "log10rho"]))
list_args_mcmh <- list(par1 = par1_fit, dr = dr, mus = mus, muc = muc, d = delta, L0 = as.numeric(fitted_params[1, "L0"]))

lim_freq <- min(binned_match_df$freq[binned_match_df$freq > 0])/10
fitted_curve_plot <- ggplot(binned_match_df, aes(x = match_length, y = freq)) +
  geom_point() +
  geom_function(fun = mc_fun, color = "blue", args = list_args_mcmh) +
  geom_function(fun = mh_fun, color = "red", args = list_args_mcmh) +
  geom_function(fun = sum_fun, color = "black", args = list_args_mcmh) +
  scale_x_log10() +
  scale_y_log10(limits = c(lim_freq, NA)) +
  labs(x = "Match length", y = "Frequency")

ggsave(paste0(results_dir, "fig2_plots/", paste0(species[1], "_vs_", species[2], "_", "plot_fig2.png")), fitted_curve_plot)
