if(interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots = c(
      input = "list",
      output = "list",
      params = "list",
      wildcards = "list",
      threads = "numeric",
      log = "list",
      resources = "list",
      config = "list",
      rule = "character",
      bench_iteration = "numeric",
      scriptdir = "character",
      source = "function"
    )
  )
  snakemake <- Snakemake(
    input = list(
      fitted_params = list.files("./results_refseq_real/fitted_params/")
    ),
    output = list(),
    params = list(
      fitted_params_dir = "results_refseq_real/fitted_params/",
      analyse_dir = "results_refseq_real/analyse_comparisons/",
      results_dir = "results_refseq_real/",
      min_r_infl = 30
        ),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list(),
    rule = "",
    bench_iteration = 1,
    scriptdir = "",
    source = function(...) {{ wd <- getwd()
      setwd(snakemake@scriptdir)
      source(...)
      setwd(wd) }}
  )
}

library(RColorBrewer)
library(stringr)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)
library(reticulate)
use_condaenv("test_florian")
parsefit <- import("parsefit")
library(purrr)
params_dir <- snakemake@params[["fitted_params_dir"]]
results_dir <- snakemake@params[["results_dir"]]
analyse_dir <- snakemake@params[["analyse_dir"]]

read_res <- function(x) {
  tmp_df <- read_csv(paste0(analyse_dir, x), col_types = cols(.default = col_double()))
  tmp_df$comp <- paste0(colnames(tmp_df)[1:2], collapse = "_vs_")
  tmp_df <- tmp_df %>%
    rename(bac_1_genomes_entropy = 1, bac_2_genomes_entropy = 2)
  tmp_df
}

read_fits <- function(x) {
  tmp_df <- read_csv(paste0(params_dir, x), col_types = cols(.default = col_double()))
  tmp_df$comp <- paste(str_split_1(x, "_")[1:3], collapse = "_")
  tmp_df
}


mc_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[2]])
}


mh_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[1]])
}


# Overall stats ----------------------------------------------------------------
list_res <- list.files(analyse_dir)[str_detect(list.files(analyse_dir), "res.csv")]
conc_res <- map_df(list_res, ~ read_res(.x)) %>%
  separate_wider_delim(cols = "comp", delim = "_vs_", names = c("bac_1", "bac_2"), cols_remove = FALSE)
conc_res <- conc_res %>%
  mutate(infl_per = (n_comp_infl)/n_comp_tot)
promising_res <- conc_res %>%
  filter(n_comp_infl != 0)
unpromising_res <- conc_res %>%
  filter(n_comp_infl == 0)


# read fitted params
list_fits_files <- list.files(params_dir)
fit_res <- map_df(list_fits_files, ~ read_fits(.x))

joined_fits_analyse <- left_join(conc_res, fit_res, by = "comp")

if(interactive())
  infl_per_vs_tau<- ggplot(joined_fits_analyse, aes(x = infl_per, y = log10tau)) +
    geom_point(alpha = 0.3)


# inflexion on whole species fits

get_infl_exist <- function(log10tau, log10rho, r, L0) {
  dr <- 0.1
  mus <- 5e-9
  muc <- 6e-11
  d <- 0.55
  mh <- mh_fun(
    r = np_array(r, dtype = "float64"),
    par1 = np_array(c(log10tau, log10rho), dtype = "float64"),
    dr = np_array(0.1, dtype = "float64"),
    mus = np_array(5e-9, dtype = "float64"),
    muc = np_array(6e-11, dtype = "float64"),
    d = np_array(0.55, dtype = "float64"),
    L0 = np_array(L0, dtype = "float64")
  )
  mc <- mc_fun(
    r = np_array(r, dtype = "float64"),
    par1 = np_array(c(log10tau, log10rho), dtype = "float64"),
    dr = np_array(0.1, dtype = "float64"),
    mus = np_array(5e-9, dtype = "float64"),
    muc = np_array(6e-11, dtype = "float64"),
    d = np_array(0.55, dtype = "float64"),
    L0 = np_array(L0, dtype = "float64")
  )
  diff <- mh - mc
  if(all(diff > 0))
    return(0)
  return(r[max(which(diff < 0))])
}


r_infl <- joined_fits_analyse %>%
  rowwise() %>%
  mutate(r_infl =  get_infl_exist(log10tau, log10rho, 1:10000, L0)) %>%
  mutate(infl_exist = ifelse(r_infl > snakemake@params[["min_r_infl"]], "yes", "no")) %>%
  ungroup()

if(interactive())
  boxplot_comp_inflexion <- ggplot(r_infl, aes(x = infl_exist, y = n_comp_tot)) +
    geom_boxplot() +
    scale_y_log10()

infl_exist <- r_infl %>%
  filter(infl_exist == "yes") %>%
  select(comp)

write_csv(infl_exist, paste0(results_dir, "inflexion_exists.csv")
