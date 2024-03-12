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
      fitted_params = list.files("./results_refseq_test/fitted_params/")
    ),
    output = list(),
    params = list(
      fitted_params_dir = "results_refseq_test/fitted_params/",
      results_dir = "results_refseq_test/",
      min_r_infl = 16,
      taxon_csv = "bacillaceae_taxon_2.csv",
      cluster_name = "genus"
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
library(purrr)

source("utils.R")


read_fits <- function(x) {
  tmp_df <- read_csv(paste0(params_dir, x), col_types = cols(.default = col_double()))
  comp <- str_split_1(x, "_vs_")
  comp[2] <- str_split_1(comp[2], "_fitted")[1]
  tmp_df$comp <- paste(comp, collapse = "_")
  tmp_df$cluster_1 <- comp[1]
  tmp_df$cluster_2 <- comp[2]
  tmp_df
}


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


per_infl_clust <- function(r_infl, clust_level, n_clusters) {
  # how many pairwise distances are actually fitted
  r_infl %>%
    filter(bacterias_in_reference(., clust_level, "cluster")) %>%
    group_by(infl_exist) %>%
    tally %>%
    filter(infl_exist == "yes") %>%
    # -1 because no distance to self, obv
    mutate(per = n/(n_clusters - 1)) %>%
    pull(per)
}

params_dir <- snakemake@params[["fitted_params_dir"]]
results_dir <- snakemake@params[["results_dir"]]
cluster_name <- snakemake@params[["cluster_name"]]
# read fitted params -----------------------------------------------------------
list_fits_files <- list.files(params_dir)
fit_res <- map_df(list_fits_files, ~ read_fits(.x))

# read taxon_csv
taxon_csv <- read_csv(snakemake@params[["taxon_csv"]])

# does inflexion exist ?
r_infl <- fit_res %>%
  rowwise() %>%
  mutate(r_infl =  get_infl_exist(log10tau, log10rho, 1:1000, L0)) %>%
  mutate(infl_exist = ifelse(r_infl > snakemake@params[["min_r_infl"]], "yes", "no")) %>%
  ungroup()

write_csv(r_infl, paste0(results_dir, "inflexion_exists.csv"))

# How much does a given cluster have regime change with others ?

present_clusters <- c(r_infl$cluster_1, r_infl$cluster_2) %>% unique()

clusters <- taxon_csv %>%
  select(all_of(cluster_name)) %>%
  distinct %>%
  pull(all_of(cluster_name))

n_clusters <- length(clusters)
per_infl <- map_dbl(clusters, ~ per_infl_clust(r_infl, .x, n_clusters))
infl_by_cluster <- tibble({{cluster_name}} := clusters, "per_infl" = per_infl)

write_csv(infl_by_cluster, paste0(results_dir, "inflexion_by_cluster.csv"))
