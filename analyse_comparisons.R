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
      species = "Bacillus,Virgibacillus",
      fitted_params_dir = "results_refseq_real/fitted_params/",
      full_mlds_dir = "results_refseq_real/full_mlds/",
      binned_mld_dir = "results_refseq_real/binned_mlds/",
      mock = "yes",
      results_dir = "results_refseq_real/",
      min_r_infl = 50
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

get_comp_number <- function(bac_1, bac_2, full_mlds_dir) {
  full_comp_file <- paste0(full_mlds_dir, bac_1, "_vs_", bac_2, "_full_mld_comp.csv")
  line_number <- str_split_1(system2("wc", args = c("-l", full_comp_file), stdout = TRUE), " ")[1] %>% as.numeric()
  return(line_number - 1)
}


bin_mld_r <- function(named_mld) {
  named_mld[is.na(named_mld)] <- 0
  correct_input <- data.frame("match_length" = names(named_mld) %>% as.double, "freq" = unname(named_mld)) %>%
    arrange(match_length)
  tmp_df <- parsefit$fit$bin_mld(
                           correct_input,
                           3,
                           35.5,
                           0.1,
                           1.0
                           )
  setNames(tmp_df$freq, tmp_df$match_length)
}


fit_params_r <- function(bin_mld, match_length) {
  res_opt <- parsefit$fit$fit_params(
                            opt_method = "dual-annealing",
                            init_pars = np_array(c(5, -5), dtype = "float64"),
                            empirical_mld = np_array(bin_mld, dtype = "float64"),
                            smal_dif = np_array(0.1, dtype = "float64"),
                            match_lengths = np_array(match_length, dtype = "float64"),
                            mus = np_array(5e-9, dtype = "float64"),
                            muc = np_array(6e-11, dtype = "float64"),
                            delta = np_array(0.55, dtype = "float64"),
                            L0 = np_array(L0, dtype = "float64")
                            )
  if(res_opt$success == TRUE) {
    list("tau" = res_opt$x[1], "rho" = res_opt$x[2])
  }
  else {
    list("tau" = 0, "rho" = 0)
  }
}


mc_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[2]])
}


mh_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(parsefit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[1]])
}


params_dir <- snakemake@params[["fitted_params_dir"]]
full_mlds_dir <- snakemake@params[["full_mlds_dir"]]
binned_mld_dir <- snakemake@params[["binned_mld_dir"]]
results_dir <- snakemake@params[["results_dir"]]
species <- str_split_1(snakemake@params[["species"]], ",")
min_r_infl <- snakemake@params[["min_r_infl"]]

output_dir <- paste0(results_dir, "analyse_comparisons/")
if(!dir.exists(output_dir))
  dir.create(output_dir)

# Read data
species_fit <- read_csv(paste0(params_dir, species[1], "_vs_", species[2], "_fitted_params.csv"))
L0 <- species_fit$L0
full_mld <- read_csv(paste0(full_mlds_dir, species[1], "_vs_", species[2], "_full_mld_comp.csv"), col_types = cols(.default = col_double(), "comp" = col_character()))[, -1] %>%
  column_to_rownames("comp") %>% as.matrix()
binned_mld <- read_csv(paste0(binned_mld_dir, species[1], "_vs_", species[2], "_binned_mld.csv")) %>%
  rename("r" = "match_length")



# Binning single mld
full_mld_list <- lapply(seq_len(nrow(full_mld)), function(i) full_mld[i, ])
names(full_mld_list) <- rownames(full_mld)
binned_comparisons_mld_list <- full_mld_list |> map(\(x) bin_mld_r(x))
binned_comparisons_df <- bind_cols(bind_rows(binned_comparisons_mld_list), "comp" = names(full_mld_list))
r <- as.double(colnames(select(binned_comparisons_df, -comp)))

# fitting a tau (and rho) to each single mld
fit_params_start_time <- Sys.time()
binned_comparisons_df <- binned_comparisons_df %>%
  rowwise(comp) %>%
  mutate(opt_pars = list(fit_params_r(c(c_across(where(is.numeric))), r)))
fit_params_stop_time <- Sys.time()
fit_params_full_time <- difftime(fit_params_stop_time, fit_params_start_time, units = "mins")
binned_comparisons_df <- binned_comparisons_df %>%
  mutate(tau = opt_pars[["tau"]], rho = opt_pars[["rho"]]) %>%
  ungroup() %>%
  select(-opt_pars)

# create data to plot all theoretical mlds
theoretical_mlds_df <- map_df(seq_len(nrow(binned_comparisons_df)), function(i) {
  tibble("r" = r,
         "mc" = mc_fun(
           r = np_array(r, dtype = "float64"),
           par1 = np_array(c(binned_comparisons_df$tau[i], binned_comparisons_df$rho[i]), dtype = "float64"),
           dr = np_array(0.1, dtype = "float64"),
           mus = np_array(5e-9, dtype = "float64"),
           muc = np_array(6e-11, dtype = "float64"),
           d = np_array(0.55, dtype = "float64"),
           L0 = np_array(L0, dtype = "float64")
         ),
         "mh" = mh_fun(
           r = np_array(r, dtype = "float64"),
           par1 = np_array(c(binned_comparisons_df$tau[i], binned_comparisons_df$rho[i]), dtype = "float64"),
           dr = np_array(0.1, dtype = "float64"),
           mus = np_array(5e-9, dtype = "float64"),
           muc = np_array(6e-11, dtype = "float64"),
           d = np_array(0.55, dtype = "float64"),
           L0 = np_array(L0, dtype = "float64")
         ),
         "tau" = binned_comparisons_df$tau[i],
         "comp" = binned_comparisons_df$comp[i],
         "rho" = binned_comparisons_df$rho[i]
         )
})

# compute the theoretical mld for the fitted params on the sum mld
sum_pars <- species_fit %>%
  select(log10tau, log10rho) %>%
  as.numeric
list_args_mcmh <- list(par1 = sum_pars, dr = 0.1, mus = 5e-9, muc = 6e-11, d = 0.55, L0 = L0)
mc_all <- tibble("r" = r, "mc_all" = mc_fun(
           r = np_array(r, dtype = "float64"),
           par1 = np_array(sum_pars, dtype = "float64"),
           dr = np_array(0.1, dtype = "float64"),
           mus = np_array(5e-9, dtype = "float64"),
           muc = np_array(6e-11, dtype = "float64"),
           d = np_array(0.55, dtype = "float64"),
           L0 = np_array(L0, dtype = "float64")
         ))
mh_all <- tibble("r" = r, "mh_all" = mh_fun(
           r = np_array(r, dtype = "float64"),
           par1 = np_array(sum_pars, dtype = "float64"),
           dr = np_array(0.1, dtype = "float64"),
           mus = np_array(5e-9, dtype = "float64"),
           muc = np_array(6e-11, dtype = "float64"),
           d = np_array(0.55, dtype = "float64"),
           L0 = np_array(L0, dtype = "float64")
         ))


negative_fitted_mc <- theoretical_mlds_df %>% filter(mc < 0) %>% pull(comp) %>% unique

lim_freq <- min(binned_mld$freq[binned_mld$freq != 0])/100
theoretical_plot <- theoretical_mlds_df %>%
  pivot_longer(c(mh, mc), names_to = "mcmh", values_to = "freq") %>%
  ggplot(aes(x = r, y = freq, color = mcmh, group = interaction(tau, mcmh))) +
  geom_line(alpha = 0.1) +
  geom_line(data = mc_all, aes(x = r, y = mc_all, group = NULL), linetype = "dashed", color = "#000066", linewidth = 2) +
  geom_line(data = mh_all, aes(x = r, y = mh_all, group = NULL), linetype = "dashed", color = "#660000", linewidth = 2) +
  geom_point(data = binned_mld, aes(x = r, y = freq, group = NULL), color = "black", alpha = 1, size = 2) +
  scale_y_log10(limits = c(lim_freq, NA)) +
  scale_x_log10() +
  annotation_logticks() +
  theme_classic() +
  scale_color_manual(values = c("mc" = "blue", "mh" = "red"))

ggsave(paste0(output_dir, species[1], "_vs_", species[2], "_fitted_single_mlds_plot.png"), theoretical_plot)


# MDS - inflexion --------------------------------------------------------------

comp_infl <- theoretical_mlds_df %>%
  mutate(diff_mhmc = mh - mc) %>%
  mutate(sup_mcmh = ifelse(diff_mhmc < 0, 1, 0)) %>%
  filter(sup_mcmh != 0) %>%
  group_by(comp) %>%
  summarise(r_inflexion = max(r)) %>%
  filter(r_inflexion > min_r_infl) %>%
  pull(comp)

infl_df <- binned_comparisons_df %>%
  select(comp, tau, rho) %>%
  mutate(inflexion = ifelse(comp %in% comp_infl, "yes", "no")) %>%
  separate_wider_delim(cols = comp, delim = "_vs_", names = species, cols_remove = FALSE)

write_csv(infl_df, paste0(output_dir, species[1], "_vs_", species[2], "_inflexion_single.csv"))

mds_full_time <- 0
if (get_comp_number(species[1], species[2], full_mlds_dir) >= 6) {
  mds_start_time <- Sys.time()
  counts_bins <- colnames(full_mld) %>% as.numeric()
  distance_comp_matrix <- matrix(0, nrow = nrow(full_mld), ncol = nrow(full_mld))
  colnames(distance_comp_matrix) <- rownames(distance_comp_matrix) <- rownames(full_mld)
  for (row in seq_len(nrow(distance_comp_matrix))) {
    for (col in seq(row, ncol(distance_comp_matrix))) {
      if (row != col) {
        distance_comp_matrix[row, col] <- ks.test(
          rep(full_mld[row, ], times = counts_bins),
          rep(full_mld[col, ], times = counts_bins)
        )$stat %>% unname()
      }
    }
  }

  mds_stop_time <- Sys.time()
  mds_full_time <- difftime(mds_stop_time, mds_start_time, units = "min")
  distance_comp_matrix <- distance_comp_matrix + t(distance_comp_matrix)
  write.csv(distance_comp_matrix, paste0(output_dir, species[1], "_vs_", species[2], "_ks_distancemat.csv"))
  distance_comp_matrix <- as.dist(distance_comp_matrix)
  mds_fit <- cmdscale(distance_comp_matrix, eig = TRUE, k = 2)
  mds_proj <- tibble(x = mds_fit$points[, 1], y = mds_fit$points[, 2], comp = rownames(mds_fit$points))


  mds_proj_infl <- mds_proj %>%
    mutate(infl = ifelse(comp %in% comp_infl, "inflexion", "pas d'inflexion"))


  mds_plot <- ggplot(mds_proj_infl, aes(x = x, y = y, color = infl)) +
    geom_point()


  ggsave(paste0(output_dir, species[1], "_vs_", species[2], "_mds_inflexion_plot.png"), mds_plot)

}

if (length(comp_infl) > 0) {
  # looking at individual mlds - with inflexion point
  # 1. plotting inflexion examples
  binned_df_infl <- binned_comparisons_df %>%
    filter(comp %in% comp_infl) %>%
    slice_sample(n = 20) %>% # not too many lines
    select(-tau, -rho) %>%
    as.data.frame() %>%
    column_to_rownames("comp") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("r") %>%
    mutate(r = as.numeric(r)) %>%
    pivot_longer(!r, values_to = "freq", names_to = "comp") %>%
    left_join(theoretical_mlds_df, by = c("r", "comp")) %>%
    pivot_longer(c(mc, mh), names_to = "mcmh", values_to = "estimations")


  # looking at individual mlds - without inflexion point
  n_wo_inflexion <- get_comp_number(species[1], species[2], full_mlds_dir) - length(comp_infl)
  n_inflexion_plotted <- min(20, length(comp_infl))
  if(n_wo_inflexion == 0) {
    plot_data_infl <- binned_df_infl %>%
      mutate(inf = "inflexion")
  } else {
    binned_df_samp <- binned_comparisons_df %>%
      filter(!comp %in% comp_infl) %>%
      slice_sample(n = min(n_wo_inflexion, n_inflexion_plotted)) %>%
      select(-tau, -rho) %>%
      as.data.frame() %>%
      column_to_rownames("comp") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("r") %>%
      mutate(r = as.numeric(r)) %>%
      pivot_longer(!r, values_to = "freq", names_to = "comp") %>%
      left_join(theoretical_mlds_df, by = c("r", "comp")) %>%
      pivot_longer(c(mc, mh), names_to = "mcmh", values_to = "estimations")

    plot_data_infl <- bind_rows(list("inflexion" = binned_df_infl, "no_inflexion_samp" = binned_df_samp), .id = "inf")
  }



  lim_freq <- min(plot_data_infl$freq[plot_data_infl$freq > 0]) / 10
  inflexion_plot <- plot_data_infl %>%
    ggplot(aes(x = r)) +
    geom_point(aes(y = freq, color = inf)) +
    geom_line(aes(y = estimations, group = interaction(comp, mcmh), color = inf, linetype = mcmh)) +
    scale_x_log10() +
    scale_y_log10(limits = c(lim_freq, NA)) +
    scale_color_manual(values = c("inflexion" = "darkolivegreen3", "no_inflexion_samp" = "coral3")) +
    scale_linetype_manual(values = c("mc" = "solid", "mh" = "dashed"))

  ggsave(paste0(output_dir, species[1], "_vs_", species[2], "_inflexion_mlds_plot.png"), inflexion_plot)

  # 2. Qui est dans l'inflexion ?
  entropy <- function(vec) {
    ptab <- prop.table(table(vec))
    -sum(ptab * log2(ptab))
  }

  genomes_infl_list <- binned_df_infl %>%
    pull(comp) %>%
    unique() %>%
    str_split("_vs_")

  genomes_infl_df <- map_df(genomes_infl_list, ~ data.frame("bac_1" = .x[1], "bac_2" = .x[2]))
  colnames(genomes_infl_df) <- species
  res_inflexion <- genomes_infl_df %>%
    pivot_longer(everything(), names_to = "Species", values_to = "Genomes") %>%
    group_by(Species) %>%
    summarise(n = n(), entropy = entropy(Genomes))
  res_proper <- bind_cols(
    res_inflexion %>%
    pivot_wider(names_from = Species, values_from = entropy) %>%
    select(-n),
    tibble(
      "n_comp_tot" = get_comp_number(species[1], species[2], full_mlds_dir),
      "n_comp_infl" = sum(res_inflexion$n),
      "mc_inf_0" = length(negative_fitted_mc),
      "fit_params_time" = fit_params_full_time,
      "mds_time" = mds_full_time
    )
  )
  write_csv(res_proper, paste0(output_dir, species[1], "_vs_", species[2], "_inflexion_res.csv"))
} else {
  res_proper <- tibble(
    "bl" = NA,
    "bla" = NA,
    get_comp_number(species[1], species[2], full_mlds_dir),
    0,
    length(negative_fitted_mc),
    fit_params_full_time,
    mds_full_time
  )
  colnames(res_proper) <- c(species, "n_comp_tot", "n_comp_infl", "mc_inf_0", "fit_params_time", "mds_time")
  write_csv(res_proper, paste0(output_dir, species[1], "_vs_", species[2], "_inflexion_res.csv"))
}
