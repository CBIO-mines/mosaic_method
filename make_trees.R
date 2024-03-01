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
  # if manual execution, the paths and parameters need to be adapted

  snakemake <- Snakemake(
    input = list(
      fitted_params = list.files("./results_refseq_test/fitted_params/"),
      inflexion_file = "results_refseq_test/inflexion_exists.csv",
      inflexion_percentage = "results_refseq_test/inflexion_by_cluster.csv"
    ),
    output = list(),
    params = list(
      taxon_csv = "bacillaceae_taxon_test.csv",
      cluster_name = "genus",
      fitted_params_dir = "results_refseq_test/fitted_params/",
      results_dir = "results_refseq_test/",
      genome_wise_inflexion = "no",
      genome_lengths = "results_refseq_test/lengths_distributions/"
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

library(ggtree)
library(stringr)
library(dplyr)
library(tibble)
library(readr)
library(phangorn)
library(tidyr)
library(ggplot2)
library(janitor)
library(purrr)

read_counts <- function(x) {
  x_file <- paste0(x, "_distribution.csv")
  tmp_df <- read_csv(paste0(length_dir, x_file), col_types = cols(col_character(), col_double()))
  tibble("species" = x, "count" = nrow(tmp_df))
}

bacterias_in_reference <- function(reference, bacterias, column_pref="bac") {
  # finds if a comparisons exists in a reference of comparison
  columns_clust <- paste0(column_pref, "_",  c(1, 2))
  if (length(bacterias) == 2) {
    bacterias_ordered <- sort(bacterias)
    cond <- reference[[columns_clust[1]]] == bacterias_ordered[1] &
      reference[[columns_clust[2]]] == bacterias_ordered[2]
  } else if (length(bacterias == 1)) {
    cond <- reference[[columns_clust[1]]] == bacterias |
      reference[[columns_clust[2]]] == bacterias
  }
  return(cond)
}


# Read data --------------------------------------------------------------------

params_dir <- snakemake@params[["fitted_params_dir"]]
taxon_csv <- snakemake@params[["taxon_csv"]]
cluster_name <- snakemake@params[["cluster_name"]]
taxon_df <- read_csv(taxon_csv)
species_list <- taxon_df %>% pull(.data[[cluster_name]]) %>% unique
results_dir <- snakemake@params[["results_dir"]]
length_dir <- snakemake@params[["genome_lengths"]]
inflexions <- read_csv(snakemake@input[["inflexion_file"]])

fitted_params_files <- list.files(params_dir)

fitted_params <- tibble(
  "log_tau" = numeric(),
  "log_rho" = numeric(),
  "L0" = numeric(),
  "bac_1" = character(),
  "bac_2" = character()
)


for(spec_par in fitted_params_files) {
  level_1 <- str_split_1(spec_par, "_vs_")[1]
  level_2 <- str_sub(str_split_1(spec_par, "_vs_")[2], 1, -19)
  tmp_df <- read_csv(file.path(params_dir, spec_par), col_types = cols(.default = col_double()))
  tmp_df <- bind_cols(tmp_df, "bac_1" = level_1, "bac_2" = level_2)
  colnames(tmp_df) <- colnames(fitted_params)
  # réfléchir à dupliquer bac 1 et 2 pour avoir les paires possibles
  # genre bind rows aussi en inversant les bac
  fitted_params <- bind_rows(fitted_params, tmp_df)
}

fitted_params_intra <- fitted_params %>%
  filter(
    grepl(paste(species_list, collapse = "|"), bac_1) &
    grepl(paste(species_list, collapse = "|"), bac_2)
  )

th_comparison_df <- t(combn(species_list, 2, simplify = TRUE))

for (i in seq_len(nrow(th_comparison_df))) {
  select_row <- fitted_params_intra %>%
    filter(
      grepl(th_comparison_df[i, 1], bac_1) & grepl(th_comparison_df[i, 2], bac_2) |
      grepl(th_comparison_df[i, 1], bac_2) & grepl(th_comparison_df[i, 2], bac_1)
    )
  if(nrow(select_row) == 0)
    print(paste("Missing comparison",
                th_comparison_df[i, 1],
                "vs",
                th_comparison_df[i, 2]
                )
          )
}

# construct tree ---------------------------------------------------------------

pseudo_distance <- matrix(NA, nrow = length(species_list), ncol = length(species_list))
colnames(pseudo_distance) <- rownames(pseudo_distance) <- species_list

for (row_i in seq_len(nrow(pseudo_distance))) {
  for (col_i in seq(row_i, ncol(pseudo_distance))) {
    if (col_i == row_i)
      next
    if (inflexions[bacterias_in_reference(inflexions, c(species_list[row_i], species_list[col_i]), "cluster"), "infl_exist"] == "no")
      next

    logtau <- fitted_params_intra %>%
      filter(
        (bacterias_in_reference(., c(species_list[row_i], species_list[col_i])))
      ) %>%
      pull(log_tau)
    pseudo_distance[row_i, col_i] <- logtau
  }
}
# removing bacteria without a single comp with inflexion
# TODO still relevant ?
## empty_bacs <- c()
## for (row_i in seq_len(nrow(pseudo_distance))) {
##   if(all(is.na(pseudo_distance[row_i, ])))
##     empty_bacs <- c(empty_bacs, row_i)
## }
## pseudo_distance <- pseudo_distance[-empty_bacs, -empty_bacs]
# filling empty cells with the mean distance over the tree
mean_pseudo_distance <- mean(pseudo_distance, na.rm = TRUE)
for (row_i in seq_len(nrow(pseudo_distance))) {
  for (col_i in seq(row_i, ncol(pseudo_distance))) {
    if(row_i == col_i)
      next
    if (is.na(pseudo_distance[row_i, col_i])) {
      pseudo_distance[row_i, col_i] <- mean_pseudo_distance
    }
  }
}
pseudo_distance[is.na(pseudo_distance) & !is.nan(pseudo_distance)] <- 0
pseudo_distance <- pseudo_distance + t(pseudo_distance)

if (snakemake@params[["genome_wise_fit"]] == "yes") {
  # is there something to do additionnaly in this case ?
}
tree_upgma <- upgma(as.dist(10^(pseudo_distance)))

# distances and tau/distance comparison ----------------------------------------

coph_distances <- cophenetic(tree_upgma) %>%
  as.data.frame() %>%
  rownames_to_column("bac_1") %>%
  pivot_longer(!bac_1, names_to = "bac_2", values_to = "distance") %>%
  filter(bac_1 != bac_2)


distance_and_fitted <- fitted_params_intra %>%
  inner_join(coph_distances) %>%
  mutate(tau = 10^log_tau) %>%
  mutate(relative_dif = abs(tau - distance)/(tau+distance))

difi <- ggplot(distance_and_fitted, aes(x = distance, y = tau)) +
  geom_point() +
  geom_function(fun = identity)

ggsave(paste0(results_dir, "fitteddistance_vs_founddistance.png"), difi)

difi_hist <- ggplot(distance_and_fitted, aes(x = relative_dif)) +
  geom_histogram(bins = 20, color = "darkblue", fill = "lightblue")

ggsave(paste0(results_dir, "hist_fitteddistance.png"), difi_hist)


# Add taxon, inflexion percentages and counts/inflexion info ------------------------------------------

family_df <- taxon_df %>%
  dplyr::rename(label = all_of(cluster_name))

counts_df <- family_df %>%
  select(label, genome) %>%
  group_by(label) %>%
  summarise(count = n())

inflexions_per <- read_csv(snakemake@input[["inflexion_percentage"]]) %>%
  rename(label = all_of(cluster_name))

label_order <- tree_upgma %>%
  as_tibble %>%
  filter(!is.na(label)) %>%
  select(label)

fam <- family_df %>%
  distinct(family, label) %>%
  inner_join(label_order) %>%
  select(label, family) %>%
  column_to_rownames("label")


p <- ggtree(tree_upgma) + geom_tiplab()
p <- revts(p) + scale_x_continuous(labels = abs)
    ## scale_x_continuous(labels=function(x) scales::comma(abs(x))) # <-- what you need is actually a function.
  ##

gh <- gheatmap(p, fam,
               colnames = FALSE,
               legend_title = "Family",
               width = 0.1,
               offset = 1.2e8
               ) +
  scale_x_ggtree() +
  theme_tree2(legend.position = "bottom",
              legend.box = "vertical", legend.margin = margin())
gh <- gh +
  geom_facet(panel = "Genome count",
             data = counts_df,
             geom = geom_col,
             aes(x = count),#, fill = Family),
             orientation = "y",
             scales = "freex")
  ## theme_tree2(legend.position=c(.05, .85))

gh <- gh +
  geom_facet(panel = "Inflexion percentage",
             data = inflexions_per,
             geom = geom_col,
             aes(x = per_infl),#, fill = Family),
             orientation = "y",
             scales = "freex")


ggsave(paste0(results_dir, "family_tree.svg"), gh)
