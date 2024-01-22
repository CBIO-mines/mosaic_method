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

  species_l <- readLines("species_list_bacillaceae_all.txt")
  species_l <- gsub("\n", "", species_l)

  snakemake <- Snakemake(
    input = list(
      fitted_params = list.files("./results_refseq_real/fitted_params/")
    ),
    output = list(),
    params = list(
      metadata = "mock_metadata.csv",
      species_list = species_l,
      fitted_params_dir = "results_refseq_real/fitted_params/",
      mock = "yes",
      results_dir = "results_refseq_real/",
      inflexion_file = "inflexion_exists.csv",
      use_inflexion = "yes",
      genomes_lengths = "results_refseq_real/lengths_distributions/"
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

read_counts <- function(x) {
  x_file <- paste0(x, "_distribution.csv")
  tmp_df <- read_csv(paste0(length_dir, x_file), col_types = cols(col_character(), col_double()))
  tibble("species" = x, "count" = nrow(tmp_df))
}

# Read data --------------------------------------------------------------------

params_dir <- snakemake@params[["fitted_params_dir"]]
species_list <- snakemake@params[["species_list"]]
metadata <- snakemake@params[["metadata"]]
results_dir <- snakemake@params[["results_dir"]]
length_dir <- snakemake@params[["genomes_lengths"]]
if(snakemake@params[["use_inflexion"]] == "yes") {
  inflexions <- read_csv(paste0(results_dir, snakemake@params[["inflexion_file"]])) %>%
    separate_wider_delim(cols = comp, delim = "_vs_", names = c("bac_1", "bac_2"), cols_remove = FALSE)
}

fitted_params_files <- list.files(params_dir)

fitted_params <- tibble(
  "log_tau" = numeric(),
  "log_rho" = numeric(),
  "L0" = numeric(),
  "bac_1" = character(),
  "bac_2" = character()
)


for(spec_par in fitted_params_files) {
  species <- str_split_1(spec_par, "_")[c(1,3)]
  tmp_df <- read_csv(file.path(params_dir, spec_par), col_types = cols(.default = col_double()))
  tmp_df <- bind_cols(tmp_df, "bac_1" = species[1], "bac_2" = species[2])
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

bacterias_in_reference <- function(reference, bacterias) {
  # finds if a comparisons exists in a reference of comparison
  bacterias_ordered <- sort(bacterias)
  cond <- reference$bac_1 == bacterias_ordered[1] &
    reference$bac_2 == bacterias_ordered[2]
  return(cond)
}

if (snakemake@params[["use_inflexion"]] != "yes") {
  for (row in rownames(pseudo_distance)) {
    for (col in colnames(pseudo_distance)) {
      if (col == row) {
        next
      }
      logtau <- fitted_params_intra %>%
        filter(
          col == bac_1 & row == bac_2 |
            col == bac_2 & row == bac_1
        ) %>%
        pull(log_tau)
      pseudo_distance[row, col] <- logtau
    }
  }
} else {
  # filtering for inflexion
  for (row_i in seq_len(nrow(pseudo_distance))) {
    for (col_i in seq(row_i , ncol(pseudo_distance))) {
      if(row_i == col_i)
        next
      ## if(species_list[col_i] == "Bacillus")
      ##   browser()
      if (all(!bacterias_in_reference(inflexions, c(species_list[row_i], species_list[col_i]))))
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
  empty_bacs <- c()
  for (col_i in seq_len(ncol(pseudo_distance))) {
    if(all(is.na(pseudo_distance[, col_i])))
      empty_bacs <- c(empty_bacs, col_i)
  }
  # filling empty cells with means
  pseudo_distance <- pseudo_distance[-empty_bacs, -empty_bacs]
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


# Add external data ------------------------------------------------------------

if (snakemake@params["mock"] == "yes") {
  family_df <- read_csv(metadata) %>% dplyr::rename(label = Species)
  counts_df <- map_df(family_df$label, ~ read_counts(.x))
  label_order <- tree_upgma %>%
    as_tibble %>%
    filter(!is.na(label)) %>%
    select(label)

  fam <- label_order %>%
    inner_join(family_df) %>%
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
  gh +
    geom_facet(panel = "Genome count",
               data = counts_df,
               geom = geom_col,
               aes(x = count),#, fill = Family),
               orientation = "y",
               scales = "freex")
    ## theme_tree2(legend.position=c(.05, .85))

  ggsave(paste0(results_dir, "family_tree.svg"), gh)

} else {
  get_taxon_name <- function(path_external, comp_df, taxon_level) {
    tryCatch({
      pattern_files <- paste(union(comp_df$bac_1, comp_df$bac_2), collapse = "|")
      res_df <- tibble("label" = character(), "{taxon_level}" := character())
      for (file in list.files(path_external, pattern = pattern_files)) {
        file_read <- read_tsv(
          file = file.path(path_external, file),
          n_max = 1
        ) %>%
          select(all_of(c(taxon_level)))
        res_df <- bind_rows(
          res_df,
          bind_cols(file_read, tibble("label" = str_split_1(file, ".csv")[1]))
        )
      }
    },
    error = function(e) {
      print(paste("An error occured with file", file))
      print(e)
    }
    )
    return(res_df)
  }

  for (taxon_level in c("Phylum", "Class", "Order", "Family")) {
    label_order <- tree_upgma %>%
      as_tibble %>%
      filter(!is.na(label)) %>%
      select(label)

    external_taxon <- get_taxon_name(metadata, distance_and_fitted, taxon_level)

    fam <- label_order %>%
      left_join(external_taxon) %>%
      column_to_rownames("label")

    p <- ggtree(tree_upgma) + geom_tiplab()
    p <- revts(p) + scale_x_continuous(labels = abs)

    gh <- gheatmap(p, fam,
                   colnames = FALSE,
                   legend_title = taxon_level,
                   width = 0.1,
                   offset = 1.2e8
                   ) +
      scale_x_ggtree() +
      theme_tree2(legend.position = "bottom",
                  legend.box = "vertical", legend.margin = margin())

    ggsave(paste0("nice_tree_", taxon_level, ".svg"), gh)
  }

}
