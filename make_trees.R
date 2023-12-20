library(ggtree)
library(stringr)
library(dplyr)
library(tibble)
library(readr)
library(phangorn)
library(tidyr)
library(ggplot2)
library(janitor)

# Read data --------------------------------------------------------------------

params_dir <- snakemake@params["fitted_params_dir"]
species_list <- snakemake@params["species_lsit"]
metadata <- snakemake@params["metadata"]
results_dir <- snakemake@params["results_dir"]

fitted_params_files <- list.files(params_dir)

fitted_params <- tibble(
  "log_tau" = numeric(),
  "log_rho" = numeric(),
  "L0" = numeric(),
  "bac_1" = character(),
  "bac_2" = character()
)


for(spec_par in fitted_params_files) {
  species <- str_split_1(spec_par, "_vs_")[1:2]
  tmp_df <- read_csv(file.path(params_dir, spec_par))
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

pseudo_distance <- matrix(0, nrow = length(species_list), ncol = length(species_list))
colnames(pseudo_distance) <- rownames(pseudo_distance) <- species_list

for(row in rownames(pseudo_distance)) {
  for(col in colnames(pseudo_distance)) {
    if (col == row)
      next
    logtau <- fitted_params_intra %>%
      filter(
      grepl(col, bac_1) & grepl(row, bac_2) |
      grepl(col, bac_2) & grepl(row, bac_1)
      ) %>%
      pull(log_tau)
    pseudo_distance[row, col] <- logtau
  }
}

tree_upgma <- upgma(as.dist(10^(pseudo_distance)))

# distances and tau/distance comparison ----------------------------------------

coph_distances <- cophenetic(tree_upgma) %>%
  as.data.frame() %>%
  rownames_to_column("bac_1") %>%
  pivot_longer(!bac_1, names_to = "bac_2", values_to = "distance") %>%
  filter(bac_1 != bac_2)


distance_and_fitted <- fitted_params_intra %>%
  left_join(coph_distances) %>%
  mutate(tau = 10^log_tau) %>%
  mutate(relative_dif = abs(tau - distance)/(tau+distance))

ggplot(distance_and_fitted, aes(x = distance, y = tau)) +
  geom_point() +
  geom_function(fun = identity)

ggplot(distance_and_fitted, aes(x = relative_dif)) +
  geom_histogram(bins = 20, color = "darkblue", fill = "lightblue")


# Add external data ------------------------------------------------------------

if (snakemake@params["mock"] == "yes") {
  family_df <- read_csv(metadata) %>% dplyr::rename(label = Species)
  label_order <- tree_upgma %>%
    as_tibble %>%
    filter(!is.na(label)) %>%
    select(label)

  fam <- label_order %>%
    left_join(family_df) %>%
    column_to_rownames("label")

  p <- ggtree(tree_upgma) + geom_tiplab()
  p <- revts(p) + scale_x_continuous(labels = abs)

  gh <- gheatmap(p, fam,
                 colnames = FALSE,
                 legend_title = "Family",
                 width = 0.1,
                 offset = 1.2e8
                 ) +
    scale_x_ggtree() +
    theme_tree2(legend.position = "bottom",
                legend.box = "vertical", legend.margin = margin())

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
