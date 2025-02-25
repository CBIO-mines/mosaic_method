library(ggtree)
library(dplyr)
library(tibble)
library(readr)
library(phangorn)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)

bacterias_in_reference <- function(reference, bacterias, column_pref="bac") {
  # finds if a comparisons exists in a reference of comparison
  columns_clust <- paste0(column_pref, "_",  c(1, 2))
  if (length(bacterias) == 2) {
    bacterias_ordered <- str_sort(bacterias, locale = "C")
    cond <- reference[[columns_clust[1]]] == bacterias_ordered[1] &
      reference[[columns_clust[2]]] == bacterias_ordered[2]
  } else if (length(bacterias == 1)) {
    cond <- reference[[columns_clust[1]]] == bacterias |
      reference[[columns_clust[2]]] == bacterias
  }
  return(cond)
}
# Read data --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

res_df <- read_csv(args[1])
taxon_df <- read_csv(args[2])
cluster_name <- args[3]
species_list <- taxon_df %>% pull(.data[[cluster_name]]) %>% unique
results_dir <- args[4]
length_dir <- args[5]
tree_annotation <- args[6]
min_r_infl <- as.numeric(args[7])

if (results_dir[length(results_dir)] != "/") {
  results_dir <- paste0(results_dir, "/")
}
if (length_dir[length(length_dir)] != "/") {
  length_dir <- paste0(length_dir, "/")
}

no_inflexion_comps <- res_df %>%
  filter(infl_exist == "no") %>%
  select(species_1, species_2)



# construct tree ---------------------------------------------------------------

pseudo_distance <- matrix(NA, nrow = length(species_list), ncol = length(species_list))
colnames(pseudo_distance) <- rownames(pseudo_distance) <- species_list

for (row_i in seq_len(nrow(pseudo_distance))) {
  for (col_i in seq(row_i, ncol(pseudo_distance))) {
    if (col_i == row_i)
      next
    if (no_inflexion_comps %>%
        bacterias_in_reference(., c(species_list[row_i], species_list[col_i]), "species") %>%
        any()){
      next
    }
    logtau <- res_df %>%
      filter(
        (bacterias_in_reference(., c(species_list[row_i], species_list[col_i]), "species"))
      ) %>%
      pull(log10tau)
    pseudo_distance[row_i, col_i] <- logtau
  }
}
# percentage of no_inflexion out of all necessary taus
missing_taus <- 2 * sum(is.na(as.dist(t(pseudo_distance))))/(length(species_list)*(length(species_list)-1))

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
# to get the proper pseudo distance
tau_distance <- as.dist(t(pseudo_distance))
stopifnot(sum(is.na(tau_distance)) == 0)

# checking treelikeness
delta_res <- delta.plot(10^tau_distance, plot = FALSE)
mean_delta <- mean(delta_res$delta.bar)
write_csv(tibble("mean_delta" = mean_delta, "missing_taus" = missing_taus), paste0(results_dir, "tree_stats.csv"))

png(paste0(results_dir, "delta_plot_hist.png"), width = 10, height = 8, units = "in", res = 300)
delta.plot(10^tau_distance, which = 1)
dev.off()


tree_upgma <- upgma(10^(tau_distance))

# distances and tau/distance comparison ----------------------------------------

coph_distances <- cophenetic(tree_upgma) %>%
  as.data.frame() %>%
  rownames_to_column("species_1") %>%
  pivot_longer(!species_1, names_to = "species_2", values_to = "distance") %>%
  filter(species_1 != species_2)

distance_and_fitted <- res_df %>%
  filter(infl_exist == "yes") %>%
  inner_join(coph_distances) %>%
  mutate(tau = 10^log10tau) %>%
  mutate(relative_dif = abs(tau - distance)/(tau+distance))

difi <- ggplot(distance_and_fitted, aes(x = distance, y = tau)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_function(fun = identity)

ggsave(paste0(results_dir, "fitteddistance_vs_founddistance.png"), difi)

difi_hist <- ggplot(distance_and_fitted, aes(x = relative_dif)) +
  geom_histogram(bins = 20, color = "darkblue", fill = "lightblue")

ggsave(paste0(results_dir, "hist_fitteddistance.png"), difi_hist)


# Add taxon, and counts/inflexion info ------------------------------------------

family_df <- taxon_df %>%
  dplyr::rename(label = all_of(cluster_name))

counts_df <- family_df %>%
  select(label, genome) %>%
  group_by(label) %>%
  summarise(count = n())

label_order <- tree_upgma %>%
  as_tibble %>%
  filter(!is.na(label)) %>%
  select(label)

fam <- family_df %>%
  distinct(.data[[tree_annotation]], label) %>%
  inner_join(label_order) %>%
  select(label, .data[[tree_annotation]]) %>%
  column_to_rownames("label")


p <- ggtree(tree_upgma) + geom_tiplab()
p <- revts(p) + scale_x_continuous(labels = abs)
    ## scale_x_continuous(labels=function(x) scales::comma(abs(x))) # <-- what you need is actually a function.
  ##

gh <- gheatmap(p, fam,
               colnames = FALSE,
               legend_title = tree_annotation,
               width = 0.1,
               offset = 1e8
               ) +
  scale_x_ggtree() +
  theme_tree2(legend.position = "bottom",
              legend.box = "vertical", legend.margin = margin())

## gh <- gh +
##   geom_facet(panel = "Genome count",
##              data = counts_df,
##              geom = geom_col,
##              aes(x = count),#, fill = Family),
##              orientation = "y"
##              )
## gh <- facet_widths(gh, widths = c(4, 1))
##   ## theme_tree2(legend.position=c(.05, .85))

## # according to ggtree doc FAQ
## gh <- gh + xlim_tree(0) + xlim_expand(c(0, 1000), "Genome count")

## d <- data.frame(.panel = c("Tree", "Genome count"),
##                 lab = c("tau/2", "count"),
##                 x = c(-1.5e8,100), y = -2)

## ghf <- gh + geom_text(aes(label=lab), data=d) +
##   coord_cartesian(clip='off') # +
##   ## theme(plot.margin=margin(6, 6, 40, 6))

## gh <- gh +
##   geom_facet(panel = "Inflexion percentage",
##              data = inflexions_per,
##              geom = geom_col,
##              aes(x = per_infl),#, fill = Family),
##              orientation = "y",
##              scales = "freex")


ggsave(paste0(results_dir, tree_annotation, "_tree_big.svg"), gh, width = 8.5, height = 6, dpi = 300)
