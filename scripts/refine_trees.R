library(ggtree)
library(dplyr)
library(tibble)
library(readr)
library(phangorn)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(TreeDist)


get_distance <- function(res_df, column = "log10tau") {
  matrix_tau <- bind_rows(
    res_df |> select(all_of(c("species_1", "species_2", column))),
    res_df |>
    select(all_of(c("species_1", "species_2", column))) |>
    rename(species_1 = species_2, species_2 = species_1) |>
    relocate(species_1)
  ) |>
    pivot_wider(
      id_cols = species_1,
      names_from = species_2,
      values_from = all_of(column)
    ) |>
    column_to_rownames("species_1") |>
    as.matrix()
  species_names <- rownames(matrix_tau)
  matrix_tau <- matrix_tau[species_names, species_names]
  diag(matrix_tau) <- 0
  distance_tau <- as.dist(matrix_tau)
  return(distance_tau)
}

plot_tree <- function(res_df, taxon_df, tree_annotation, leaves_name, plot=TRUE, ret_dist=FALSE) {
  outgroup <- c(
    ## "Staphylococcus saprophyticus",
    ## "Staphylococcus epidermidis",
    ## "Staphylococcus haemolyticus",
    ## "Staphylococcus aureus"
    "Lactiplantibacillus paraplantarum"
  )
  staph <- c(
    "Staphylococcus saprophyticus",
    "Staphylococcus epidermidis",
    "Staphylococcus haemolyticus",
    "Staphylococcus aureus"
  )
  matrix_tau <- bind_rows(
    res_df |> select(species_1, species_2, log10tau),
    res_df |> select(species_1 = species_2, species_2 = species_1, log10tau)
  ) |>
    pivot_wider(
      id_cols = species_1,
      names_from = species_2,
      values_from = log10tau
    ) |>
    column_to_rownames("species_1") |>
    as.matrix()
  species_names <- rownames(matrix_tau)
  matrix_tau <- matrix_tau[species_names, species_names]
  diag(matrix_tau) <- 0
  distance_tau <- as.dist(matrix_tau)


  tree <- upgma(10^distance_tau)
  ## tree <- nj(10^distance_tau)
  ## tree <- fastme.bal(10^distance_tau)
  ## tree <- root(tree, outgroup)
  if (!plot && ret_dist)
    return(list(tree, distance_tau))
  if (!plot && !ret_dist)
    return(plot)

  coph_distances <- cophenetic(tree) %>%
    as.data.frame() %>%
    rownames_to_column("species_1") %>%
    pivot_longer(!species_1, names_to = "species_2", values_to = "distance") %>%
    filter(species_1 != species_2)

  distance_and_fitted <- res_df %>%
    filter(infl_exist == "yes") %>%
    inner_join(coph_distances) %>%
    mutate(tau = 10^log10tau) %>%
    rowwise() %>%
    mutate(staph = if_else(species_1 %in% staph && species_2 %in% staph, "both", if_else(species_1 %in% staph, "one", if_else(species_2 %in% staph, "one", "none")))) %>%
    mutate(ref_dist =  tau) %>%
    ungroup() %>%
    mutate(relative_dif = abs(ref_dist - distance) / (ref_dist + distance))

  difi <- ggplot(distance_and_fitted, aes(x = distance, y = ref_dist, color=staph)) +
    geom_point() +
    geom_function(fun = identity)
  difi <- difi + scale_x_log10() +
    scale_y_log10()


  family_df <- taxon_df %>%
    dplyr::rename(label = all_of(leaves_name))

  label_order <- tree %>%
    as_tibble() %>%
    filter(!is.na(label)) %>%
    select(label)

  fam <- family_df %>%
    distinct(.data[[tree_annotation]], label) %>%
    inner_join(label_order) %>%
    select(label, .data[[tree_annotation]]) %>%
    column_to_rownames("label")


  p <- ggtree(tree) + geom_tiplab()
  p <- revts(p) + scale_x_continuous(labels = abs)

  gh <- gheatmap(
    p, fam,
    colnames = FALSE,
    legend_title = tree_annotation,
    width = 0.1,
    offset = 2e8
  )
  gh <- gh +
    scale_x_ggtree() +
    theme_tree2(legend.position = "bottom",
                legend.box = "vertical", legend.margin = margin())
  list(gh, difi)
}

plot_reference <- function(tmp_taxon, taxonomy, gtdb_tree, plot=TRUE) {
  studied_species <- tmp_taxon |>
    pull("species.gtdb") |>
    unique()

  accessions_gtdb <- taxonomy |>
    filter(species %in% studied_species) |>
    pull(assembly_accession)

  species_accession <- taxonomy |>
    filter(species %in% studied_species) |>
    select(assembly_accession, species)

  pruned <- keep.tip(gtdb_tree, tip = accessions_gtdb)
  tibble_tree <- as_tibble(pruned) |>
    mutate(label = str_replace_all(label, "'", ""))

  ## genomes <- tibble_tree |>
  ##   select(label) |>
  ##   left_join(species_accession, by = c("label" = "assembly_accession")) |>
  ##   drop_na()

  annotated_tree <- tibble_tree |>
    left_join(species_accession, by = c("label" = "assembly_accession"))
  ## annotated_tree <- full_join(tibble_tree, genomes, by = "label")
  if (!plot)
    return(annotated_tree)


  heatmap_ann <- genomes |>
    column_to_rownames("label") |>
    select(family)

  p <- ggtree(tidytree::as.treedata(annotated_tree)) +
    geom_tiplab(aes(label=species)) +
    theme_tree()
  p
}

plot_delta <- function(old_res_df, new_res_df, plot = TRUE) {
  old_res_df <- old_res_df |>
    mutate(tau = 10^log10tau)
  fit_dist <- get_distance(new_res_df, "tau")
  old_dist <- get_distance(old_res_df, "tau")

  old_delta <- delta.plot(old_dist, plot=FALSE)
  old_delta_mean <- mean(old_delta$delta.bar)
  fit_delta <- delta.plot(fit_dist, plot=FALSE)
  fit_delta_mean <- mean(fit_delta$delta.bar)
  delta_hist_df <- tibble("delta"=seq(0, 1, length.out = 20), "pairwise"=old_delta$counts, "global"= fit_delta$counts) |>
    pivot_longer(!delta, names_to = "fitting_procedure", values_to = "counts")
  stopifnot(rownames(as.matrix(old_dist)) == rownames(as.matrix(fit_dist)))
  average_delta_df <- tibble("tip" = rownames(as.matrix(old_dist)), "pairwise" = old_delta$delta.bar, "global" = fit_delta$delta.bar) |>
    pivot_longer(!tip, names_to = "fitting_procedure", values_to = "delta_bar")

  if(plot == FALSE) {
    return(list(delta_hist_df, average_delta_df))
  }
  pl_hist <- ggplot(delta_hist_df, aes(x = delta, y=counts, color = fitting_procedure))+
    geom_freqpoly(stat="identity") +
    scale_color_manual(values=c("pairwise" = "#ff7f0e", "global" = "#1f77b4"),
                       labels = c("pairwise" = paste0("\u03b8 ratio from literature, \u03b4 = ", round(old_delta_mean, 2)),
                                  "global" = paste0("Fitted \u03b8, \u03b4 = ", round(fit_delta_mean, 2)))) +
    theme_bw() +
    theme(
      legend.position = c(0.6, 0.6),
      legend.background = element_rect(fill = "white"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10)
    )
  pl_hist
}

drop_species <- function(res_df, species_to_drop) {
  filter_sp <- regex(paste(species_to_drop, collapse = "|"))
  res_df |> filter(!(str_detect(species_1, filter_sp) | str_detect(species_2, filter_sp)))
}
drop_species_taxon  <- function(taxon_df, species_to_drop, leaves_name) {
  filter_sp <- regex(paste(species_to_drop, collapse = "|"))
  taxon_df |> filter(!str_detect(.data[[leaves_name]], filter_sp))
}

# inferring muc fitted trees ----------------------------------------------------------------------------------------------
# Bacillales
outgroups <- c(
  "Staphylococcus saprophyticus",
  "Staphylococcus epidermidis",
  "Staphylococcus haemolyticus",
  "Staphylococcus aureus",
  "Lactiplantibacillus paraplantarum",
  "Paenibacillus polymyxa",
  "Paenibacillus polymyxa_B",
  "Lysinibacillus fusiformis"
)
res_df <- read_csv("/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/results.csv")
taxon_df <- read_csv("/home/paulimer/Data/bacillales_species/species_taxon_w_staph_w_outg.csv")
tree_annotation <- "family.gtdb"
leaves_name <- "species.gtdb"
grid_tree <- plot_tree(res_df, taxon_df, tree_annotation, leaves_name, FALSE)
treeio::write.tree(as.phylo(grid_tree), "global_tree.nwk")
# comparing with not fitting muc
old_res_df <- read_csv("/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/results.csv")
old_res_df <- old_res_df |>
  mutate(tau = 10^log10tau)
## fit_dist <- get_distance(res_df, "tau")
## fit_species <- colnames(as.matrix(fit_dist))
## old_dist <- get_distance(old_res_df, "tau")
## old_delta <- delta.plot(old_dist, plot=FALSE)
## old_species <- colnames(as.matrix(old_dist))
## old_delta_mean <- mean(old_delta$delta.bar)
## fit_delta <- delta.plot(fit_dist, plot=FALSE)
## fit_delta_mean <- mean(fit_delta$delta.bar)

# ANI
ani_dist <- get_distance(old_res_df, "average_divergence")
delta_species <- colnames(as.matrix(ani_dist))
ani_delta <- delta.plot(ani_dist, plot = FALSE)
ani_delta_mean <- mean(ani_delta$delta.bar)

# plotting delta plot
delta_hist_df <- tibble("delta"=seq(0, 1, length.out = 20), "pairwise"=old_delta$counts, "global"= fit_delta$counts) |>
  pivot_longer(!delta, names_to = "fitting_procedure", values_to = "counts")
pl_hist_2 <- plot_delta(old_res_df, res_df)
delta_res <- plot_delta(old_res_df, res_df, plot=FALSE)
write_csv(delta_res[[1]], "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/pair_global_hist_delta.csv")
write_csv(delta_res[[2]], "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/pair_global_average_delta.csv")
ggsave("/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/pair_global_delta.svg", pl_hist_2, height = 3, width = 4, dpi = 300)


# plotting by species
old_df <- tibble("species"=old_species, "average_delta"= old_delta$delta.bar)
fit_df <- tibble("species"=fit_species, "average_delta"= fit_delta$delta.bar)
ani_df <- tibble("species"=delta_species, "average_delta"= ani_delta$delta.bar)
delta_df <- inner_join(old_df, fit_df, by="species", suffix = c("_old", "_new")) |>
  pivot_longer(!species, names_to = "old_new", values_to = "average_delta") |>
  mutate(old_new = str_extract(old_new, "...$"))
pl <- ggplot(delta_df, aes(x = species, y = average_delta, fill = old_new))+
  geom_col(position = "dodge2") +
  theme(axis.text.x = element_text(angle = 90))

progress_df <- inner_join(old_df, fit_df, by="species", suffix = c("_old", "_new")) |>
  mutate(progress = (average_delta_old - average_delta_new) / average_delta_old) |>
  mutate(outgroup = if_else(species %in% outgroups, TRUE, FALSE))
pl_g <- ggplot(progress_df, aes(x = species, y = progress, fill = outgroup))+
  geom_col(position = "dodge2") +
  theme(axis.text.x = element_text(angle = 90))

# comparing ANI distances delta


# Enterobacteriales
new_res_df_entero <- read_csv("/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/results.csv")
old_res_df_entero <- read_csv("/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/results.csv")
taxon_df_entero <- read_csv("/home/paulimer/Data/latest_entero_genomes_no_plasmid/misha_taxon_w_annotation_v3.csv")
new_res_df_entero <- new_res_df_entero |>
  mutate(log10tau = log10(tau))
tree_annotation <- "family.gtdb"
leaves_name <- "species.gtdb"
grid_tree <- plot_tree(new_res_df_entero, taxon_df_entero, tree_annotation, leaves_name, FALSE)
treeio::write.tree(as.phylo(grid_tree), "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/global_entero_tree.nwk")
delta_entero <- plot_delta(old_res_df_entero, new_res_df_entero)
delta_res_entero <- plot_delta(old_res_df_entero, new_res_df_entero, plot=FALSE)
write_csv(delta_res_entero[[1]], "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/pair_global_hist_delta.csv")
write_csv(delta_res_entero[[2]], "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/pair_global_average_delta.csv")
ggsave("/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/pair_global_delta.svg", delta_entero, height = 3, width = 4, dpi = 300)

# Methanobacteriales
new_res_df_methano <- read_csv("/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/results.csv")
old_res_df_methano <- read_csv("/home/paulimer/Documents/results_bacteria_mlds/results_methano/results.csv")
taxon_df_methano <- read_csv("/home/paulimer/Documents/Encadrement/Stage/M2/Stage_bacteries/Archaea/Methanobacteriota/genomes/taxon_Methanobacteriota.csv")
new_res_df_methano <-   new_res_df_methano |>
  mutate(log10tau = log10(tau))
tree_annotation <- "family.gtdb"
leaves_name <- "species.gtdb"
grid_tree <- plot_tree(new_res_df_methano, taxon_df_methano, tree_annotation, leaves_name, FALSE)
treeio::write.tree(as.phylo(grid_tree), "/home/paulimer/Documents/results_bacteria_mlds/results_methano/global_methano_tree.nwk")
delta_methano <- plot_delta(old_res_df_methano, new_res_df_methano)
delta_res_methano <- plot_delta(old_res_df_methano, new_res_df_methano, plot=FALSE)
write_csv(delta_res_methano[[1]], "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/pair_global_hist_delta.csv")
write_csv(delta_res_methano[[2]], "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/pair_global_average_delta.csv")
ggsave("/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/pair_global_delta.svg", delta_methano, height = 3, width = 4, dpi = 300)


# grouped family replot
res_df <- read_csv("/home/paulimer/Documents/results_bacteria_mlds/temp_res/results_bacillales_species_w_staph_family/results.csv")
taxon_df <- read_csv("/home/paulimer/Data/bacillales_species/species_taxon_w_staph.csv")
tree_annotation <- "family.gtdb"
leaves_name <- "family.gtdb"
tmp_gh <- plot_tree(res_df, taxon_df, "order.gtdb", "family.gtdb")[[1]]
ggsave(paste0("/home/paulimer/Documents/results_bacteria_mlds/temp_res/results_bacillales_species_w_staph_family/mosaic_tree_order.svg"), tmp_gh, width = 8.5, height = 6, dpi = 300)


# selecting genera for plotting only them ------------------------------------------------------------------------------------
res_dir <- "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/"
res_df <- read_csv(paste0(res_dir, "results.csv"))
taxon_df <- read_csv("/home/paulimer/Data/bacillales_species/species_taxon_w_staph_w_outg.csv")
tree_annotation <- "family.gtdb"
leaves_name <- "species.gtdb"
gtdb_tree <- ape::read.tree("/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_v5_lastz/bac120_r226.tree")
gtdb_tree_genomes <- gtdb_tree$tip.label
taxonomic_hierarchy <- c("domain", "phylum", "class", "order", "family", "genus", "species")
taxonomy <- read_tsv("/home/paulimer/Data/bacilla_genomes/bac120_taxonomy_r226.tsv.gz", col_names = FALSE) %>%
  rename(assembly_accession = X1, taxonomic_lineage = X2) %>%
  separate_wider_delim(taxonomic_lineage, delim = ";", names = taxonomic_hierarchy) %>%
  mutate(across(all_of(taxonomic_hierarchy), ~str_remove_all(., ".*__"))) |>
  filter(assembly_accession %in% gtdb_tree_genomes)


fams <- taxon_df |>
  select(family.gtdb) |>
  unique()

selected_fam <- "Bacillaceae"
selected_outgroup <- "Bacillus_A anthracis"
species_to_drop <- taxon_df |>
  filter(family.gtdb != selected_fam & species.gtdb != selected_outgroup) |>
  select(species.gtdb) |>
  unique() |>
  mutate(species.gtdb = paste0("^", species.gtdb, "$")) |>
  pull(species.gtdb)
## species_to_drop<- c("blablabla")

tmp_taxon <- drop_species_taxon(taxon_df, species_to_drop, "species.gtdb")
write_csv(tmp_taxon, "bacillaceae_taxon.csv")
tmp_res <- drop_species(res_df, species_to_drop)


reference <- plot_reference(tmp_taxon, taxonomy, gtdb_tree, plot = FALSE)
reference <- reference |>
  mutate(label = case_when(
    grepl("^RS_GCF|^GB_GCA", label) ~ species,
    grepl("^[0-9]", label) ~ str_extract(label, "^[0-9.]+"),
    .default = label
  ))
treeio::write.tree(as.phylo(reference), "ref_tree_bacillaceae_g.nwk")
mosaic <- plot_tree(tmp_res, tmp_taxon, "genus.gtdb", "species.gtdb", plot=FALSE)
treeio::write.tree(as.phylo(mosaic), "mosaic_tree_bacillaceae_g.nwk")

distance <- TreeDistance(as.phylo(reference), as.phylo(mosaic))


ggsave("mosaic_tree.png", mosaic[[1]], dpi=300, width = 8, height = 8)
ggsave("mosaic_diag.png", mosaic[[2]], dpi=300, width = 8, height = 8)
ggsave("reference_tree.png", p, dpi=300, width = 8, height = 4)




# subsampling part ---------------------------------------------------------------------
species <- unique(taxon_df |> select(all_of(leaves_name)))
staph <- c(
  "^Staphylococcus saprophyticus$",
  "^Staphylococcus epidermidis$",
  "^Staphylococcus haemolyticus$",
  "^Staphylococcus aureus$"
)
suspicious_outgroup <- c(
  "^Priestia megaterium_C$",
  "^Priestia megaterium_D$",
  "^Priestia megaterium$",
  "^Priestia zanthoxyli$",
  "^Priestia flexa$",
  "^Priestia aryabhattai$"
)
taxon_df <- drop_species_taxon(taxon_df, staph, "species.gtdb")
res_df <- drop_species(res_df, staph)

subsets <- list()
for (i in seq_len(length(suspicious_outgroup)))
  subsets[[i]] <-  combn(suspicious_outgroup, i, simplify = FALSE)
subsets <- list_flatten(subsets)

results_dir <- "/home/paulimer/Documents/results_bacteria_mlds/temp_res/results_bacillales_species_w_staph/"
for (subset in subsets) {
  tmp_res <- drop_species(res_df, subset)
  tmp_taxon <- drop_species_taxon(taxon_df, subset, "species.gtdb")
  tmp_gh <- plot_tree(tmp_res, tmp_taxon, "family.gtdb", "species.gtdb")[[1]]
  ggsave(paste0(results_dir, "suspicious_no",paste0(subset, collapse = "-"), "tree.svg"), tmp_gh, width = 8.5, height = 6, dpi = 300)
}

tmp_fig <- plot_tree(res_df, taxon_df, "family.gtdb", "species.gtdb")
ggsave(paste0(results_dir, "nj_rooted_tree.svg"), tmp_gh, width = 8.5, height = 6, dpi = 300)
