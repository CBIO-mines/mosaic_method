library(tidyverse)
library(janitor)
library(jsonlite)
library(RColorBrewer)
library(viridis)
# read data
ncbi_meta <- fromJSON("../ncbi_metadata/dataset_summary.json", flatten = TRUE)$reports
tetramers <- read_csv("./tetra_stats_df_all_recode.csv")

# preprocess
tetramers_ncbi <- tetramers %>%
  filter(grepl("GCF", fasta)) %>%
  rowwise() %>% mutate(accession = substr(fasta, 1, 15)) %>% ungroup() %>%
  left_join(ncbi_meta, by = "accession")
tetramers_img <- tetramers %>%
  filter(!grepl("GCF", fasta)) %>%
  separate_wider_delim(fasta, ".", names = c("ncbi_biosample_accession", "ext"), cols_remove = FALSE) %>%
  select(-ext)
tetramers <- bind_rows(tetramers_ncbi, tetramers_img)
## filter_out_weird_esch <- c("SAMN02949651.fasta", "SAMN02949652.fasta", "SAMN17703722.fasta", "SAMN07305427.fasta", "GCF_022049065.2_ASM2204906v2_genomic.fna")
## tetramers <- tetramers %>% filter(!fasta %in% filter_out_weird_esch)

# calulate PCA
pca_cols <- c(colnames(tetramers)[str_detect(colnames(tetramers), "[ATGC][ATGC][ATGC][ATGC]")], "fasta")
mat <- tetramers %>% select(all_of(pca_cols)) %>%
  as.data.frame() %>%
  column_to_rownames("fasta")
pc <- prcomp(mat)
eig <- (pc$sdev)^2
variance <- eig * 100 / sum(eig)
PCAdata <- as.data.frame(pc$x) %>%
  select(PC1, PC2, PC3) %>%
  rownames_to_column("fasta") %>%
  inner_join(tetramers, by = "fasta")

pca_data <- list("data" = PCAdata, "variance" = variance)

# try to PCA species and not genera
## list_ext_files <- list.files("../external_data/")
## all_ext_data <- map_df(list_ext_files, ~read_tsv(paste0("../external_data/", .x)))
# goal is to do a plot species wide and not genera wide
# let's start with bacillus
cereus_clade <- c(
  "Bacillus cereus",
  "Bacillus cytotoxicus",
  "Bacillus paramycoides",
  "Bacillus thuringiensis",
  "Bacillus luti",
  "Bacillus albus",
  "Bacillus mobilis",
  "Bacillus wiedmannii",
  "Bacillus tropicus",
  "Bacillus anthracis",
  "Bacillus paranthracis",
  "Bacillus pacificus",
  "Bacillus toyonensis",
  "Bacillus proteolyticus",
  "Bacillus nitratireducens",
  "Bacillus mycoides",
  "Bacillus pseudomycoides"
)
subtilis_clade <- c(
  "Bacillus subtilis",
  "Bacillus gibsonii",
  "Bacillus vallismortis",
  "Bacillus cabrialesii",
  "Bacillus tequilensis",
  "Bacillus halotolerans",
  "Bacillus mojavensis",
  "Bacillus atrophaeus",
  "Bacillus siamensis",
  "Bacillus velezensis",
  "Bacillus amyloliquefaciens",
  "Bacillus nakamurai",
  "Bacillus xiamenensis",
  "Bacillus stratosphericus",
  "Bacillus celllulasensis",
  "Bacillus aerophilus",
  "Bacillus altitudinis",
  "Bacillus australimaris",
  "Bacillus safensis",
  "Bacillus pumilus",
  "Bacillus zhangzhouensis",
  "Bacillus licheniformis",
  "Bacillus haynesii",
  "Bacillus aerius",
  "Bacillus paralicheniformis",
  "Bacillus swezeyi",
  "Bacillus sonorensis",
  "Bacillus glycinifermentans",
  "Bacillus gobiensis"
)

species_count_min <- pca_data$data %>% filter(species == "Bacillus") %>%
  group_by(average_nucleotide_identity.submitted_species) %>%
  summarise(n = n()) %>%
  filter(n>15) %>%
  pull(average_nucleotide_identity.submitted_species)
pca_data_bacillus <- pca_data$data %>%
  mutate(bacillus_or_no = ifelse(species == "Bacillus", "yes", "no")) %>%
  mutate(bacillaceae_or_no = ifelse(species %in% c("EscherichiaColi", "Lactiplantibacillus", "Staphylococcus"), "no", "yes")) %>%
  mutate(cereus_or_no = ifelse(average_nucleotide_identity.submitted_species %in% cereus_clade, "yes", "no")) %>%
  mutate(subtilis_or_no = ifelse(average_nucleotide_identity.submitted_species %in% subtilis_clade, "yes", "no")) %>%
  mutate(
    Species_bac = ifelse(cereus_or_no == "yes", "Cereus_clade",
      ifelse(subtilis_or_no == "yes", "Subtilis_clade",
        ## ifelse(average_nucleotide_identity.submitted_species %in% species_count_min, average_nucleotide_identity.submitted_species,
          ifelse(bacillus_or_no == "yes", "Other_Bacillus",
            ifelse(bacillaceae_or_no == "yes", "Other_Bacillaceae", "Other_Bacteria")
          )
        )
      )
    ## )
  ) %>%
  group_by(Species_bac)
# color settings

manual_levels <- c("Other_Bacillus", "Other_Bacillaceae", "Other_Bacteria")
all_bacs <- sort(unique(pca_data_bacillus$Species_bac))
colors_by_bact <- setNames(
  c(c("#DB4537", "gray9", "grey"), viridis_pal()(length(all_bacs) - length(manual_levels))),
  c(manual_levels, setdiff(all_bacs, manual_levels))
)

plot_res <- pca_data_bacillus %>%
  ggplot(
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(aes(color = Species_bac), alpha = 0.5) +
  ylab(paste0("PC2: ", round(pca_data$variance[2], 1), "% variance")) +
  xlab(paste0("PC1: ", round(pca_data$variance[1], 1), "% variance")) #+
  ## coord_fixed() #+
  ## scale_color_manual(values = colors_by_bact)

ggsave(plot = plot_res, filename = "./PCA_bacillaceae_pres.png", width = 8, height = 4, dpi = 300)
## ggsave(plot = plot_res, filename = "./PCA_subsamp_wo_plasmids_ggplot.png")

# plot assembly_stats.gc_percent
plot_gc <- ggplot(
  pca_data_bacillus,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(aes(color = assembly_stats.gc_percent), alpha = 0.5) +
  ylab(paste0("PC2: ", round(pca_data$variance[2], 1), "% variance")) +
  xlab(paste0("PC1: ", round(pca_data$variance[1], 1), "% variance")) +
  coord_fixed() +
  scale_color_viridis()
ggsave(plot = plot_gc, filename = "./PCA_subsamp_wo_plasmids_gc_colors.png")
# kmeans clustering

kinput <- pca_data_bacillus %>%
  select(PC1, PC2)
kres <- tibble(n_clusts = (length(all_bacs) - 10):length(all_bacs)) %>%
  mutate(
           kclust = map(
                             n_clusts,
                             ~kmeans(kinput, centers = .x, iter.max = 20, nstart = 5)
                           ),
           augmented = map(kclust, broom::augment, pca_data_bacillus %>% select(PC1, PC2, fasta)),
           tidied = map(kclust, broom::tidy),
           glanced = map(kclust, broom::glance)
         ) %>%
  select(-kclust)

point_assignments <- kres %>%
  select(n_clusts, augmented) %>%
  tidyr::unnest(augmented)
cluster_info <- kres %>%
  select(n_clusts, tidied) %>%
  tidyr::unnest(tidied)
model_stats <- kres %>%
  select(n_clusts, glanced) %>%
  tidyr::unnest(glanced)


ggplot() +
  geom_point(
    data = point_assignments, aes(x = long, y = lat, color = .cluster)
  ) +
  geom_point(
    data = cluster_info, aes(x = long, y = lat), size = 4, shape = "x"
  ) +
  facet_wrap(~n_clusts)

# Elbow chart
ggplot(data = model_stats, aes(n_clusts, tot.withinss)) +
  geom_line() +
  scale_x_continuous(limits = c(1, 12), breaks = seq(1, 12, 1)) +
  ggtitle("Total within sum of squares, by # clusters")

# LEs tailles
ggplot(
  pca_data$data,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(aes(color = assembly_stats.total_sequence_length)) +
  ylab(paste0("PC2: ", round(pca_data$variance[2], 1), "% variance")) +
  xlab(paste0("PC1: ", round(pca_data$variance[1], 1), "% variance")) +
  coord_fixed()



# read external_data
lacti_meta <- read_tsv("../external_data/Lactiplantibacillus.csv") %>% clean_names()

all_lacti_data <- PCAdata %>%
  filter(species == "Lactiplantibacillus") %>%
  select(-species) %>%
  separate_wider_delim(fasta, ".", names = c("ncbi_biosample_accession", "ext")) %>%
  select(-ext) %>%
  left_join(lacti_meta, by = "ncbi_biosample_accession") %>%
  mutate(gc_per = gc_count_assembled/genome_size_assembled)

all_lacti_data  %>% select(species, PC1, gc_per, horizontally_transferred_percent, gene_count_assembled)
all_lacti_data %>% slice_min(PC1) %>% select(ncbi_assembly_accession)



# plot cropping
plots_files <- list.files(".", "*.png")
walk(plots_files, knitr::plot_crop)

# plasmid plot

# number of species/genomes
hist_genomes <- pca_data_bacillus %>%
  filter(Species_bac != "Other_Bacteria") %>%
  group_by(Species_bac) %>%
  summarize(n = n())
hist_genomes$genera <- c("Bacillus", "Other_Bacillaceae", "Bacillus", "Bacillus")
ggplot(hist_genomes, aes(x = Species_bac, y = n, fill = genera)) +
  geom_bar(stat = "identity") +
  ylab("Genome count") +
  xlab("Species group")
