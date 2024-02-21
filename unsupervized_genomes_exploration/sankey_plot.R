library(tidyverse)
library(janitor)
library(jsonlite)
library(RColorBrewer)
library(viridis)
library(ggsankey)
# read data
ncbi_meta <- fromJSON("../ncbi_metadata/dataset_summary.json", flatten = TRUE)$reports
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


species_columns <- c(
  "assembly_info.biosample.description.organism.organism_name",
  "ncbi_meta$assembly_info.biosample.description.organism.tax_id",
  "organism.tax_id",
  "organism.organism_name",
  "checkm_info.checkm_species_tax_id"
)
nodes_colnames <- c(
  "tax_id",
  "parent_tax_id",
  "rank",
  "embl_code",
  "division_id",
  "inherited_div_flag",
  "genetic_code_id",
  "inherited_GC_flag",
  "mitochondrial_genetic_code_id",
  "inherited_MGC_flag",
  "GenBank_hidden_flag",
  "hidden_subtree_root_flag",
  "comments",
  "plastid_genetic_code_id",
  "inherited_PGC_flag",
  "specified_species",
  "hydrogenosome_genetic_code_id",
  "inherited_HGC_flag"
)


nodes_df <- read_delim("nodes.dmp", delim = "\t|\t", col_names = nodes_colnames)
names_df <- read_delim("names.dmp", delim = "\t|\t", col_names = c("tax_id", "name_txt", "unique_name", "name_class"))

ncbi_parent <- left_join(
  ncbi_meta,
  nodes_df %>% select(tax_id, parent_tax_id, rank),
  by = join_by(organism.tax_id == tax_id)
)
ncbi_parent %>%
  filter(rank == "no rank") %>%
  left_join(nodes_df, by = join_by(parent_tax_id == tax_id)) %>%
  select(starts_with("rank"))

recursive_join <- function(ncbi_meta, taxonomy_nodes, iteration, join_col) {
  # TODO
  # or rather use rankedlineage.dmp
  ncbi_res <- left_join(
    ncbi_meta,
    taxonomy_nodes %>%
      select(tax_id, parent_tax_id, rank) %>%
      rename_with(.fn = ~ paste0(.x, ".", iteration), .cols = c(parent_tax_id, rank)),
    by = join_by(parent_tax_id == join_col),
  )
  rerun_ranks <- c("no rank", "species", "strain", "subspecies", "species group")
  if (any(ncbi_res$rank))

}
