library(tidyverse)
library(jsonlite)


res_file <- "./bacillus_res_ani.tsv"
res_matrix_file <- "./bacillus_res_ani.tsv.matrix"
fastani_res <- read_tsv(res_file, col_names = c("query_genome_file", "reference_genome_file", "ANI", "n_bidirectional_fragment_mappings", "n_total_query_fragments")) %>%
  separate_wider_regex(query_genome_file, c(".*", "/", "query_accession" = ".{15}", "_", "query_assembly_name" = ".*", "_.*")) %>%
  separate_wider_regex(reference_genome_file, c(".*", "/", "reference_accession" = ".{15}", "_", "reference_assembly_name" = ".*", "_.*"))
ncbi_meta <- fromJSON("../ncbi_metadata/dataset_summary.json", flatten = TRUE)$reports

# analyze metadata
ncbi_meta %>%
  filter(grepl("Bacillus", average_nucleotide_identity.submitted_ani_match.organism_name)) %>%
  group_by(average_nucleotide_identity.match_status) %>%
  summarise(n = n())

# ani_res to distance matrix ----------------------------------------------------
all_paths <- bind_rows(read_tsv("./test_ani_query.txt", col_names = FALSE), read_tsv("./test_ani_ref.txt", col_names = FALSE)) %>%
  separate_wider_regex(X1, c("path1" = ".*", "/", "genome" = ".*")) %>%
  separate_wider_regex(genome, c("accession" = ".{15}", "_", "assembly_name" = ".*", "_.*"))
all_genomes <- all_paths$accession
min_ani <- min(fastani_res$ANI)

# TODO: input validation - make sure file has expected format & throw errors if it doesn't
# read in the first row to determine the matrix dimensions
matrix_dim <-
  as.numeric(utils::read.table(res_matrix_file, nrows = 1, as.is = TRUE))
# read in all the data from the lower triangle (exclude the first which is the matrix dim)
distance_matrix <- utils::read.table(res_matrix_file,
  fill = TRUE,
  skip = 1,
  # plus one to add a last column to make a square matrix
  col.names = c(as.character(1:(matrix_dim + 1)))
)

colnames(distance_matrix) <- c("rows", distance_matrix$X1)
rownames(distance_matrix) <- distance_matrix[, "rows"]
distance_matrix <- distance_matrix[, -1]
# save base matrix for test
ori_matrix <- distance_matrix

# remove paths/filenames from row and colnames
rownames_to_split <- rownames(distance_matrix)
colames_to_split <- colnames(distance_matrix)
genomes_rows <- str_extract(rownames_to_split, ".*/(GCF_[0-9]{9}.[0-9])", group = TRUE)
genomes_cols <- str_extract(rownames_to_split, ".*/(GCF_[0-9]{9}.[0-9])", group = TRUE)
colnames(distance_matrix) <- genomes_cols
rownames(distance_matrix) <- genomes_rows

# add t() to have a square matrix
transposed_triangular <- t(distance_matrix)
# replace NAs with 0 where transposed is not NA (NA + 85 = NA -> 0 + 85 = 85)
distance_matrix[!is.na(transposed_triangular)] <- 0
transposed_triangular[!is.na(ori_matrix)] <- 0
# finally sum
distance_matrix <- distance_matrix + transposed_triangular
# sanity check
stopifnot(isSymmetric(as.matrix(distance_matrix)))
test_distance_matrix <- distance_matrix
test_distance_matrix[upper.tri(distance_matrix)] <- NA
stopifnot(all(test_distance_matrix == ori_matrix, na.rm = TRUE))

# identity on the diagonal
diag(distance_matrix) <- 100

distance_matrix_matrix <- as.matrix(distance_matrix)
distance_matrix_matrix <- (100 - distance_matrix_matrix) / 100

# replace NAs ?? lets say min_ani - 5
min_ani_dist <- (95 - min_ani) / 100
dist_mat_final <- distance_matrix_matrix
dist_mat_final[is.na(dist_mat_final)] <- min_ani_dist
clust_res <- hclust(as.dist(dist_mat_final), method = "single")
png("hclust_ani.png")
plot(clust_res, labels = FALSE)
dev.off()
