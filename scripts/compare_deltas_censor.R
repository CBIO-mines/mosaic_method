library(dplyr)
library(tibble)
library(readr)
library(phangorn)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)


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

bacillales_path <- "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/"
censor_file_pattern <- "censor_x_reparam_global_results.csv"
censor_bacillales <- c("0.0", "1.7", "4.9", "7.9", "11.0", "14.0", "17.0")

# first few does not change much, then degrades
for (censor in censor_bacillales) {
  df_path <- paste0(bacillales_path, str_replace(censor_file_pattern, "x", censor))
  df <- read_csv(df_path, show_col_types = FALSE)
  dist <- get_distance(df, "theta")
  delta <- delta.plot(dist, plot=FALSE)
  delta_mean <- mean(delta$delta.bar)
  print(paste(censor, delta_mean))
}
