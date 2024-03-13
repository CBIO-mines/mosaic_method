library(stringr)
library(reticulate)
use_condaenv("test_florian")
fit <- import("fit")


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

mc_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(fit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[2]])
}


mh_fun <- function(r, par1, dr, mus, muc, d, L0) {
  r <- np_array(r)
  return(fit$fit$theoretical_mld(par1, dr, r, mus, muc, d, L0)[[1]])
}
