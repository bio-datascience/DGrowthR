library(DGrowthR)
library(tidyverse)


# Read dgobject
dg.logged <- readRDS("data/brenz24-01-preprocess_data/dgobj_logged.rds")


# Unique wells
unique_wells <-  unique(dg.logged@metadata$well)

# Create contrast list
generate_vector <- function(item) {
  c("genotype_well", paste0("dCBASS_", item), paste0("wt_", item))
}

contrasts_list <- lapply(unique_wells, generate_vector)


rdf <- multiple_comparisons(
  dg.logged,
  comparison_list = contrasts_list,
  predict_n_steps=100,
  permutation_test=TRUE,
  n_permutations=500,
  n_cores=20,
  save_perm_stats = TRUE
)

saveRDS(rdf$result_df, "permutation_results/brenzinger_500p.rds")
saveRDS(rdf$permutation_df, "permutation_results/brenzinger_dgresults.rds")

laGP::deleteGPseps()
