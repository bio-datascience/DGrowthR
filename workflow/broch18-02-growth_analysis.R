library(tidyverse)
library(DGrowthR)

dgobj.ecbw.logged <- readRDS("../Source_Data/Brochado_2018/dgobj_ecbw_logged.rds")

parsed <- dgobj.ecbw.logged@metadata %>%
  select(combination.id) %>%
  distinct() %>%
  separate(combination.id, into = c("left_raw", "right_raw"), sep = ":", remove = FALSE) %>%
  mutate(
    drug_left  = str_replace(left_raw, "\\.(\\d+|NA)$", ""),
    conc_left  = str_extract(left_raw, "(\\d+|NA)$"),
    drug_right = str_replace(right_raw, "\\.(\\d+|NA)$", ""),
    conc_right = str_extract(right_raw, "(\\d+|NA)$")
  )

special_drugs <- c("NoDurg", "Control", "SingleDrug")

parsed_conc1 <- parsed %>%
  filter(
    (conc_left == "1" | drug_left %in% special_drugs) &
      (conc_right == "1" | drug_right %in% special_drugs)
  )

parsed_conc1 <- parsed_conc1 %>%
  mutate(
    left_is_special  = drug_left %in% special_drugs,
    right_is_special = drug_right %in% special_drugs,

    combo_type = case_when(
      left_is_special & right_is_special ~ "control",
      left_is_special | right_is_special ~ "single_drug",
      TRUE ~ "combination"
    ),
    combo_label = case_when(
      combo_type == "control" ~ "Control",
      combo_type == "single_drug" & left_is_special  ~ drug_right,
      combo_type == "single_drug" & right_is_special ~ drug_left,
      combo_type == "combination" ~ paste(pmin(drug_left, drug_right), "+", pmax(drug_left, drug_right))
    )
  )

label_lookup <- parsed_conc1 %>%
  select(combination.id, combo_label, combo_type)

control_reps <- c("Control.rep.1", "Control.rep.2", "Control.rep.9",
                  "Control.rep.10", "Control.rep.11")

dgobj.ecbw.logged@metadata <- dgobj.ecbw.logged@metadata %>%
  left_join(label_lookup, by = "combination.id") %>%
  mutate(
    combination.label = case_when(
      is.na(combo_label) ~ "Other",
      combo_type == "control" ~ "Other",
      combo_type == "single_drug" & str_starts(combination.id, "Control\\.NA:") &
        !(query %in% control_reps) ~ "Other",
      TRUE ~ combo_label
    )
  ) %>%
  select(-combo_label, -combo_type)

active_labels <- unique(dgobj.ecbw.logged@metadata$combination.label)
active_labels <- active_labels[active_labels != "Other"]
combo_labels  <- active_labels[str_detect(active_labels, " \\+ ")]
single_labels <- setdiff(active_labels, combo_labels)

comparison_list <- list()
for (combo in combo_labels) {
  drugs <- str_split(combo, " \\+ ")[[1]]

  for (drug in drugs) {
    if (drug %in% single_labels) {
      comparison_list <- c(comparison_list,
                           list(c("combination.label", combo, drug)))
    }
  }
}

rdf <- multiple_comparisons(
  dgobj.ecbw.logged,
  comparison_list = comparison_list[1:2928],
  predict_n_steps = 50,
  downsample_every_n_timepoints = 3,
  permutation_test = TRUE,
  n_permutations = 500,
  n_cores = 20,
  save_perm_stats = TRUE
)

# Save results
saveRDS(rdf$result_df, "tmp_permutations/broch18_results_1-2928.rds")
saveRDS(rdf$permutation_df, "tmp_permutations/broch18_permutations_1-2928.rds")

laGP::deleteGPseps()

saveRDS(dgobj.ecbw.logged, "tmp_permutations/dgobj_broch18_mc.rds")
