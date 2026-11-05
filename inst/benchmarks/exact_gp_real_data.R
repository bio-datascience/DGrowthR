# Compare DGrowthR GP backends on bundled real 3-replicate examples.
#
# Run from the package root after installing/loading DGrowthR development
# sources. This script intentionally uses the package's measured example data,
# not synthetic curves.

suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
})

if (file.exists("src/exact_gaussian_gp.so")) {
  dyn.load("src/exact_gaussian_gp.so")
}

source("R/gp_regression.R")

load("data/example_od_data.rda")
load("data/example_metadata.rda")

time_col <- example_od_data[[1]]
curve_ids <- colnames(example_od_data)[-1]
curve_auc <- vapply(curve_ids, function(id) sum(example_od_data[[id]]), numeric(1))

curve_info <- data.frame(
  curve_id = curve_ids,
  auc = curve_auc,
  stringsAsFactors = FALSE
) %>%
  left_join(example_metadata, by = "curve_id")

select_triplet <- function(ids, auc, label) {
  ids <- ids[order(auc)]
  n <- length(ids)

  switch(label,
         low_auc = ids[1:3],
         mid_auc = {
           mid <- floor(n / 2)
           ids[(mid - 1):(mid + 1)]
         },
         high_auc = ids[(n - 2):n])
}

examples <- curve_info %>%
  group_by(Treatment) %>%
  group_modify(~{
    bind_rows(lapply(c("low_auc", "mid_auc", "high_auc"), function(label) {
      ids <- select_triplet(.x$curve_id, .x$auc, label)
      data.frame(
        subset = label,
        curve_ids = paste(ids, collapse = ","),
        subset_auc = mean(.x$auc[match(ids, .x$curve_id)]),
        stringsAsFactors = FALSE
      )
    }))
  }) %>%
  ungroup()

make_curve_data <- function(ids) {
  ids <- strsplit(ids, ",", fixed = TRUE)[[1]]
  data.frame(
    timepoint = rep(time_col, length(ids)),
    timepoint_n = rep(seq_along(time_col), length(ids)),
    od = as.numeric(as.matrix(example_od_data[, ids]))
  )
}

fit_backend <- function(df, backend, strategy) {
  invisible(gc())
  t0 <- proc.time()[["elapsed"]]
  fit <- fit_predict_gp(
    df = df,
    t_steps = length(unique(df$timepoint)),
    gp_backend = backend,
    gp_replicate_strategy = strategy,
    optimize_hyperparameters = TRUE,
    prepare_dataframe = TRUE,
    estimate_derivatives = FALSE,
    delete = TRUE
  )

  list(fit = fit, elapsed_ms = 1000 * (proc.time()[["elapsed"]] - t0))
}

benchmark_one <- function(i) {
  ex <- examples[i, ]
  df <- make_curve_data(ex$curve_ids)
  lagp <- fit_backend(df, "lagp", "raw")
  exact_raw <- fit_backend(df, "exact_gaussian", "raw")
  exact_aggregate <- fit_backend(df, "exact_gaussian", "aggregate")

  mean_pred <- function(x) x$fit$prediction$mean

  data.frame(
    Treatment = ex$Treatment,
    subset = ex$subset,
    subset_auc = ex$subset_auc,
    lagp_ms = lagp$elapsed_ms,
    exact_raw_ms = exact_raw$elapsed_ms,
    exact_aggregate_ms = exact_aggregate$elapsed_ms,
    exact_raw_rmse = sqrt(mean((mean_pred(exact_raw) - mean_pred(lagp))^2)),
    exact_aggregate_rmse = sqrt(mean((mean_pred(exact_aggregate) - mean_pred(lagp))^2)),
    exact_raw_convergence = exact_raw$fit$mle$conv,
    exact_aggregate_convergence = exact_aggregate$fit$mle$conv
  )
}

results <- bind_rows(lapply(seq_len(nrow(examples)), benchmark_one))

print(results)
print(summarise(
  results,
  n = n(),
  median_lagp_ms = median(lagp_ms),
  median_exact_raw_ms = median(exact_raw_ms),
  median_exact_aggregate_ms = median(exact_aggregate_ms),
  median_exact_raw_rmse = median(exact_raw_rmse),
  max_exact_raw_rmse = max(exact_raw_rmse),
  median_exact_aggregate_rmse = median(exact_aggregate_rmse),
  max_exact_aggregate_rmse = max(exact_aggregate_rmse)
))
