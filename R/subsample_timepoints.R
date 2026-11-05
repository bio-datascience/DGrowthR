#' Subsample DGrowthR timepoints for faster inference
#'
#' @description
#' Creates a copy of a `DGrowthR` object with fewer observed OD timepoints.
#' This is intended for computationally expensive inference workflows such as
#' permutation tests, while keeping the full-resolution object available for
#' exploratory analysis and growth-parameter estimation.
#'
#' @param object A `DGrowthR` object.
#' @param keep_fraction Fraction of unique timepoints to keep. Ignored when
#' `n_timepoints` is supplied.
#' @param n_timepoints Number of unique timepoints to keep.
#' @param method Subsampling method. `"information"` greedily selects
#' timepoints that maximize GP mutual information between selected and
#' unselected locations following the sensor-placement criterion of Krause,
#' Singh, and Guestrin (2008). `"even"` keeps approximately evenly spaced
#' timepoints. `"growth_parameters"` prioritizes GP-derived landmark timepoints
#' already present in `object@growth_parameters`. `"gp_derivative"` prioritizes
#' high-derivative regions from `object@gpfit_info$gpfit_data`.
#' @param group_by Optional metadata column used to select dynamic landmarks per
#' group. This is useful when curve families have different dynamic windows.
#' @param gpfit_id_column Metadata column whose values match GP fit identifiers.
#' Defaults to `group_by` when `group_by` is supplied.
#' @param always_keep Numeric vector of timepoints that should always be kept.
#' The nearest observed timepoint is used for values not exactly present.
#' @param early_timepoints Numeric vector of early timepoints to keep. The
#' nearest observed timepoint is used for values not exactly present.
#' @param dynamic_fraction Fraction of the timepoint budget used for dynamic
#' landmarks when `method` is `"growth_parameters"` or `"gp_derivative"`.
#' @param information_lengthscale Optional Gaussian-kernel lengthscale used for
#' `method = "information"`. If `NULL`, a robust default is estimated from the
#' observed time range.
#' @param information_nugget Observation-noise nugget used for the
#' information-selection objective.
#' @param reset_model_slots If `TRUE`, clear model and embedding slots whose
#' values were generated on the full time grid.
#'
#' @return A `DGrowthR` object with `od_data` subset to the selected timepoints.
#' The selected timepoints are stored in `attr(result@od_data,
#' "subsampled_timepoints")`.
#'
#' @export
setGeneric(
  name = "subsample_timepoints",
  function(object,
           keep_fraction = 0.1,
           n_timepoints = NULL,
           method = c("information", "even", "growth_parameters", "gp_derivative"),
           group_by = NULL,
           gpfit_id_column = group_by,
           always_keep = NULL,
           early_timepoints = NULL,
           dynamic_fraction = 1 / 3,
           information_lengthscale = NULL,
           information_nugget = 1e-8,
           reset_model_slots = TRUE) {
    standardGeneric("subsample_timepoints")
  }
)

#' @rdname subsample_timepoints
#' @export
setMethod(
  f = "subsample_timepoints",
  signature = "DGrowthR",
  definition = function(object,
                        keep_fraction = 0.1,
                        n_timepoints = NULL,
                        method = c("information", "even", "growth_parameters", "gp_derivative"),
                        group_by = NULL,
                        gpfit_id_column = group_by,
                        always_keep = NULL,
                        early_timepoints = NULL,
                        dynamic_fraction = 1 / 3,
                        information_lengthscale = NULL,
                        information_nugget = 1e-8,
                        reset_model_slots = TRUE) {
    method <- match.arg(method)

    if (!is.null(n_timepoints) && (!is.numeric(n_timepoints) || length(n_timepoints) != 1)) {
      stop("n_timepoints must be a single numeric value.", call. = FALSE)
    }
    if (!is.null(n_timepoints) && n_timepoints < 2) {
      stop("n_timepoints must be at least 2.", call. = FALSE)
    }
    if (is.null(n_timepoints)) {
      if (!is.numeric(keep_fraction) || length(keep_fraction) != 1 ||
          keep_fraction <= 0 || keep_fraction > 1) {
        stop("keep_fraction must be a single value in (0, 1].", call. = FALSE)
      }
    }
    if (!is.numeric(dynamic_fraction) || length(dynamic_fraction) != 1 ||
        dynamic_fraction < 0 || dynamic_fraction > 1) {
      stop("dynamic_fraction must be a single value in [0, 1].", call. = FALSE)
    }
    if (!is.null(information_lengthscale) &&
        (!is.numeric(information_lengthscale) || length(information_lengthscale) != 1 ||
         information_lengthscale <= 0)) {
      stop("information_lengthscale must be NULL or a single positive numeric value.", call. = FALSE)
    }
    if (!is.numeric(information_nugget) || length(information_nugget) != 1 ||
        information_nugget < 0) {
      stop("information_nugget must be a single non-negative numeric value.", call. = FALSE)
    }

    observed_timepoints <- sort(unique(object@od_data$timepoint))
    if (length(observed_timepoints) < 2) {
      stop("object@od_data must contain at least two unique timepoints.", call. = FALSE)
    }

    if (is.null(n_timepoints)) {
      n_timepoints <- ceiling(length(observed_timepoints) * keep_fraction)
    }
    n_timepoints <- max(2, min(length(observed_timepoints), as.integer(n_timepoints)))

    nearest_observed_timepoint <- function(x) {
      if (length(x) == 0 || all(is.na(x))) {
        return(numeric(0))
      }
      x <- x[!is.na(x)]
      vapply(
        x,
        function(value) observed_timepoints[which.min(abs(observed_timepoints - value))],
        numeric(1)
      )
    }

    evenly_spaced_timepoints <- function(n) {
      observed_timepoints[unique(round(seq(1, length(observed_timepoints), length.out = n)))]
    }

    priority <- unique(c(
      nearest_observed_timepoint(always_keep),
      nearest_observed_timepoint(early_timepoints)
    ))

    selected <- priority
    n_dynamic <- max(1, floor(n_timepoints * dynamic_fraction))

    if (identical(method, "information")) {
      selected <- .subsample_information_timepoints(
        observed_timepoints = observed_timepoints,
        n_timepoints = n_timepoints,
        initial_timepoints = priority,
        lengthscale = information_lengthscale,
        nugget = information_nugget
      )
    }

    if (identical(method, "growth_parameters")) {
      growth_parameters <- object@growth_parameters
      landmark_columns <- intersect(
        c(
          "max_growth_time",
          "max_growth_rate_timepoint",
          "max_death_rate_timepoint",
          "lag_phase_timepoint",
          "stationary_phase_timepoint"
        ),
        colnames(growth_parameters)
      )

      if (!is.data.frame(growth_parameters) || nrow(growth_parameters) == 0 ||
          !"gpfit_id" %in% colnames(growth_parameters) || length(landmark_columns) == 0) {
        stop(
          "method = 'growth_parameters' requires object@growth_parameters with gpfit_id and GP landmark timepoint columns. ",
          "Run estimate_growth_parameters() first or use method = 'even'.",
          call. = FALSE
        )
      }

      landmark_data <- growth_parameters %>%
        dplyr::select(dplyr::all_of(c("gpfit_id", landmark_columns))) %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(landmark_columns),
          names_to = ".landmark",
          values_to = "timepoint"
        ) %>%
        dplyr::filter(!is.na(timepoint)) %>%
        dplyr::mutate(gpfit_id = as.character(gpfit_id))

      if (!is.null(group_by)) {
        group_map <- .subsample_gpfit_group_map(object, group_by, gpfit_id_column)
        group_count <- max(1, dplyr::n_distinct(group_map$.group_value))
        per_group_n <- max(1, ceiling(n_dynamic / group_count))

        landmark_data <- landmark_data %>%
          dplyr::inner_join(group_map, by = "gpfit_id") %>%
          dplyr::count(.group_value, timepoint, name = ".score") %>%
          dplyr::group_by(.group_value) %>%
          dplyr::arrange(dplyr::desc(.score), timepoint, .by_group = TRUE) %>%
          dplyr::slice_head(n = per_group_n) %>%
          dplyr::ungroup()
      } else {
        landmark_data <- landmark_data %>%
          dplyr::count(timepoint, name = ".score") %>%
          dplyr::arrange(dplyr::desc(.score), timepoint) %>%
          dplyr::slice_head(n = n_dynamic)
      }

      selected <- unique(c(selected, nearest_observed_timepoint(landmark_data$timepoint)))
    }

    if (identical(method, "gp_derivative")) {
      gpfit_data <- object@gpfit_info$gpfit_data
      if (is.null(gpfit_data) || !is.data.frame(gpfit_data) ||
          !all(c("timepoint", "first_derivative", "gpfit_id") %in% colnames(gpfit_data))) {
        stop(
          "method = 'gp_derivative' requires object@gpfit_info$gpfit_data with timepoint, first_derivative, and gpfit_id columns. ",
          "Run estimate_growth_parameters(..., save_gp_data = TRUE) or use method = 'growth_parameters'/'even'.",
          call. = FALSE
        )
      }

      dynamic_data <- gpfit_data %>%
        dplyr::filter(!is.na(first_derivative), !is.na(timepoint)) %>%
        dplyr::mutate(
          gpfit_id = as.character(gpfit_id),
          .score = abs(first_derivative)
        )

      if (!is.null(group_by)) {
        group_map <- .subsample_gpfit_group_map(object, group_by, gpfit_id_column)
        group_count <- max(1, dplyr::n_distinct(group_map$.group_value))
        per_group_n <- max(1, ceiling(n_dynamic / group_count))

        dynamic_data <- dynamic_data %>%
          dplyr::inner_join(group_map, by = "gpfit_id") %>%
          dplyr::group_by(.group_value) %>%
          dplyr::slice_max(.score, n = per_group_n, with_ties = FALSE) %>%
          dplyr::ungroup()
      } else {
        dynamic_data <- dynamic_data %>%
          dplyr::slice_max(.score, n = n_dynamic, with_ties = FALSE)
      }

      selected <- unique(c(selected, nearest_observed_timepoint(dynamic_data$timepoint)))
    }

    if (!identical(method, "information")) {
      selected <- unique(c(selected, evenly_spaced_timepoints(n_timepoints)))
    }

    if (length(selected) > n_timepoints) {
      selected_priority <- selected[selected %in% priority]
      selected_other <- selected[!selected %in% priority]
      if (length(selected_priority) >= n_timepoints) {
        selected <- selected_priority[seq_len(n_timepoints)]
      } else {
        selected <- c(
          selected_priority,
          selected_other[seq_len(n_timepoints - length(selected_priority))]
        )
      }
    }

    selected <- sort(unique(selected))

    result <- object
    slot(result, "od_data") <- object@od_data %>%
      dplyr::filter(timepoint %in% selected) %>%
      dplyr::arrange(curve_id, timepoint) %>%
      dplyr::group_by(curve_id) %>%
      dplyr::mutate(timepoint_n = dplyr::row_number()) %>%
      dplyr::ungroup()

    attr(result@od_data, "subsampled_timepoints") <- selected

    if (reset_model_slots) {
      slot(result, "fpca") <- list()
      slot(result, "umap_coord") <- data.frame()
      slot(result, "cluster_assignment") <- data.frame()
      slot(result, "gpfit_info") <- list()
      slot(result, "growth_comparison") <- list()
      slot(result, "growth_parameters") <- data.frame()
    }

    methods::validObject(result)
    result
  }
)

#' Greedy GP mutual-information timepoint selection
#'
#' @param observed_timepoints Sorted observed timepoints.
#' @param n_timepoints Number of timepoints to keep.
#' @param initial_timepoints Timepoints selected before the greedy pass.
#' @param lengthscale Optional Gaussian-kernel lengthscale.
#' @param nugget Observation-noise nugget.
#'
#' @return A numeric vector of selected timepoints.
#'
#' @references
#' Krause, A., Singh, A., & Guestrin, C. (2008). Near-optimal sensor placements
#' in Gaussian processes: Theory, efficient algorithms and empirical studies.
#' Journal of Machine Learning Research, 9, 235-284.
#'
#' @keywords internal
.subsample_information_timepoints <- function(observed_timepoints, n_timepoints,
                                              initial_timepoints = numeric(0),
                                              lengthscale = NULL,
                                              nugget = 1e-8) {
  n_observed <- length(observed_timepoints)
  if (is.null(lengthscale)) {
    time_range <- diff(range(observed_timepoints, na.rm = TRUE))
    lengthscale <- ifelse(time_range > 0, time_range^2 / 3, 1)
  }

  dist_sq <- outer(observed_timepoints, observed_timepoints, "-")^2
  covariance <- exp(-dist_sq / lengthscale)
  initial_indices <- match(initial_timepoints, observed_timepoints)
  selected_indices <- unique(initial_indices[!is.na(initial_indices)])

  if (length(selected_indices) >= n_timepoints) {
    return(sort(observed_timepoints[selected_indices[seq_len(n_timepoints)]]))
  }

  candidate_indices <- setdiff(seq_len(n_observed), selected_indices)
  while (length(selected_indices) < n_timepoints && length(candidate_indices) > 0) {
    objective_values <- vapply(
      candidate_indices,
      function(candidate_index) {
        .subsample_mutual_information(
          covariance = covariance,
          selected_indices = c(selected_indices, candidate_index),
          nugget = nugget
        )
      },
      numeric(1)
    )

    best_position <- which.max(objective_values)
    selected_indices <- c(selected_indices, candidate_indices[best_position])
    candidate_indices <- candidate_indices[-best_position]
  }

  sort(observed_timepoints[selected_indices])
}

#' Calculate GP mutual information between selected and unselected locations
#'
#' @param covariance Covariance matrix on the full observed grid.
#' @param selected_indices Integer indices of selected locations.
#' @param nugget Observation-noise nugget.
#'
#' @return Numeric mutual information value.
#'
#' @keywords internal
.subsample_mutual_information <- function(covariance, selected_indices, nugget) {
  all_indices <- seq_len(nrow(covariance))
  unselected_indices <- setdiff(all_indices, selected_indices)

  if (length(selected_indices) == 0 || length(unselected_indices) == 0) {
    return(0)
  }

  K_ss <- covariance[selected_indices, selected_indices, drop = FALSE]
  K_uu <- covariance[unselected_indices, unselected_indices, drop = FALSE]
  K_us <- covariance[unselected_indices, selected_indices, drop = FALSE]
  diag(K_ss) <- diag(K_ss) + nugget

  conditional_covariance <- K_uu - K_us %*% solve(K_ss, t(K_us))
  conditional_covariance <- (conditional_covariance + t(conditional_covariance)) / 2

  0.5 * (
    .subsample_logdet_psd(K_uu) -
      .subsample_logdet_psd(conditional_covariance)
  )
}

#' Numerically stable log determinant for positive semidefinite matrices
#'
#' @param matrix A numeric symmetric matrix.
#'
#' @return Numeric log determinant.
#'
#' @keywords internal
.subsample_logdet_psd <- function(matrix) {
  matrix <- (matrix + t(matrix)) / 2
  jitter <- 0
  for (i in seq_len(8)) {
    chol_result <- tryCatch(
      chol(matrix + diag(jitter, nrow(matrix))),
      error = function(cond) NULL
    )
    if (!is.null(chol_result)) {
      return(2 * sum(log(diag(chol_result))))
    }
    jitter <- if (i == 1) sqrt(.Machine$double.eps) else jitter * 10
  }

  -Inf
}

#' Map GP fit IDs to metadata groups for subsampling
#'
#' @param object A `DGrowthR` object.
#' @param group_by Metadata column used for grouping.
#' @param gpfit_id_column Metadata column containing GP fit IDs.
#'
#' @return A data frame with `.group_value` and `gpfit_id`.
#'
#' @keywords internal
.subsample_gpfit_group_map <- function(object, group_by, gpfit_id_column) {
  if (!group_by %in% colnames(object@metadata)) {
    stop("group_by is not a column in object@metadata.", call. = FALSE)
  }
  if (is.null(gpfit_id_column) || !gpfit_id_column %in% colnames(object@metadata)) {
    stop("gpfit_id_column must be a metadata column when group_by is supplied.", call. = FALSE)
  }

  object@metadata %>%
    dplyr::transmute(
      .group_value = as.character(.data[[group_by]]),
      gpfit_id = as.character(.data[[gpfit_id_column]])
    ) %>%
    dplyr::distinct()
}
