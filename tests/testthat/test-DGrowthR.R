# Either use the provided example od_data
od_data <- DGrowthR::example_od_data
metadata <- DGrowthR::example_metadata
contrast <- c("Treatment", "Cancer_0.05", "DMSO_0.05")

# Create an instance of the DGrowthR class
test_object <- DGrowthRFromData(od_data = od_data, metadata = metadata)

test_that("DGrowthR S4: class can be created with valid inputs", {
 
  DGrowthR_obj <- DGrowthRFromData(od_data = od_data, metadata = metadata)

  expect_s4_class(DGrowthR_obj, "DGrowthR")
  expect_true(is.list(slot(DGrowthR_obj, "growth_comparison")))
  expect_true(is.list(slot(DGrowthR_obj, "gpfit_info")))
  expect_true(is.list(slot(DGrowthR_obj, "growth_comparison")))
  expect_true(is.list(slot(DGrowthR_obj, "fpca")))
  expect_true(is.data.frame(slot(DGrowthR_obj, "umap_coord")))
  expect_false(slot(DGrowthR_obj, "log_od"))
})

test_that("DGrowthR S4: class cannot be created with invalid inputs", {
  # Sample invalid data (missing 'Time' column)
  od_data_invalid <- data.frame(
    sample_1 = c(1, 2, 3),
    sample_2 = c(1, 2, 3)
  )

  expect_error(
    DGrowthRFromData(od_data_invalid, metadata),
    'No column called "Time" in od_data'
  )
})


# Test that we can do a growth comparison
test_that("growth_comparison works", {
  result <- growth_comparison(object=test_object, predict_n_steps = 10, comparison_info = contrast)
  expect_true(is(result, "DGrowthR"))
})

test_that("estimate_growth_parameters works", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(test_object_pre, save_gp_data=FALSE)
  expect_true(nrow(test_object_params@growth_parameters) == nrow(metadata))
})

test_that("estimate_growth_parameters by group works", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(test_object_pre, model_covariate="Treatment", predict_n_steps = 10)
  expect_true(nrow(test_object_params@growth_parameters) == 5)
})


test_that("estimate_growth_parameters with posterior sampling works", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(test_object_pre, 
                                                   model_covariate="Treatment", 
                                                   sample_posterior_gpfit=TRUE, 
                                                   predict_n_steps = 10,
                                                   sample_n_curves=5,
                                                   save_gp_data = FALSE)
  expect_true("max_growth_time.mean" %in% colnames(test_object_params@growth_parameters))
})

laGP::deleteGPs()
laGP::deleteGPseps()
