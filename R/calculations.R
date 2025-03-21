#------------------------------------------------------------------------------------------------------------------
#' @title Calculate the Empirical p_value
#'
#' @description This function calculates the Empirical p-value (emp_p_value) based on a given object.
#'
#' @param object A DGrowthR object which contains the Gaussian Process regression results.
#'
#' @return A numeric value representing the calculated Empirical p-value.
#'
#' @export
setGeneric(
  name = "calculate_emp_p_value",
  def = function(object) {
    standardGeneric("calculate_emp_p_value")
  }
)

#' @rdname calculate_emp_p_value
#' @export
setMethod(
  f = "calculate_emp_p_value",
  signature = "DGrowthR",
  definition = function(object) {

    # Gather observed test statstic
    bf <- lrt_statistic(object)
    # Gather permutations
    perms <- perm_test_result(object)

    # Calculate empirical p-value
    n_larger <- sum(perms >= bf)
    emp_p_value <- (n_larger + 1) / (length(perms) + 1)
    return(emp_p_value)
  }
)


##------------------------------------------------------------------------------------------------------------------
#' Determine derivatives
#'
#' @param df A data.frame with timepoint and mean fields.
#' 
#' @return A data.frame with the original mean data, and the first and second derivatives
#'
#' @keywords internal 
.estimate_derivatives <- function(df){
  
  #print(df)
  
  # Ensure that data is ordered
  df <- df %>% 
    arrange(timepoint)
  
  # Estimate first derivative
  first_derivative <- diff(df$mean) / diff(df$timepoint)
  df$first_derivative <- c(first_derivative, NA)
  
  # Estimate second derivative
  second_derivative <- diff(first_derivative) /diff(df$timepoint[-1])
  df$second_derivative <- c(second_derivative, NA, NA)
  
  return(df)
  
}
