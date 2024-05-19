source("code/Prepocessing.R")
library(rTPC)
library(minpack.lm)  # Required for nls_multstart which uses the Levenberg-Marquardt algorithm
library(nls.multstart)
# Replace with the name of the package that includes nls_multstart

# Revised version of the model fitting function
fit_2nd_order_polynomial_3_pars <- function(dataset) {
  if (!check_data_availability(dataset)) {
    return(NULL)  # Skip this dataset if it does not meet criteria
  }
  
  a_start <- max(min(dataset$OriginalTraitValue, na.rm = TRUE), 1)  # Avoid Inf and ensure positive
  b_start <- 0.1  # Reasonable start assuming Temp has been scaled or normalized if necessary
  d_start <- 0.003
  
  model_formula <- log(OriginalTraitValue) ~ log(a + b * Temp + d * Temp^2)  # Verify this formula
  
  # Safeguard against unrealistic initial values causing model fitting issues
  if (is.infinite(a_start)) {
    warning("Infinite 'a_start' encountered; adjusting.")
    a_start <- 1
  }
  
  fit <- nls_multstart(
    formula = model_formula,
    data = dataset,
    iter = 1000,
    start_lower = c(a = 0.5 * a_start, b = 0.5 * b_start, d = 0.5 * d_start),
    start_upper = c(a = 1.5 * a_start, b = 1.5 * b_start, d = 1.5 * d_start),
    lower = c(a = -Inf, b = -Inf, d = -Inf),
    upper = c(a = Inf, b = Inf, d = Inf),
    control = nls.lm.control(
      maxiter = 1024, maxfev = 100000
    ),
    supp_errors = "Y"
  )
  
  return(fit)
}

# Function to check data availability and correctness
check_data_availability <- function(df) {
  necessary_columns <- c("OriginalTraitValue", "Temp")
  has_columns <- all(necessary_columns %in% names(df))
  has_data <- all(sapply(df[necessary_columns], function(col) sum(!is.na(col)) > 0))
  
  if (has_columns && has_data) {
    return(TRUE)
  } else {
    warning("Missing necessary columns or data in a subset")
    return(FALSE)
  }
}

# Apply this check
valid_subsets <- lapply(cleaned_subsets, check_data_availability)
cleaned_subsets <- cleaned_subsets[unlist(valid_subsets)]  # Keep only valid subsets

# Example usage
# Function to apply the fitting process to each subset
apply_model_to_subsets <- function(subsets) {
  models <- lapply(subsets, function(df) {
    if("OriginalTraitValue" %in% names(df) && "Temp" %in% names(df)) {
      fit_2nd_order_polynomial_3_pars(df)
    } else {
      warning("Skipped a subset due to missing columns.")
      return(NULL)
    }
  })
  models <- Filter(Negate(is.null), models)  # Remove NULL entries
  return(models)
}

# Apply to relevant_data_subsets
model_fits <- apply_model_to_subsets(cleaned_subsets)

# Step 1: Extract relevant information from model fits
model_summaries <- lapply(model_fits, function(model) {
  if (!is.null(model)) {
    summary_data <- coef(summary(model))
    aic <- AIC(model)
    return(data.frame(a = summary_data[1, "Estimate"],
                      b = summary_data[2, "Estimate"],
                      d = summary_data[3, "Estimate"],
                      se_a = summary_data[1, "Std. Error"],
                      se_b = summary_data[2, "Std. Error"],
                      se_d = summary_data[3, "Std. Error"],
                      AIC = aic,
                      stringsAsFactors = FALSE))
  }
})

# Step 2: Organize the information into a summary dataframe
summary_df <- do.call(rbind, model_summaries)

# Step 3: Write the summary dataframe to a CSV file
write.csv(summary_df, "Results/model_summary1.csv", row.names = FALSE)

#print the summary dataframe
print(summary_df)

########################################
# EXTENDED BRIERE MODEL (5 parameters) #
# (Cruz-Loya et al., bioRxiv, 2020)    #
#                                      #
# Parameters: a, T_min, T_max, b, d    #
########################################
# Revised version of the model fitting function
fit_extended_Briere_5_pars <- function(dataset)
{
  
  # Set the starting value of a arbitrarily to 1.
  a_start <- 1
  
  # Set the starting value of T_min 1 degree below the minimum 
  # temperature in the dataset.
  T_min_start <- min(dataset$temp) - 1
  
  # Set the starting value of T_max 1 degree above the maximum 
  # temperature in the dataset.
  T_max_start <- max(dataset$temp) + 1
  
  # Set the starting value of b arbitrarily to 1.
  b_start <- 1
  
  # Set the starting value of d arbitrarily to 1.5
  d_start <- 1.5
  
  function_to_be_fitted <- function(a, T_min, T_max, b, d, temp)
  {
    if (T_min >= T_max || any(temp < T_min) || any(temp > T_max))
    {
      return(rep(1e10, length(temp)))
    } else
    { 
      return(
        log(
          a * temp * ((temp - T_min)^b) * ((T_max - temp)^d)
        )
      )
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        a, T_min, T_max, b, d, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        a = 0.5 * a_start,      T_min = 0.5 * T_min_start,
        T_max = 0.5 * T_max_start,    b = 0.5 * b_start,
        d = 0.5 * d_start
      ),
      start_upper = c(
        a = 1.5 * a_start,      T_min = 1.5 * T_min_start,
        T_max = 1.5 * T_max_start,    b = 1.5 * b_start,
        d = 1.5 * d_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, -20, 0, 0, 0),
      upper = c(Inf, 150, 150, Inf, Inf)
    )
  )
  
  return(fit)
}

# Apply the fitting process to relevant data subsets
model_fits_extended_Briere <- apply_model_to_subsets(cleaned_subsets)

# Step 1: Extract relevant information from model fits
# Filter out NULL entries and only process valid model fits
valid_model_fits <- Filter(Negate(is.null), model_fits_extended_Briere)

# Process valid model fits
model_summaries_extended_Briere <- lapply(valid_model_fits, function(model) {
  summary_data <- coef(summary(model))
  aic <- AIC(model)
  
  # Check if summary_data has enough rows for all parameters
  if (nrow(summary_data) >= 5) {
    return(data.frame(a = summary_data[1, "Estimate"],
                      T_min = summary_data[2, "Estimate"],
                      T_max = summary_data[3, "Estimate"],
                      b = summary_data[4, "Estimate"],
                      d = summary_data[5, "Estimate"],
                      se_a = summary_data[1, "Std. Error"],
                      se_T_min = summary_data[2, "Std. Error"],
                      se_T_max = summary_data[3, "Std. Error"],
                      se_b = summary_data[4, "Std. Error"],
                      se_d = summary_data[5, "Std. Error"],
                      AIC = aic,
                      stringsAsFactors = FALSE))
  } else {
    # Handle case where summary_data doesn't have enough rows
    return(NULL)
  }
})

# Filter out NULL entries again
model_summaries_extended_Briere <- Filter(Negate(is.null), model_summaries_extended_Briere)

# Combine model summaries into a single dataframe
summary_df_extended_Briere <- do.call(rbind, model_summaries_extended_Briere)

# Write the summary dataframe to a CSV file
write.csv(summary_df_extended_Briere, "Results/model_summary_extended_Briere.csv", row.names = FALSE)

# Print the summary dataframe
print(summary_df_extended_Briere)





