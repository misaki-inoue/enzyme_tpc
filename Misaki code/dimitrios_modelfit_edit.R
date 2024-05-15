enzyme <- read.csv("Database_FYP.csv")
gdh <- enzyme %>% 
  filter(originalid == "GDH001") %>%
  select(interactor1temp, originaltraitvalue) 
colnames(gdh) <- c("temp", "trait_value")
gdh_plot <- ggplot(gdh, aes(temp, trait_value)) + geom_point()
gdh_plot
dataset <- data.table::copy(gdh)

# Set the universal gas constant.
R <- 1.987

# Set the minimum trait measurement at the rise of the TPC as the 
# starting value for B_0.
B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

# If the dataset has measurements before the thermal optimum, 
# a starting value for DH_A can be set as the slope of the following 
# regression for that subset of the data:
#
# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
#
# Otherwise, just set DH_A to 15000.
T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

dataset_before_peak <- dataset[dataset$temp < T_pk,]
if ( 
  nrow(dataset_before_peak) < 3 || 
  length(unique(dataset_before_peak$temp)) < 3 || 
  length(unique(dataset_before_peak$trait_value)) < 3 
)
{
  DH_A_start <- 15000
} else
{
  y_vals <- dataset_before_peak$trait_value
  x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
  
  # Take the absolute value.
  DH_A_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
}

lm_dha <- lm(y_vals ~ x_vals)
plot(x_vals, y_vals)
abline(lm_dha, col='red')

# If the dataset has measurements after the thermal optimum, 
# a starting value for DH_H can be set as roughly the slope of the 
# following regression for that subset of the data:
#
# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
#
# Otherwise, just set DH_H to 50000.
dataset_after_peak <- dataset[dataset$temp > T_pk,]
if ( 
  nrow(dataset_after_peak) < 3 || 
  length(unique(dataset_after_peak$temp)) < 3 || 
  length(unique(dataset_after_peak$trait_value)) < 3 
)
{
  DH_H_start <- 50000
} else
{
  y_vals <- dataset_after_peak$trait_value
  x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
  
# Take the absolute value.
DH_H_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])

}

lm_dhh <- lm(y_vals ~ x_vals)
plot(x_vals, y_vals)
abline(lm_dhh, col='red')

# Set the starting value of DH_L arbitrarily to 70000.
DH_L_start <- 70000

# T_H50 should be close to the thermal optimum (T_pk), so set the 
# starting value of T_H50 = T_pk.
T_H50_start <- T_pk + 273.15

# Set the starting value of T_L50 arbitrarily to the minimum 
# temperature in the dataset.
T_L50_start <- min(dataset$temp) + 273.15

function_to_be_fitted <- function(B_0, DH_A, DH_L, T_L50, DH_H, T_H50, temp)
{
  temp <- temp + 273.15
  
  # Set the universal gas constant.
  R <- 1.987
  
  # Set the reference temperature to 0Â°C.
  T_ref <- 273.15
  
  
  result <- (B_0 * (temp/T_ref) * exp( (DH_A/R) * ((1/T_ref) - (1/temp)) ) ) / ( 1 + exp( (-DH_L/R) * ((1/T_L50) - (1/temp)) ) + exp( (DH_H/R) * ((1/T_H50) - (1/temp)) ) )
  
  if ( T_L50 >= T_H50 || is.na(result) || any(result <= 0) )
  {
    return(rep(1e10, length(temp)))
  } else
  {
    return(result)
  }
}

fit <- NULL

try(
  fit <- nls_multstart(
    trait_value ~ function_to_be_fitted(
      B_0, DH_A, DH_L, T_L50, DH_H, T_H50, temp = temp
    ),
    data = dataset,
    iter = 500,
    start_lower = c(
      B_0 = 0.5 * B_0_start,		DH_A = 0.5 * DH_A_start,
      DH_L = 0.5 * DH_L_start,	T_L50 = 0.5 * T_L50_start,
      DH_H = 0.5 * DH_H_start,	T_H50 = 0.5 * T_H50_start
    ),
    start_upper = c(
      B_0 = 1.5 * B_0_start,		DH_A = 1.5 * DH_A_start,
      DH_L = 1.5 * DH_L_start,	T_L50 = 1.5 * T_L50_start,
      DH_H = 1.5 * DH_H_start,	T_H50 = 1.5 * T_H50_start
    ),
    supp_errors = 'Y',
    convergence_count = FALSE,
    control = nls.lm.control(
      ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
      maxfev = 100000
    ),
    lower = c(0, 0, 0, -20 + 273.15, 0, 0 + 273.15),
    upper = c(Inf, 300000, 300000, 150 + 273.15, 300000, 150 + 273.15)
  )
)

return(fit)

plot(dataset$temp,dataset$trait_value)
new_data <- data.frame(temp = seq(min(dataset$temp), max(dataset$temp), 0.5))
preds <- augment(fit, newdata = new_data)

gdh_dimi <- ggplot(dataset, aes(temp, trait_value)) +
  geom_point(size=3) +
  geom_line(aes(temp, .fitted), preds, col = '#4f4789', linewidth=1)+
  ylab("Relative activity (%)") +
  xlab("Temperature (ºC)") +
  ggtitle("Dimitrios's Model fitting (Sharpe-Schoolfield)")
gdh_dimi

# Formula: trait_value ~ function_to_be_fitted(B_0, DH_A, DH_L, T_L50, DH_H, 
#                                              T_H50, temp = temp)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)
# B_0   8.616e-25  2.479e-18   0.000    1.000
# DH_A  3.000e+05  1.221e+10   0.000    1.000
# DH_L  7.483e-02  2.389e+10   0.000    1.000
# T_L50 2.532e+02  1.555e+13   0.000    1.000
# DH_H  3.000e+05  1.221e+10   0.000    1.000
# T_H50 3.068e+02  2.551e+05   0.001    0.999
# 
# Residual standard error: 57.2 on 9 degrees of freedom
# 
# Number of iterations to convergence: 3 
# Achieved convergence tolerance: 1.49e-08
