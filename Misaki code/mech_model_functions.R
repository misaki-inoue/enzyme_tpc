##########################################
# VAN'T HOFF MODEL (4 parameters)        #
# (Portner et al., Biogeosciences, 2010) #
#                                        #
# Parameters: a, b, d, e                 #
##########################################
fit_Vant_Hoff_4_pars <- function(dataset)
{
  
  # Set the starting value of a arbitrarily to 1.
  a_start <- 1
  
  # Set the starting value of b arbitrarily to 0.6.
  b_start <- 0.6
  
  # Set the starting value of d arbitrarily to 1.
  d_start <- 1
  
  # Set the starting value of e arbitrarily to 2.
  e_start <- 2
  
  function_to_be_fitted <- function(a, b, d, e, temp)
  {
    temp <- temp + 273.15
    
    if ( b == 0 || d == 0 || e == 0 )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(
        log(
          a * exp(-b/temp) * (temp^d) * exp(e*temp) 
        )
      )
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        a, b, d, e, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        a = 0.5 * a_start,	b = 0.5 * b_start,
        d = 0.5 * d_start,	e = 0.5 * e_start
      ),
      start_upper = c(
        a = 1.5 * a_start,	b = 1.5 * b_start,
        d = 1.5 * d_start,	e = 1.5 * e_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, -Inf, -Inf, -Inf),
      upper = c(Inf, Inf, Inf, Inf)
    )
  )
  
  return(fit)
}

####################################
# HINSHELWOOD MODEL (4 parameters) #
# (Hinshelwood, 1946)              #
#                                  #
# Parameters: a, E_1, b, E_2       #
####################################
fit_Hinshelwood_4_pars <- function(dataset)
{
  
  # Set the universal gas constant.
  R <- 0.001987
  
  # Set the starting value of a to the minimum trait value at the 
  # low temperature end.
  a_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
  
  # Set the starting value of b to the minimum trait value at the 
  # high temperature end.
  b_start <- min(dataset$trait_value[dataset$temp == max(dataset$temp)])
  
  # If the dataset has measurements before the thermal optimum, 
  # a starting value for E_1 can be set as roughly the slope of the 
  # following regression for that subset of the data:
  #
  # ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15)))
  #
  # Otherwise, just set E_1 to 0.6.
  T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
  
  dataset_before_peak <- dataset[dataset$temp < T_pk,]
  if ( 
    nrow(dataset_before_peak) < 3 || 
    length(unique(dataset_before_peak$temp)) < 3 || 
    length(unique(dataset_before_peak$trait_value)) < 3 
  )
  {
    E_1_start <- 0.6
  } else
  {
    y_vals <- log(dataset_before_peak$trait_value)
    x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
    
    # Take the absolute value.
    E_1_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
  }
  
  # If the dataset has measurements after the thermal optimum, 
  # a starting value for E_2 can be set as the slope of the following 
  # regression for that subset of the data:
  #
  # ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
  #
  # Otherwise, just set E_2 to 3.
  dataset_after_peak <- dataset[dataset$temp > T_pk,]
  if ( 
    nrow(dataset_after_peak) < 3 || 
    length(unique(dataset_after_peak$temp)) < 3 || 
    length(unique(dataset_after_peak$trait_value)) < 3 
  )
  {
    E_2_start <- 3
  } else
  {
    y_vals <- log(dataset_after_peak$trait_value)
    x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
    
    # Take the absolute value.
    E_2_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
  }
  
  function_to_be_fitted <- function(a, E_1, b, E_2, temp)
  {
    temp <- temp + 273.15
    
    # Set the universal gas constant.
    R <- 0.001987
    
    result <- a * exp( -E_1 / (R * temp) ) - b * exp( -E_2 / (R * temp) )
    
    if ( any(result <= 0 ) || E_1 == 0 || E_2 == 0 )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(log(result))
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        a, E_1, b, E_2, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        a = 0.5 * a_start,	E_1 = 0.5 * E_1_start,
        b = 0.5 * b_start,	E_2 = 0.5 * E_2_start
      ),
      start_upper = c(
        a = 1.5 * a_start,	E_1 = 1.5 * E_1_start,
        b = 1.5 * b_start,	E_2 = 1.5 * E_2_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, 0, 0, 0),
      upper = c(Inf, Inf, Inf, Inf)
    )
  )
  
  return(fit)
}

####################################################
# JOHNSON - LEWIN MODEL (4 parameters)             #
# (Johnson & Lewin, J. Cell. Comp. Physiol., 1946) #
#                                                  #
# Parameters: B_0, E, T_pk, E_D                    #
####################################################
fit_Johnson_Lewin_4_pars <- function(dataset)
{
  
  # Set the Boltzmann constant.
  k <- 8.617 * 10^-5
  
  # Set the minimum trait measurement at the rise of the TPC as the 
  # starting value for B_0.
  B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
  
  # Set the starting value of T_pk to the temperature at which the 
  # trait reaches its maximum value.
  T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
  
  # If the dataset has measurements before the thermal optimum, 
  # a starting value for E can be set as roughly the slope of the 
  # following regression for that subset of the data:
  #
  # ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
  #
  # Otherwise, just set E to 0.6 eV.
  
  dataset_before_peak <- dataset[dataset$temp < T_pk_start,]
  if ( 
    nrow(dataset_before_peak) < 3 || 
    length(unique(dataset_before_peak$temp)) < 3 || 
    length(unique(dataset_before_peak$trait_value)) < 3 
  )
  {
    E_start <- 0.6
  } else
  {
    y_vals <- log(dataset_before_peak$trait_value)
    x_vals <- 1/(k * (dataset_before_peak$temp + 273.15))
    
    # Take the absolute value.
    E_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
    
    if ( E_start >= 10 )
    {
      E_start <- 0.6
    }
  }
  
  # If the dataset has measurements after the thermal optimum, 
  # a starting value for E_D can be set as the slope of the following 
  # regression for that subset of the data:
  #
  # ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
  #
  # Otherwise, just set E_D to 3 eV.
  dataset_after_peak <- dataset[dataset$temp > T_pk_start,]
  if ( 
    nrow(dataset_after_peak) < 3 || 
    length(unique(dataset_after_peak$temp)) < 3 || 
    length(unique(dataset_after_peak$trait_value)) < 3 
  )
  {
    E_D_start <- 3
  } else
  {
    y_vals <- log(dataset_after_peak$trait_value)
    x_vals <- 1/(k * (dataset_after_peak$temp + 273.15))
    
    # Take the absolute value.
    E_D_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
    
    if ( E_D_start >= 50 )
    {
      E_D_start <- 3
    }
  }
  
  if ( E_start >= E_D_start )
  {
    E_start <- 0.9 * E_D_start
  }
  
  function_to_be_fitted <- function(B_0, E, T_pk, E_D, temp)
  {
    temp <- temp + 273.15
    
    # Set the Boltzmann constant.
    k <- 8.617 * 10^-5
    
    if ( E == E_D )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(
        log(
          B_0 * temp * exp((-E/k) * (1/temp)) / ( 1 + ( E/( E_D - E ) ) * exp( (E_D/k) * (1/T_pk - 1/temp) ) )
        )
      )
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        B_0, E, T_pk, E_D, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        B_0 = 0.5 * B_0_start,				E = 0.5 * E_start,
        T_pk = 0.5 * T_pk_start + 273.15,	E_D = 0.5 * E_D_start
      ),
      start_upper = c(
        B_0 = 1.5 * B_0_start,				E = 1.5 * E_start,
        T_pk = 1.5 * T_pk_start + 273.15,	E_D = 1.5 * E_D_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, 0, 273.15, 0),
      upper = c(Inf, 10, 273.15 + 150, 50)
    )
  )
  
  return(fit)
}

###################################################
# SHARPE - SCHOOLFIELD MODEL (6 parameters)       #
# (Schoolfield et al., J. Theor. Biol., 1981)     #
#                                                 #
# Parameters: B_0, DH_A, DH_L, T_L50, DH_H, T_H50 #
###################################################
fit_Sharpe_Schoolfield_6_pars <- function(dataset)
{
  
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
    y_vals <- log(dataset_before_peak$trait_value)
    x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
    
    # Take the absolute value.
    DH_A_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
  }
  
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
    y_vals <- log(dataset_after_peak$trait_value)
    x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
    
    # Take the absolute value.
    DH_H_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
  }
  
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
    
    # Set the reference temperature to 0°C.
    T_ref <- 273.15
    
    result <- (B_0 * (temp/T_ref) * exp( (DH_A/R) * ((1/T_ref) - (1/temp)) ) ) / ( 1 + exp( (-DH_L/R) * ((1/T_L50) - (1/temp)) ) + exp( (DH_H/R) * ((1/T_H50) - (1/temp)) ) )
    
    if ( T_L50 >= T_H50 || any(is.nan(result)) || any(result <= 0) )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(log(result))
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        B_0, DH_A, DH_L, T_L50, DH_H, T_H50, temp = temp
      ),
      data = dataset,
      iter = 1000,
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
}

#############################################
# ROSS - RATKOWSKY MODEL (5 parameters)     #
# (Ratkowsky et al., J. Theor. Biol., 2005) #
#                                           #
# Parameters: a, DH_A, n, DH, DC_p          #
#############################################
fit_Ross_Ratkowsky_5_pars <- function(dataset)
{
  
  # Set the universal gas constant.
  R <- 8.314
  
  # Set the starting value of a arbitrarily to 1.
  a_start <- 1
  
  # Set the starting value of DH_A arbitrarily to 50000.
  DH_A_start <- 50000
  
  # Set the starting value of n arbitrarily to 300.
  n_start <- 300
  
  # Set the starting value of DH arbitrarily to 5000.
  DH_start <- 5000
  
  # Set the starting value of DC_p arbitrarily to 60.
  DC_p_start <- 60
  
  function_to_be_fitted <- function(a, DH_A, n, DH, DC_p, temp)
  {
    temp <- temp + 273.15
    
    # Set the universal gas constant.
    R <- 8.314
    
    D <- 1 + exp( -n * ( (DH - temp * 18.1 + DC_p * ( (temp - 373.6) - temp * log(temp/385.2) ) ) / (R * temp) ) )
    result <- ( a * temp * exp((-DH_A/(R * temp)) ) ) / D
    
    if ( any(is.nan(result)) || any(result <= 0 ) || a == 0 )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(log(result))
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        a, DH_A, n, DH, DC_p, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        a = 0.5 * a_start,		DH_A = 0.5 * DH_A_start,
        n = 0.5 * n_start,		DH = 0.5 * DH_start,
        DC_p = 0.5 * DC_p_start
      ),
      start_upper = c(
        a = 1.5 * a_start,		DH_A = 1.5 * DH_A_start,
        n = 1.5 * n_start,		DH = 1.5 * DH_start,
        DC_p = 1.5 * DC_p_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, 10000, 10, 3000, 37),
      upper = c(Inf, 120000, 20000, 7000, 118)
    )
  )
  
  return(fit)
}

#########################################
# HOBBS MODEL (4 parameters)            #
# (Hobbs et al., ACS Chem. Biol., 2013) #
#                                       #
# Parameters: B_0, DH, DC_p, DS         #
#########################################
fit_Hobbs_4_pars <- function(dataset)
{
  
  # Set the minimum trait measurement at the rise of the TPC as the 
  # starting value for B_0.
  B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
  
  # Set the starting value of DH arbitrarily to 30000.
  DH_start <- 30000
  
  # Set the starting value of DC_p arbitrarily to 1000.
  DC_p_start <- 1000
  
  # Set the starting value of DS arbitrarily to 17.
  DS_start <- 17
  
  function_to_be_fitted <- function(B_0, DH, DC_p, DS, temp)
  {		
    temp <- temp + 273.15
    
    # Set the reference temperature to 0°C.
    T_ref <- 273.15
    
    # Set the universal gas constant.
    R <- 8.314
    
    # Set the Boltzmann constant.
    k <- 1.38 * 1e-23
    
    # Set the Planck constant.
    h <- 6.626 * 1e-34
    
    result <- B_0 * ( (k * temp) / h ) * exp( ( - ( DH - DC_p * (temp - T_ref) ) / (R * temp) ) + ( (DS - DC_p * log(temp/T_ref)) / R ) ) 
    
    if ( any(is.nan(result)) || any(result <= 0 ) || B_0 == 0 || DH == 0 || DC_p == 0 || DS == 0 )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(log(result))
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        B_0, DH, DC_p, DS, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        B_0 = 0.5 * B_0_start,		DH = 0.5 * DH_start,
        DC_p = 0.5 * DC_p_start,	DS = 0.5 * DS_start
      ),
      start_upper = c(
        B_0 = 1.5 * B_0_start,		DH = 1.5 * DH_start,
        DC_p = 1.5 * DC_p_start,	DS = 1.5 * DS_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, 0, 0, 0),
      upper = c(Inf, Inf, Inf, Inf)
    )
  )
  
  return(fit)
}

####################################################
# ENZYME-ASSISTED ARRHENIUS MODEL (5 parameters)   #
# (DeLong et al., Ecol. Evol., 2017)               #
#                                                  #
# Parameters: a, E_b, E_DH, T_m, E_DCp             #
####################################################
fit_Enzyme_Assisted_Arrhenius_5_pars <- function(dataset)
{
  
  # Set the starting value of a arbitrarily to 1.
  a_start <- 1
  
  # Set the starting value of E_b arbitrarily to 0.01.
  E_b_start <- 0.01
  
  # Set the starting value of E_DH arbitrarily to 2.
  E_DH_start <- 2
  
  # Set the starting value of E_DCp arbitrarily to 0.1.
  E_DCp_start <- 0.1
  
  # Set the starting value of T_m 1 degree above the highest 
  # temperature in the dataset.
  T_m_start <- max(dataset$temp) + 1
  
  function_to_be_fitted <- function(a, E_b, E_DH, E_DCp, T_m, temp)
  {
    temp <- temp + 273.15
    T_m <- T_m + 273.15
    
    # Set the Boltzmann constant.
    k <- 8.617 * 10^-5
    
    if ( any(temp > T_m) )
    {
      return(rep(1e10, length(temp)))
    } else
    {		
      return(
        log(
          a * exp( (- (E_b - (E_DH * (1 - (temp/T_m) ) + E_DCp * (temp - T_m - temp * log(temp/T_m) ) ) ) ) / (k * temp) ) 
        )
      )
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        a, E_b, E_DH, E_DCp, T_m, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        a = 0.5 * a_start,			E_b = 0.5 * E_b_start,
        E_DH = 0.5 * E_DH_start,	E_DCp = 0.5 * E_DCp_start,
        T_m = 0.5 * T_m_start
      ),
      start_upper = c(
        a = 1.5 * a_start,			E_b = 1.5 * E_b_start,
        E_DH = 1.5 * E_DH_start,	E_DCp = 1.5 * E_DCp_start,
        T_m = 1.5 * T_m_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, 0, -Inf, 0, 0),
      upper = c(Inf, Inf, Inf, Inf, 150)
    )
  )
  
  return(fit)
}

################################
# RITCHIE MODEL (4 parameters) #
# (Ritchie, Sci. Rep., 2018)   #
#                              #
# Parameters: d_0, E_D, DE, w  #
################################
fit_Ritchie_4_pars <- function(dataset)
{
  
  # Set the minimum trait measurement at the rise of the TPC as the 
  # starting value for d_0.
  d_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
  
  # Set the starting value of E_D arbitrarily to 20.
  E_D_start <- 20
  
  # Set the starting value of DE arbitrarily to 30.
  DE_start <- 30
  
  # Set the starting value of w arbitrarily to 11.
  w_start <- 11
  
  function_to_be_fitted <- function(d_0, E_D, DE, w, temp)
  {
    temp <- temp + 273.15
    
    R <- 8.318 * 1e-3
    
    if ( any( ( DE/(R * temp) ) <= w ) )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(
        log(
          R * d_0 * exp( (-E_D) / (R * temp) ) * ((DE/(R * temp)) - w)
        )
      )
    }
  }
  
  fit <- NULL
  
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        d_0, E_D, DE, w, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        d_0 = 0.5 * d_0_start,	E_D = 0.5 * E_D_start,
        DE = 0.5 * DE_start,	w = 0.5 * w_start
      ),
      start_upper = c(
        d_0 = 1.5 * d_0_start,	E_D = 1.5 * E_D_start,
        DE = 1.5 * DE_start,	w = 1.5 * w_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, 0, 0, -Inf),
      upper = c(Inf, Inf, Inf, Inf)
    )
  )
  
  return(fit)
}

