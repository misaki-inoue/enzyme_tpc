library(ggplot2)
library(dplyr)
library(minpack.lm)
library(nls.multstart)
library(rTPC)
library(broom)
library(MuMIn)
library(geiger)
library(tidyr)

setwd("C:/Users/mi620/OneDrive - Imperial College London/Year 4/FYP/R_stuff")
setwd("C:/Users/44785/OneDrive - Imperial College London/Year 4/FYP/R_stuff")

enzyme <- read.csv("Database_FYP.csv")

first_thirty <- enzyme %>% # Only first thirty enzymes are selected to reduce computing time
  distinct(originalid) %>%
  head(30)

thirty_enz <- enzyme %>%
  filter(originalid %in% first_thirty$originalid) %>%
  select(originalid, interactor1temp, originaltraitvalue)

# Load functions all 8 functions from "mech_model_functions" before running multi_fit()
# Note to self: Need to add R-squared to the output dataframe

multi_fit <- function(df, model_func, model_name, num_col) {
  plot.new()
  stats <- data.frame(model = character(), enzyme = character(), AICc = numeric(), BIC = numeric())
  
  enz_id <- df %>%
    distinct(originalid)
  
  tpc_to_fit <- df %>%
    filter(originalid %in% enz_id$originalid) %>%
    select(originalid, interactor1temp, originaltraitvalue)
  
  colnames(tpc_to_fit) <- c('originalid','temp','trait_value')
  
  unique_enzyme_ids <- unique(tpc_to_fit$originalid)
  num_enzymes <- length(unique_enzyme_ids)
  num_row <- ceiling(num_enzymes/num_col)
  
  par(mfrow = c(num_row, num_col), mar = c(2, 2, 2, 1))
  
  for (i in 1:num_enzymes) {
    enzyme_id <- unique_enzyme_ids[i]
    enzyme_data <- subset(tpc_to_fit, originalid == enzyme_id) %>%
      select(temp, trait_value)
    
    fit <- NULL
    fit <- model_func(enzyme_data)
    
    if (!is.null(fit)) {
      aicc_val <- AICc(fit)
      
      loglik <- logLik(fit)
      n   <- attributes(loglik)$nobs 
      p   <- attributes(loglik)$df
      dev <- -2*as.numeric(loglik)
      bic_val  <- dev +  p * log(n)
      
      stats[i,] <- c(model_name, enzyme_id, aicc_val, bic_val)
      
      temp_pts <- data.frame(temp = seq(min(enzyme_data$temp), max(enzyme_data$temp), 0.5))
      preds <- augment(fit, newdata = temp_pts)
      
      plot(enzyme_data$temp, log(enzyme_data$trait_value), 
           type = "p", 
           main = enzyme_id, 
           xlab = "Temperature", 
           ylab = "Activity")
      lines(preds$temp, preds$.fitted, col = "red")
    } else {
      warning(paste("Fit for enzyme", enzyme_id, "is NULL"))
      stats[i,] <- c(model_name ,enzyme_id, NA, NA)
      plot(enzyme_data$temp, log(enzyme_data$trait_value), 
           type = "p", 
           main = enzyme_id, 
           xlab = "Temperature", 
           ylab = "Activity")
    }
    print(enzyme_id)
  }
  return(stats)
} 


vant <- multi_fit(df=thirty_enz, fit_Vant_Hoff_4_pars, "Vant_Hoff",5)
hinshel <- multi_fit(df=thirty_enz, fit_Hinshelwood_4_pars, "Hinshelwood", 5)
johnson <- multi_fit(df=thirty_enz, fit_Johnson_Lewin_4_pars, "Johnson_Lewin", 5)
shape <- multi_fit(df=thirty_enz, fit_Sharpe_Schoolfield_6_pars, "Sharpe_Schoolfield", 5)
ratkowsky <- multi_fit(thirty_enz, fit_Ross_Ratkowsky_5_pars, "Ratkowksy", 5)
hobbs <- multi_fit(thirty_enz, fit_Hobbs_4_pars, "Hobbs", 5)
eaar <- multi_fit(thirty_enz, fit_Enzyme_Assisted_Arrhenius_5_pars, "EAAR", 5)
ritchie <- multi_fit(thirty_enz, fit_Ritchie_4_pars,"Ritchie", 5)
