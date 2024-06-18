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

enzyme <- read.csv("BacterialData.csv")

# Load functions all 7 functions from "mech_model_functions" before running multi_fit()

multi_fit <- function(df, model_func, model_name) {
  plot.new()
  subplot_counter = 0
  stats <- data.frame(model = character(), 
                      enzyme = character(), 
                      Rsquared = numeric(), 
                      AIC = numeric(),
                      BIC = numeric())
  #AICc = numeric(), BIC = numeric()
  
  enz_id <- df %>%
    distinct(originalid)
  
  tpc_to_fit <- df %>%
    filter(originalid %in% enz_id$originalid) %>%
    select(originalid, interactor1temp, originaltraitvalue)
  
  colnames(tpc_to_fit) <- c('originalid','temp','trait_value')
  
  unique_enzyme_ids <- unique(tpc_to_fit$originalid)
  num_enzymes <- length(unique_enzyme_ids)
  
  par(mfrow = c(6, 5), mar = c(2, 2, 1, 1), oma=c(0, 0, 2, 0))
  
  
  for (i in 1:num_enzymes) {
    enzyme_id <- unique_enzyme_ids[i]
    enzyme_data <- subset(tpc_to_fit, originalid == enzyme_id) %>%
      select(temp, trait_value)
    
    fit <- NULL
    
    #Skip fitting if fitting function produce some sort of error
    tryCatch({
      fit <- model_func(enzyme_data)
    }, error = function(e) {
      warning(paste("Error fitting enzyme", enzyme_id, ":", e$message))
    })
    
    if (!is.null(fit)) {
      subplot_counter <- subplot_counter + 1
      
      aic_val <- AIC(fit)
      bic_val <- BIC(fit)
      
      stats[i,] <- c(model_name, enzyme_id, r_squared, aic_val, bic_val) #aicc_val, bic_val
      print(paste("Fit for", enzyme_id, "is completed.", i, "out of", num_enzymes))
      
      temp_pts <- data.frame(temp = seq(min(enzyme_data$temp), max(enzyme_data$temp), 0.5))
      preds <- augment(fit, newdata = temp_pts)
      plot(enzyme_data$temp, log(enzyme_data$trait_value),
           type = "p",
           main = enzyme_id,
           xlab = "Temperature",
           ylab = "Activity")
      lines(preds$temp, preds$.fitted, col = "red")
      if (subplot_counter == 1 || subplot_counter %% 30 == 1) {
        title(main = model_name, outer = TRUE, cex.main = 2)
      }
    }
    else {
      subplot_counter <- subplot_counter + 1
      print(paste("Fit for", enzyme_id, "is NULL.", i, "out of", num_enzymes))
      stats[i,] <- c(model_name ,enzyme_id, NA, NA, NA)
      plot(enzyme_data$temp, log(enzyme_data$trait_value),
           type = "p",
           main = enzyme_id,
           xlab = "Temperature",
           ylab = "Activity")
      if (subplot_counter == 1 || subplot_counter %% 30 == 1) {
        title(main = model_name, outer = TRUE, cex.main = 2)
      }
    }
  }
  return(stats)
}

hinshel <- multi_fit(df=data_to_fit, fit_Hinshelwood_4_pars, "Hinshelwood")
johnson <- multi_fit(df=data_to_fit, fit_Johnson_Lewin_4_pars, "Johnson-Lewin")
sharpe <- multi_fit(df=data_to_fit, fit_Sharpe_Schoolfield_6_pars, "Sharpe-Schoolfield")
ratkowsky <- multi_fit(data_to_fit, fit_Ross_Ratkowsky_5_pars, "Ratkowsky")
hobbs <- multi_fit(data_to_fit, fit_Hobbs_4_pars, "Hobbs")
eaar <- multi_fit(data_to_fit, fit_Enzyme_Assisted_Arrhenius_5_pars, "EAAR")
ritchie <- multi_fit(data_to_fit, fit_Ritchie_4_pars,"Ritchie")

stats <- bind_rows(vant, hinshel, johnson, sharpe, ratkowsky, hobbs, eaar, ritchie)
write.csv(stats, "seven_models_stats_v1.csv", row.names = FALSE)

