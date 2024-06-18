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
setwd("C:/Users/thinkpad/OneDrive - Imperial College London/Year 4/FYP/R_stuff")

enzyme <- read.csv("BacterialData.csv")

gpp <- enzyme %>%
  filter(originalid == "GPP001") %>%
  select(interactor1temp, originaltraitvalue) 
colnames(gpp) <- c('temp','trait_value')
azo <- enzyme %>%
  filter(originalid == "AZO001") %>%
  select(interactor1temp, originaltraitvalue) 
colnames(azo) <- c('temp','trait_value')
seven_models <- list(fit_Hinshelwood_4_pars, 
                  fit_Johnson_Lewin_4_pars,
                  fit_Sharpe_Schoolfield_6_pars,
                  fit_Ross_Ratkowsky_5_pars,
                  fit_Hobbs_4_pars,
                  fit_Enzyme_Assisted_Arrhenius_5_pars,
                  fit_Ritchie_4_pars)
seven_models_names <- c("Hinshelwood", 
                     "Johnson-Lewin",
                     "Sharpe-Schoolfield",
                     "Ross-Ratkowsky",
                     "Hobbs",
                     "Enzyme-Assisted-Arrhenius",
                     "Ritchie")

make_fitted_plot <- function(df, enzyme_id, model_list, name_list) {
  pred_values <- data.frame(Models = as.character(),
                            temp = as.numeric(),
                            .fitted = as.numeric())
  for (i in seq_along(model_list)) {
    model_fun <- model_list[[i]]
    fit <- model_fun(df)
    temp_pts <- data.frame(temp = seq(min(df$temp), max(df$temp), 0.5))
    preds <- augment(fit, newdata = temp_pts)
    preds$Models <- name_list[i]
    pred_values <- bind_rows(pred_values, preds) 
  }
  pred_values <- pred_values %>%
    mutate(linetype = if_else(Models == "Johnson-Lewin", "best", "other"))
  cbPalette <- c("#E69F00", "#56B4E9", "#F0E442", "#009E73", "#CC79A7","#0072B2", "#D55E00")
  ggplot(data = pred_values, aes(temp, .fitted)) +
    geom_line(aes(col=Models, linewidth=linetype, linetype=linetype)) +
    scale_linetype_manual(values = c(1,5),guide="none") +
    scale_linewidth_manual(values = c(1.5,0.8), guide="none") +
    #scale_colour_manual(values=cbPalette)+
    scale_color_brewer(palette="Paired")+
    geom_point(data = df, aes(temp, log(trait_value)))+
    labs(x = "Temperature (°C)", 
         y = "log(% Enzyme Activity)")+
    theme_bw(base_size = 18)+
    theme(legend.position = "bottom")
}

make_fitted_plot(azo, "AZO001", seven_models, seven_models_names)
make_fitted_plot(gpp, "GPP001", seven_models, seven_models_names)

# Use function below to plot a single model onto single TPC
fit_and_plot <- function(df, model_func, model_name) {
  enz_id <- df %>%
    distinct(originalid)
  
  tpc_to_fit <- df %>%
    filter(originalid %in% enz_id$originalid) %>%
    select(originalid, interactor1temp, originaltraitvalue)
  
  colnames(tpc_to_fit) <- c('originalid','temp','trait_value')
  
  unique_enzyme_ids <- unique(tpc_to_fit$originalid)
  num_enzymes <- length(unique_enzyme_ids)
  
  for (i in 1:num_enzymes) {
    enzyme_id <- unique_enzyme_ids[i]
    enzyme_data <- subset(tpc_to_fit, originalid == enzyme_id) %>%
      select(temp, trait_value)
    fit <- model_func(enzyme_data)
    
    temp_pts <- data.frame(temp = seq(min(enzyme_data$temp), max(enzyme_data$temp), 0.5))
    preds <- augment(fit, newdata = temp_pts)
    plot(enzyme_data$temp, log(enzyme_data$trait_value),
         type = "p",
         main = paste(enzyme_id, model_name),
         xlab = "Temperature",
         ylab = "Activity")
    lines(preds$temp, preds$.fitted, col = "red")
  }
}

fit_and_plot(azo, fit_Johnson_Lewin_4_pars, "Johnson-Lewin")

       