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

cbPalette <- c("#E69F00", "#56B4E9","#F0E442", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")
ggplot(data = pred_values, aes(temp, .fitted)) +
  geom_line(aes(col=Models, linewidth=linetype)) +
  scale_linetype_manual(values = c(1,2)) +
  scale_linewidth_manual(values = c(2.5,0.7), guide="none") +
  scale_colour_manual(values=cbPalette) +
  labs(x = "Temperature", y = "log(% Enzyme Activity)",
       title = paste("Enzyme ID:", enzyme_id)) +
  geom_point(data = azo, aes(temp, log(trait_value)))+
  theme_bw(base_size = 12)

fit <- fit_Ross_Ratkowsky_5_pars(adh)
temp_pts <- data.frame(temp = seq(min(adh$temp), max(adh$temp), 0.5))
preds <- augment(fit, newdata = temp_pts)
plot(adh$temp, log(adh$trait_value),
     type = "p",
     main = "ADH Sharpe",
     xlab = "Temperature",
     ylab = "Activity")
lines(preds$temp, preds$.fitted, col = "red")
ggplot(data = preds, aes(temp, .fitted)) +
  geom_line()

       