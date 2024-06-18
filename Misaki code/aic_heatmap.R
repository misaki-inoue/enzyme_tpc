library(ggplot2)
library(dplyr)
library(tidyr)
library(geiger)
library(tibble)

setwd("C:/Users/mi620/OneDrive - Imperial College London/Year 4/FYP/R_stuff")
setwd("C:/Users/thinkpad/OneDrive - Imperial College London/Year 4/FYP/R_stuff")

stats <- read.csv("seven_models_stats_v1.csv")
stats_bic <- read.csv("seven_models_stats_BIC.csv")

plot.new()
par(mfrow=c(1,1))
# AIC distribution
hist(stats$AIC, main= "Distribution of AICc values")

stats$AIC <- as.numeric(stats$AIC)
stats$AIC[is.infinite(stats$AIC)] <- NA
stats$AIC[is.nan(stats$AIC)] <- NA

corrections <- 
  c("Sharpe_Schoolfield" = "Sharpe-Schoolfield",
    "Johnson_Lewin" =  "Johnson-Lewin",
    "Ratkowksy" = "Ross-Ratkowsky")
stats <- stats %>% mutate(model = str_replace_all(model, corrections))

ggplot(stats, aes(enzyme, model, fill= AIC)) + 
  geom_tile() +
  scale_fill_continuous(na.value = "grey90") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x=element_blank(),
        panel.background = element_blank())

#remove 6 enzyme TPC that were not able to be fitted for all models
remove_list <- c("AMA001", "GAC001", "GLB001", "HEP001", "SPS001", "XIS001")
cln_stats <- stats %>% filter(!enzyme %in% remove_list)
ggplot(cln_stats, aes(enzyme, model, fill= AIC)) + 
  geom_tile() +
  theme_grey(base_size = 18)+
  scale_fill_continuous(na.value = "grey80") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 12),
        axis.ticks.x=element_blank(),
        panel.background = element_blank()) +
  scale_fill_gradient(low = "red", high = "white")+
  labs(x="Enzyme ID",y="Model")

#Calculate and plot AIC weights
enzymes_with_na <- cln_stats %>%
  filter(is.na(Rsquared)) %>%
  select(enzyme) %>% pull(enzyme)
nona_stats <- cln_stats %>% filter(!enzyme %in% enzymes_with_na)

calculate_aicw <- function(df) {
  unique_enzyme_ids <- unique(df$enzyme)
  output_df <- data.frame()
  
  for (i in seq_along(unique_enzyme_ids)) {
    enzyme_id <- unique_enzyme_ids[i]
    AIC.scores <- subset(df, enzyme == enzyme_id) %>% .$AIC
    model_names <- subset(df, enzyme == enzyme_id) %>% .$model
    names(AIC.scores) <- model_names
    aicw_df <- aicw(AIC.scores) %>% tibble::rownames_to_column("model")
    aicw_df$enzyme <- enzyme_id
    
    output_df <- bind_rows(output_df, aicw_df)
  }
  return(output_df)
}

aic_weights <- calculate_aicw(nona_stats)

ggplot(aic_weights, aes(enzyme, model, fill= w)) + 
  geom_tile() + 
  theme_grey(base_size = 18)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.ticks.x=element_blank(),
        panel.background = element_blank())+
  labs(fill="AIC weights", x="Enzyme ID",y="Model") +
  scale_fill_gradient(low = "white", high = "red")

#Number of TPCs with AIC weights higher than threshold
aicw_50 <- aic_weights %>% filter(w > 0.50)
aicw_75 <- aic_weights %>% filter(w > 0.75)
aicw_80 <- aic_weights %>% filter(w > 0.80)
aicw_90 <- aic_weights %>% filter(w > 0.90)

aicw_50 %>% filter(model == "Johnson_Lewin") %>% nrow()
aic_weights %>%
  filter(model == "Sharpe_Schoolfield") %>% pull(w) %>% mean()

#R-squared vs. AIC values
ggplot(stats_bic, aes(x = AIC, y = BIC)) +
  geom_point() +
  facet_wrap(~ model, scales = "free") +
  labs(x = "AIC", y = "BIC")

# Bad model fitting for Johnson Lewin

azo <- data_to_fit %>%
  filter(originalid == "AZO001")

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
  
temp_pts <- data.frame(temp = seq(min(azo$temp), max(azo), 0.5))
preds <- augment(fit, newdata = temp_pts)
plot(adh$temp, log(adh$trait_value),
     type = "p",
     main = paste(df,model_name),
     xlab = "Temperature",
     ylab = "Activity")
lines(preds$temp, preds$.fitted, col = "red")

azo <- data_to_fit %>%
  filter(originalid == "AZO001") 

fit_and_plot(azo, fit_Johnson_Lewin_4_pars, "Johnson-Lewin")
