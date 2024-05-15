install.packages("ggplot2")
install.packages("tidyverse")
install.packages("minpack.lm")
install.packages('nls.multstart')
install.packages("rTPC")
install.packages("broom")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(minpack.lm)
library(nls.multstart)
library(rTPC)
library(broom)

setwd("C:/Users/mi620/OneDrive - Imperial College London/Year 4/FYP/R_stuff")
setwd("C:/Users/44785/OneDrive - Imperial College London/Year 4/FYP/R_stuff")

enzyme <- read.csv("Database_FYP.csv")

gdh <- enzyme %>% 
  filter(originalid == "GDH001") %>%
  select(interactor1temp, originaltraitvalue) 

colnames(gdh) <- c("temp", "rate")  #rate for rTPC test

##### rTPC package #####

mod = 'sharpeschoolfull_1981'
test_data = gdh

# get start vals
start_vals <- get_start_vals(test_data$temp, test_data$rate, model_name = 'sharpeschoolfull_1981')
start_vals

# get limits
low_lims <- get_lower_lims(test_data$temp, test_data$rate, model_name = 'sharpeschoolfull_1981')
upper_lims <- get_upper_lims(test_data$temp, test_data$rate, model_name = 'sharpeschoolfull_1981')
low_lims
upper_lims

#fit model
fit_rtpc <- nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
                     data = test_data,
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')

fit_rtpc

# Nonlinear regression model
# model: rate ~ sharpeschoolfull_1981(temp = temp, r_tref, e, el, tl,     eh, th, tref = 15)
# data: data
# r_tref         e        el        tl        eh        th 
# 4.568e+01 5.822e-01 9.473e-06 6.000e+01 1.151e+00 5.536e+01 
# residual sum-of-squares: 568.2
# 
# Number of iterations to convergence: 56 
# Achieved convergence tolerance: 1.49e-08


data_rtpc <- data.frame(temp = seq(min(test_data$temp), max(test_data$temp), 0.5))
preds_rtpc <- augment(fit_rtpc, newdata = data_rtpc)


gdh_rtpc <- ggplot(gdh, aes(temp, rate)) +
  geom_point(size=3) +
  geom_line(aes(temp, .fitted), preds_rtpc, col = '#4f4789', size=1)+
  ylab("Relative activity (%)") +
  xlab("Temperature (ºC)") +
  ggtitle("rTPC Model fittting (Sharpe-Schoolfield)")
gdh_rtpc



res=data.frame(residuals = residuals(fit))
res_hist <- ggplot(res, aes(residuals))+
  geom_histogram(fill="#fce762", color="black", bins=7)+
  ylab("Frequency")+
  xlab("Residuals")
res_hist
