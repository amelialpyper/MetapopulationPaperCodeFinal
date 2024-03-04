library(tidyverse)
library(readr)
library(readxl)
library(ggplot2)
library(lubridate)
library(rstan)
library(gridExtra)
library(parallel)
library(deSolve)
library(dplyr)
library(rstanarm)
library(writexl)

## load data ---------------------------------------------------------------------------------------------------------------------------------------------

covid_data_leeds <- read_excel("Data/covid_data_leeds_york.xlsx",
                               sheet = "utla_2023-02-23_leeds", col_types = c("skip",
                                                                              "skip", "date", "numeric", "numeric",
                                                                              "numeric", "numeric"))

covid_data_leeds <- subset(covid_data_leeds,date>= "2020-03-17" & date <= "2020-7-1") #subset the data for the date range we require

#assign each column a label we will use going forward
Deaths_cumulative <- covid_data_leeds$cumDeaths28DaysByDeathDate
Deaths <- covid_data_leeds$newDeaths28DaysByDeathDate
Inf_cumulative <- covid_data_leeds$cumCasesBySpecimenDate
Infected <- covid_data_leeds$newCasesBySpecimenDate

## compile model ---------------------------------------------------------------------------------------------------------------------------------------------
#setting up stan options
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
cl <- parallel::makeCluster(3, setup_strategy = "sequential")
set.seed(3) # for reproductibility


N <- 795430  #total leeds population 2020, based on population projection value from Leeds Observatory 

# times
n_days <- length(Deaths_cumulative)
t <- seq(1, n_days, by = 1)
t0 = 0
t <- t

#initial conditions
i0 <- Infected[1] #no of total cases on 17/03
e0 <- 1
s0 <- N - i0 - e0
r0 <- 0
d0 <- Deaths_cumulative[1] # no of total deaths on 17/03

y0 = c(S = s0,E = e0, I = i0, R = r0, D=d0)
date_switch <- "2020-03-23" #date that full lockdown measure came into force
tswitch <- covid_data_leeds %>% filter(date < date_switch) %>% nrow() + 1 #Turn this date into a numerical value
# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, newdeaths = Deaths_cumulative,tswitch=tswitch)

model <- stan_model("SEIRD_leedsyork_switch.stan")

## sampling ---------------------------------------------------------------------------------------------------------------------------------------------

# test_leeds <- sampling(model,
# data = data_sir,
# iter = 1000,
# chains = 1,
# cores = 1,
# control = list(max_treedepth = 15))

seird_phase1_preported <- sampling(model,
                                   data = data_sir,
                                   iter = 3000,
                                   chains = 4,
                                   cores = 4,
                                   control = list(max_treedepth = 15))

#show the parameter estimates

## print output ---------------------------------------------------------------------------------------------------------------------------------------------

pars=c('beta','sigma','nu','delta', "R0", "recovery_time","phi",'latent_period','eta','mu','xi')
print(seird_phase1_preported, pars = pars, digits = 5)

#save(SEIRD_leeds_results, file = "SEIRD_leeds_results_actual.RData") ## run only this line after R script has finished running
