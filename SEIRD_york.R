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

covid_data_york <- read_excel("Data/covid_data_leeds_york.xlsx", 
                              sheet = "utla_2023-02-23_york", col_types = c("skip", 
                                                                            "skip", "date", "numeric", "numeric", 
                                                                            "numeric", "numeric"))
covid_data_york <- subset(covid_data_york,date>= "2020-03-17" & date <= "2020-7-1")

Deaths_cumulative_york <- covid_data_york$`Cumulative Deaths`
Deaths_york <- covid_data_york$`New Deaths`
Inf_cumulative_york <- covid_data_york$`Cumlative Cases`
Infected_york <- covid_data_york$`New Cases`

## compile model ---------------------------------------------------------------------------------------------------------------------------------------------

rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
cl <- parallel::makeCluster(3, setup_strategy = "sequential")
set.seed(3) # for reproductibility


N <- 211116 #total uk population

# times
n_days <- length(Deaths_cumulative_york)
t <- seq(1, n_days, by = 1)
t0 = 0
t <- t

#initial conditions
i0 <- Infected_york[1] #no of total cases on 14/03
e0 <- 1
s0 <- N - i0 - e0
r0 <- 0
d0 <- Deaths_cumulative_york[1] # no of total deaths on 14/03

y0 = c(S = s0,E = e0, I = i0, R = r0, D=d0)
date_switch <- "2020-03-23"
tswitch <- covid_data_york %>% filter(date < date_switch) %>% nrow() + 1
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, newdeaths = Deaths_cumulative_york,tswitch=tswitch)

model <- stan_model("SEIRD_leedsyork_switch.stan")

## sampling ---------------------------------------------------------------------------------------------------------------------------------------------

# test_york <- sampling(model,
#                 data = data_sir,
#                 iter = 1000,
#                 chains = 1,
#                 cores = 1,
#                 control = list(max_treedepth = 15))

seird_phase1_preported_york <- sampling(model,
                                        data = data_sir,
                                        iter = 3000,
                                        chains = 4,
                                        cores = 4,
                                        control = list(max_treedepth = 15, adapt_delta = 0.999))

## print output ---------------------------------------------------------------------------------------------------------------------------------------------

pars=c('beta','sigma','nu','delta', "R0", "recovery_time","phi",'latent_period','eta','mu','xi')
print(seird_phase1_preported_york, pars = pars, digits = 5)
#save(SEIRD_york_results, file = "SEIRD_york_results.RData") ## run only this line after R script has finished running

