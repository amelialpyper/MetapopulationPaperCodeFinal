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
covid_data_york <- read_excel("Data/covid_data_leeds_york.xlsx", 
                              sheet = "utla_2023-02-23_york", col_types = c("skip", 
                                                                            "skip", "date", "numeric", "numeric", 
                                                                            "numeric", "numeric"))

covid_data_york <- subset(covid_data_york,date>= "2020-03-17" & date <= "2020-7-1")

Deaths_cum_york <- covid_data_york$`Cumulative Deaths`
Deaths_york <- covid_data_york$`New Deaths`
Inf_cum_york <- covid_data_york$`Cumlative Cases`
Infected_york <- covid_data_york$`New Cases`


rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
cl <- parallel::makeCluster(4, setup_strategy = "sequential")

set.seed(3) # for reproductibility

N <- 211116 #total uk population

# times
n_days <- length(Deaths_cum_york)
t <- seq(1, n_days, by = 1)
t0 = 0
t <- t

#initial conditions
i0 <- Infected_york[1] #no of cases on 17/03
e0 <- 1
s0 <- N - i0 - e0
r0 <- 0
d0 <- Deaths_cum_york[1] # no of total deaths on 17/03

y0 = c(S = s0,E = e0, I = i0, R = r0, D=d0)
date_switch <- "2020-03-23"
tswitch <- covid_data_york %>% filter(date < date_switch) %>% nrow() + 1
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, newdeaths = Deaths_cum_york,tswitch=tswitch)


model <- stan_model("seirddeaths_switch.stan")


# test_york <- sampling(model,
#                 data = data_sir,
#                 seed = 3,
#                 iter = 1000,
#                 chains = 1,
#                 cores = 1,
#                 control = list(max_treedepth = 15))



seird_phase1_preported_york <- sampling(model,
                                        data = data_sir,
                                        seed = 3,
                                        iter = 3000,
                                        chains = 4,
                                        cores = 4,
                                        control = list(max_treedepth = 15, adapt_delta = 0.999))

pars=c('beta','sigma','nu','delta', "R0", "recovery_time","phi",'latent_period','eta','mu','xi')
print(seird_phase1_preported_york, pars = pars, digits = 5)

#save(seird_york_paper_results, file = "seird_york_paper_results_actual.RData") ## Run only this line after R script has completed to save results to file
