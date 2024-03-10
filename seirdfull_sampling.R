## Import Data---------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(readr)
library(readxl)
library(ggplot2)
library(rstan)
library(gridExtra)
library(parallel)
library(deSolve)
library(dplyr)
library(writexl)

covid_data_leeds <- read_excel("Data/covid_data_leeds_york.xlsx", 
                               sheet = "utla_2023-02-23_leeds", col_types = c("skip", 
                                                                              "skip", "date", "numeric", "numeric", 
                                                                              "numeric", "numeric"))
covid_data_york <- read_excel("Data/covid_data_leeds_york.xlsx", 
                              sheet = "utla_2023-02-23_york", col_types = c("skip", 
                                                                            "skip", "date", "numeric", "numeric", 
                                                                            "numeric", "numeric"))

covid_data_leeds <- subset(covid_data_leeds,date>= "2020-03-17" & date <= "2020-7-1")
covid_data_york <- subset(covid_data_york,date>= "2020-03-17" & date <= "2020-7-1")


Deaths_cum_leeds <- covid_data_leeds$cumDeaths28DaysByDeathDate
Deaths_leeds <- covid_data_leeds$newDeaths28DaysByDeathDate
Inf_cum_leeds <- covid_data_leeds$cumCasesBySpecimenDate
Infected_leeds <- covid_data_leeds$newCasesBySpecimenDate 

Deaths_cum_york <- covid_data_york$`Cumulative Deaths`
Deaths_york <- covid_data_york$`New Deaths`
Inf_cum_york <- covid_data_york$`Cumlative Cases`
Infected_york <- covid_data_york$`New Cases`

rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
cl <- parallel::makeCluster(3, setup_strategy = "sequential")
set.seed(3) # for reproductibility


N1 <- 795430 #total leeds population https://worldpopulationreview.com/world-cities/leeds-population
N2 <- 211116 #york estimate for 2020, 2019+2023/2 + 2019 divide by 2 again

# times
n_days <- length(Deaths_cum_leeds)
t <- seq(1, n_days, by = 1)
t0 = 0
t <- t

#initial conditions
i01 <- Infected_leeds[1] #no of total cases on 14/03
e01 <- 1
s01 <- N1 - i01 - e01
r01 <- 0
d01 <- Deaths_cum_leeds[1] # no of total deaths on 14/03
i02 <- Infected_york[1] #no of total cases on 14/03
e02 <- 1
s02 <- N2 - i02 - e02
r02 <- 0
d02 <- Deaths_cum_york[1] # no of total deaths on 14/03

y0 = c(S1 = s01,E1 = e01, I1 = i01, R1 = r01, D1=d01,S2 = s02,E2 = e02, I2 = i02, R2 = r02, D=d02)
date_switch <- "2020-03-23" # date of introduction of control measures
tswitch <- covid_data_leeds %>% filter(date < date_switch) %>% nrow() + 1 # convert time to number

# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N1 = N1, N2 = N2, newdeaths1 = Deaths_cum_leeds, newdeaths2 = Deaths_cum_york,
                 tswitch=tswitch, alphaP12 = 0.0006, alphaP21 = 0.0056,alphaA12 = 0.0034, alphaA21 = 0.0253)


model <- stan_model("~/Library/CloudStorage/OneDrive-UniversityofLeeds/PhD Prep/Thesis/R/Bayesian Inference/seird_full_switch_actual.stan")

# test_full <- sampling(model,
#                 data = data_sir,
#                 iter = 2000,
#                 chains = 1,
#                 cores = 1,
# control = list(max_treedepth = 15))

fullmeta <- sampling(model,
                            data = data_sir,
                            seed = 3,
                            iter = 6000,
                            chains = 2,
                            cores = 2,
                            control = list(max_treedepth = 15))

pars=c('betaL','betaY','sigma','sigmaY','delta','deltaY','gamma','gammaY','nu','nuY', "recovery_time","recovery_timeY","phi",'phiY','latent_period','latent_periodY','eta','mu','xi','eta2','mu2','xi2')
print(fullmeta, pars = pars, digits = 4)

#fullmeta_params <- as.data.frame(get_posterior_mean(seirdfull_paper_results$fullmeta,pars)[,3])
#write_xlsx(fullmeta_params, "Parameters/seird_full_paper.xlsx", col_names = FALSE)

#save(seirdfull_paper_results, file = "seirdfull_paper_results_actual.RData") ## Run only this line after R script has completed to save results to file
