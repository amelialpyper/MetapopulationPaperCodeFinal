functions {real switch_eta(real t, real t1, real eta, real mu, real xi) {
    return(eta+(1-eta)/(1+exp(xi*(t-t1-mu))));
}
    real switch_etaY(real t, real t1, real eta2, real mu2, real xi2) {
    return(eta2+(1-eta2)/(1+exp(xi2*(t-t1-mu2))));
}
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S1 = y[1];
      real E1 = y[2];
      real I1 = y[3];
      real R1 = y[4];
      real D1 = y[5];
      
      real S2 = y[6];
      real E2 = y[7];
      real I2 = y[8];
      real R2 = y[9];
      real D2 = y[10];
      
      real N1 = x_i[1];
      real N2 = x_i[2];
      
      real tswitch = x_r[1];
      real alphaP12 = x_r[2];
      real alphaP21 = x_r[3];
      real alphaA12 = x_r[4];
      real alphaA21 = x_r[5];

      real betaL = theta[1]; //transmission rate for city 1
      real betaY = theta[2]; //transmission rate for city 2
      real sigma = theta[3]; //death rate for city 1
      real sigmaY = theta[4]; //death rate for city 2
      real nu = theta[5] ; // 1/latent period for city 1
      real nuY = theta[6] ; // 1/latent period for city 2
      real gamma = theta[7]; //transport related infection rate for city 1
      real gammaY = theta[8]; //transport related infection rate for city 2
      real eta = theta[9]; // reduction in transmission rate after quarantine, city 1
      real mu = theta[10]; // shift of quarantine implementation, city 1
      real xi = theta[11]; // slope of quarantine implementation, city 1
      real eta2 = theta[12]; // reduction in transmission rate after quarantine, city 2
      real mu2 = theta[13]; // shift of quarantine implementation, city 2
      real xi2 = theta[14]; // slope of quarantine implementation, city 2
      real delta = theta[15]; //recovery rate for city 1
      real deltaY = theta[16]; //recovery rate for city 2 
      real p_reported = theta[17];
      real p_reportedY = theta[18];
      
      real forcing_function = switch_eta(t,tswitch,eta,mu,xi); // switch function
      real forcing_functionY = switch_eta(t,tswitch,eta2,mu2,xi2); 
      real betaL_eff = betaL * forcing_function; // beta decreased to take control measures into account
      real betaY_eff = betaY * forcing_functionY; 
      
      real dS1_dt = -betaL_eff * I1 * S1 / N1 - betaL_eff*S1*I2*(alphaA21)/N2 - betaY_eff*S1*I2*(alphaA12)/N2 - betaY_eff*S1*I1*(alphaA12)^2/N1 - gamma*S1*alphaP12*(alphaP12*I1/N1 + alphaP21*I2/N2);
      real dE1_dt =  betaL_eff * I1 * S1 / N1 + betaL_eff*S1*I2*(alphaA21)/N2 + betaY_eff*S1*I1*(alphaA12)^2/N1 + betaY_eff*S1*I2*(alphaA12)/N2  + gamma*S1*alphaP12*(alphaP12*I1/N1 + alphaP21*I2/N2) - nu * E1;
      real dI1_dt = nu * E1 - delta * I1- sigma * I1;
      real dR1_dt =  delta * I1;
      real dD1_dt =  sigma * I1 ;
      
      real dS2_dt = -betaY_eff * I2 * S2 / N2 - betaY_eff*S2*I1*(alphaA12)/N1 - betaL_eff*S2*I1*(alphaA21)/N1 - betaL_eff*S2*I2*(alphaA21)^2/N2 - gammaY*S2*alphaP21*(alphaP12*I1/N1 + alphaP21*I2/N2);
      real dE2_dt =  betaY_eff * I2 * S2 / N2 + betaY_eff*S2*I1*(alphaA12)/N1 + betaL_eff*S2*I1*(alphaA21)/N1 + betaL_eff*S2*I2*(alphaA21)^2/N2  + gammaY*S2*alphaP21*(alphaP12*I1/N1 + alphaP21*I2/N2) - nuY * E2;
      real dI2_dt = nuY * E2 - deltaY * I2 - sigmaY * I2;
      real dR2_dt =  deltaY * I2;
      real dD2_dt =  sigmaY * I2 ;
    
      
      return {dS1_dt, dE1_dt, dI1_dt, dR1_dt, dD1_dt,dS2_dt, dE2_dt, dI2_dt, dR2_dt, dD2_dt}; //define our SIRD model for stan, equation (2)
  }
}
data {
  int<lower=1> n_days; //number of observed days
  real y0[10]; //number of initial conditions 
  real t0; //initial time point
  real ts[n_days]; //time points observed
  int<lower=1> N1; //population
  int<lower=1> N2;
  int<lower=0> newdeaths1[n_days];
  int<lower=0> newdeaths2[n_days];
  real tswitch;
  real<lower=0> alphaA12;
  real<lower=0> alphaA21;
  real<lower=0> alphaP12;
  real<lower=0> alphaP21;
}
transformed data {
  real x_r[5] = {tswitch, alphaA12, alphaA21,alphaP12, alphaP21}; //real variables used to evaluate the function, which only depend on fixed data
  int x_i[2] = { N1, N2 }; //integer values used to evaluate the function, which only depend on fixed data
}
parameters {
  real<lower=0> betaL;
  real<lower=0> betaY;
  real<lower=0> phi_inv;
  real<lower=0> phiY_inv;
  real<lower=0> gamma;
  real <lower=0> gammaY;
  real<lower=0> nu;
  real<lower=0> nuY;
  real<lower=0> sigma;
  real<lower=0> sigmaY;
  real<lower=0> delta;
  real<lower=0> deltaY;
  real<lower=0, upper=1> p_reported;
  real<lower=0, upper=1> p_reportedY;// probability for an infected person to be reported (i.e counted as a case)
 real<lower=0,upper=1> eta; // reduction in transmission due to control measures (in proportion of beta)
  real<lower=0> mu; // shift of quarantine implementation (strictly positive as it can only occur after tswitch)
  real<lower=0,upper=1> xi_raw; // slope of quarantine implementation (strictly positive as the logistic must be downward)
   real<lower=0,upper=1> eta2; // reduction in transmission due to control measures (in proportion of beta)
  real<lower=0> mu2; // shift of quarantine implementation (strictly positive as it can only occur after tswitch)
  real<lower=0,upper=1> xi2_raw; // slope of quarantine implementation (strictly positive as the logistic must be downward)
 //our model paramaters and the inverse of our overdispersion paramater
}
transformed parameters{

  real y[n_days, 10];
  real phi = 1. / phi_inv;
  real phiY = 1. / phiY_inv;
real xi = xi_raw + 0.5;
real xi2 = xi2_raw + 0.5;

  {
    real theta[18];
    theta[1] = betaL;
    theta[2] = betaY;
    theta[3] = sigma;
    theta[4] = sigmaY;
    theta[5] = nu;
    theta[6] = nuY;
    theta[7]= gamma;
    theta[8]= gammaY;
    theta[9] = eta;
    theta[10] = mu;
    theta[11]=xi;
    theta[12] = eta2;
    theta[13] = mu2;
    theta[14]=xi2;
    theta[15] = delta;
    theta[16] = deltaY;
    theta[17] = p_reported;
    theta[18] = p_reportedY;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i); //solving the SEIRD ode with initial conditions
  } 
}

model {
//defining our prior distributions for our parmameters, informative wide priors so as not to skew the 
//posterior distribution, a narrower priors for italy sigma to allow for mixing of chains.
  //priors
 betaL ~ lognormal(log(0.5),.25);
 betaY ~ lognormal(log(0.5),.25);
 gamma ~ lognormal(log(0.5),.25);
 gammaY ~ lognormal(log(0.5),.25);
  delta ~ lognormal(log(1./7),.4);
  deltaY ~ lognormal(log(1./7), 0.4);
  nu ~ lognormal(log(1.0/5.1),0.25);
  nuY ~ lognormal(log(1.0/5.1),0.25);
  sigma ~ lognormal(log(0.015),.6);
   sigmaY ~ lognormal(log(0.015),.6);
   p_reported ~ beta(1, 2);
   p_reportedY ~ beta(1, 2);
phi_inv ~ exponential(0.01); 
phiY_inv ~ exponential(0.01); 
eta ~ beta(2.5, 4);
  mu ~ exponential(1./5);
  xi_raw ~ beta(1, 1);
  eta2 ~ beta(2.5, 4);
  mu2 ~ exponential(1./5);
  xi2_raw ~ beta(1, 1);
//sampling distribution
//col(matrix x, int n) - The n-th column of matrix x. Here the number of dead people
newdeaths1 ~ neg_binomial_2(col(to_matrix(y), 5)*p_reported, phi); //equation (6)
newdeaths2 ~ neg_binomial_2(col(to_matrix(y), 10)*p_reportedY, phiY);
 
}
generated quantities {
  real beta_switchL[n_days];
  real beta_switchY[n_days];
  
  real recovery_time = 1 / delta;
  real recovery_timeY = 1 / deltaY;
  real latent_period = 1 / nu;
  real latent_periodY = 1 / nuY;
  real pred_newdeaths1[n_days];
  real pred_newdeaths2[n_days];
 pred_newdeaths1 = neg_binomial_2_rng(col(to_matrix(y), 5)*p_reported, phi); // predicted deaths estimate for the posterior predictive check.
 pred_newdeaths2 = neg_binomial_2_rng(col(to_matrix(y), 10)*p_reportedY, phiY);
    for (i in 1:n_days)
    beta_switchY[i] = switch_eta(i, tswitch, eta2, mu2, xi2) * betaY;
       for (i in 1:n_days)
    beta_switchL[i] = switch_eta(i, tswitch, eta, mu, xi) * betaL;
}
