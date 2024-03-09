functions {real switch_eta(real t, real t1, real eta, real mu, real xi) {
    return(eta+(1-eta)/(1+exp(xi*(t-t1-mu))));
}

  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real E = y[2];
      real I = y[3];
      real R = y[4];
      real D = y[5];
      
      real tswitch = x_r[1];
      real N = x_i[1];
      
  
      real beta = theta[1]; //Transmission rate
      real sigma = theta[2]; //Deaths Rate
      real nu = theta[3] ; // 1/Latent Period
      real delta = theta[4]; //Recovery Rate
      real eta = theta[5]; // reduction in transmission rate after quarantine
      real mu = theta[6]; // shift of quarantine implementation
      real xi = theta[7]; // slope of quarantine implementation

      real forcing_function = switch_eta(t,tswitch,eta,mu,xi); // switch function
      real beta_eff = beta * forcing_function; // beta decreased to take control measures into account
      
      real dS_dt = -beta_eff * I * S / N;
      real dE_dt =  beta_eff * I * S / N - nu * E;
      real dI_dt = nu * E - delta * I - sigma *I;
      real dR_dt =  delta * I;
      real dD_dt =  sigma * I ;
    
      
      return {dS_dt, dE_dt, dI_dt, dR_dt, dD_dt}; //define our SIRD model for stan, equation (2)
  }
}
data {
  int<lower=1> n_days; //number of observed days
  real y0[5]; //number of initial conditions 
  real t0; //initial time point
  real ts[n_days]; //time points observed
  int<lower=1> N; //population
  int<lower=0> newdeaths[n_days];
  real tswitch; //date that lockdown was introduced converted to a number
}
transformed data {
  real x_r[1] = {tswitch}; //real variables used to evaluate the function, which only depend on fixed data
  int x_i[1] = { N }; //integer values used to evaluate the function, which only depend on fixed data
}
parameters {
  real<lower=0> beta;
  real<lower=0> phi_inv; //Inverse of the Over Dispersion 
  real<lower=0> delta;
  real<lower=0> nu;
  real<lower=0> sigma;
  real<lower=0, upper=1> p_reported; // probability for an infected person to be reported (i.e counted as a case)
 real<lower=0,upper=1> eta; // reduction in transmission due to control measures (in proportion of beta)
  real<lower=0> mu; // shift of quarantine implementation (strictly positive as it can only occur after tswitch)
  real<lower=0,upper=1> xi_raw; // slope of quarantine implementation (strictly positive as the logistic must be downward)
 //our model paramaters and the inverse of our overdispersion paramater
}
transformed parameters{

  real y[n_days, 5];
  real phi = 1. / phi_inv;
real xi = xi_raw + 0.5;

  {
    real theta[7];
    theta[1] = beta;
    theta[2] = sigma;
    theta[3] = nu;
    theta[4]= delta;
    theta[5] = eta;
    theta[6] = mu;
    theta[7]=xi;


    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i); //solving the SEIRD ode with initial conditions
  } 
}

model {
//defining our prior distributions for our parmameters, informative wide priors so as not to skew the 
//posterior distribution, a narrower priors for italy sigma to allow for mixing of chains.
  //priors
 beta ~ lognormal(log(0.5),.25);
  delta ~ lognormal(log(1./7),.4);
  nu ~ lognormal(log(1.0/5.1),0.25);
  sigma ~ lognormal(log(0.015),.6);
   p_reported ~ beta(1, 2);
phi_inv ~ exponential(0.01); 
eta ~ beta(2.5, 4);
  mu ~ exponential(1./5);
  xi_raw ~ beta(1, 1);

//sampling distribution
//col(matrix x, int n) - The n-th column of matrix x. Here the number of dead people
newdeaths ~ neg_binomial_2(col(to_matrix(y), 5)*p_reported, phi); 
 
}
generated quantities {
  real R0 = beta / (delta + sigma); //equation used to calculate R0
  real Reff[n_days];
  real beta_switch[n_days];

  real recovery_time = 1 / delta; 
  real latent_period = 1 / nu;
  real pred_newdeaths[n_days];
 pred_newdeaths = neg_binomial_2_rng(col(to_matrix(y), 5)*p_reported, phi); // predicted deaths estimate for the posterior predictive check.
for (i in 1:n_days)
    Reff[i] = switch_eta(i, tswitch, eta, mu, xi) * beta / (delta+sigma);
    for (i in 1:n_days)
    beta_switch[i] = switch_eta(i, tswitch, eta, mu, xi) * beta;
  
}
