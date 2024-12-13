data {
  int<lower=1> N; //number of subjects
  int<lower=1> P; //number of covariates
  int<lower=1> T; //number of total observations = dates(first patient) + dates(second) + ...
  int<lower=1> subj[T]; //subject id vector
  //int<lower=1> counter[N]; //exact number of observations for every patient
  real y[T]; //outcome
  matrix[T,P] X; //predictors
}

parameters {
  vector[P] beta; //fixed intercept and slope
  vector[N] b; //subject intercepts (note: fixed for each patient, not time-dependent)
  real<lower=0> eta; //sd for subject intercepts
  real<lower=0> sigma_e; //error sd
  real mub;
  
}

model {
  //priors
  mub ~ normal(0, 2); //subj random effects
  eta ~ inv_gamma(3,2); //variance of intercepts
  b ~ normal(mub, eta); //subj random effects
  beta ~ normal(0,5);
  sigma_e ~ inv_gamma(3,2);
  
  vector[T] mu;
  mu = X * beta + b[subj];
  
  y ~ normal(mu, sigma_e);
}

generated quantities {

  vector[T] y_hat; 

   y_hat = X * beta + b[subj];

}
