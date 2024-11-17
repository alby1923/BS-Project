data {
  int<lower=1> N; //number of subjects
  int<lower=1> P; //number of covariates
  int<lower=1> T; //number of total observations = dates(first patient) + dates(second) + ...
  int<lower=1> subj[T]; //subject id vector
  int y[T]; //outcome
  matrix[T,P] X; //predictors
}

parameters {
  vector[P] beta; //fixed intercept and slope
  vector[N] b; //subject intercepts (note: fixed for each patient, not time-dependent)
  real<lower=0> eta; //sd for subject intercepts
  real<lower=0> sigma_e; //error sd
  // had to add vector[N] mub;
}

model {
  vector[T] mu;
  //priors
  
  mub ~ normal(0, 2); //subj random effects
  eta ~ inv_gamma(3,2); //variance of intercepts
  b ~ normal(mub, eta); //subj random effects
  //???w ~ normal(0, sigma_w); //item random effects
  beta ~ normal(0,5);
  sigma_e ~ normal(0,5);
  
  //likelihood
  for (i in 1:T){
    mu = x[i]*beta + b[subj[i]]; //same as below
  }
  y ~ normal(mu, sigma_e);
}
generated quantities {
  vector[T] log_lik; //log-likelihood
  for (i in 1:T) {
    //generate predicted value
    real y_hat = x[i]*beta + b[subj[i]]; //x or X ???? Also, with X I get this error
    //Ill-typed arguments supplied to assignment operator =: lhs has type vector and rhs has type real
    //calculate log-likelihood
    log_lik[i] = normal_lpdf(y[i] | y_hat, e);
    //normal_lpdf is the log of the normal probability density function
  }
}
